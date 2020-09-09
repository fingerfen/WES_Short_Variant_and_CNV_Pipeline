## Import necessary packages
import numpy as np
import pandas as pd
import sys
import scipy.stats as stats
import sqlite3
import argparse


## Define arg parser to take in inputs
parser = argparse.ArgumentParser(description="""Individual Sample Analysis: Using GO terms and 
                                                sample's SNP data to perform metaP Fisher test""") 


parser.add_argument('--database', required=True,
                    help='path to the sqlite3 database file that holds the sample names of the two cohorts')
parser.add_argument('--goterms', required=True,
                    help='path to the file containing go term and associated gene\
                    1st column is gene name, 2nd column is go_term')
parser.add_argument('--gofunctions', required=True,
                    help='path to the file containing go term and its functions')


args = parser.parse_args()

## Read in files
human_GO = pd.read_csv(args.goterms, delimiter='\t')
go_def = pd.read_csv(args.gofunctions, delimiter='\t')


############  Generate a list of UNIQUE go terms #################
u_list_GO = pd.Series(np.unique(human_GO['GO']))
u_list_GO = u_list_GO.loc[u_list_GO != 'all']
u_list_GO_gene = pd.Series(u_list_GO.copy())


for u_GO_index,u_GO in enumerate(u_list_GO):
    bob = list((human_GO.loc[human_GO['GO'] == u_GO])['SYMBOL'])
    u_list_GO_gene.at[u_GO_index] = bob
    print(u_GO_index)

##################################################################

conn = sqlite3.connect(args.database)

## Read in the "id_table" that holds the list of all sample IDs
df = pd.read_sql_query("SELECT * from id_table", conn)
matrix_temp = pd.read_sql_query("SELECT * from short_variant_mutation_matrix", conn)
#Make a copy of the gene column
gene_assoc_vcf = matrix_temp['GENE'].copy()
#Then drop the gene column from the table so the table returns to normal GATK vcf format
matrix_temp = matrix_temp.drop(columns='GENE')

conn.close()


mutation_info = matrix_temp.iloc[:,:9]
mutation_matrix = matrix_temp.iloc[:,9:]

filtered_col_names = pd.DataFrame(mutation_matrix.columns, columns=['ID'], dtype='int64')

## Get the table with all of the ID names
merged_id = df.merge(filtered_col_names, on='ID' )
id_column = merged_id['ID'].astype('str')


print('Done importing files')
##Evaluation of individual SNP GO


## Create 2 lists. One to hold the pval_id for each sample
## The other one is to hold the oddratio_id for each sample
## This is made so that these lists could be used to name the columns of the resulting matrix
pval_id = [id + '_pval' for id in id_column.unique()]
oddratio_id = [id + '_or' for id in id_column.unique()]

########################################################################################
## The code below basically try to set up a contingency table to perform the fisher test
## The contingency table needs 4 numbers, thus the LEN() function appears 4 times below
## After the 4 len() functions are called, the contigency table is done and the fisher test get called
########################################################################################


## Create an empty dataframe to store the results
short_variant_GO = pd.DataFrame(0, index=u_list_GO, columns = list(sum(zip(pval_id,oddratio_id),())))

## Get the number of genes in the human
gene_in_human = len(human_GO['SYMBOL'].unique())
sample_index = 0

for sample in id_column:
    print(sample)
    sample_mutation_matrix = matrix_temp.loc[matrix_temp[str(sample)] != 0]
    ########## CHANGE BED TO HERE
    sample_gene_assoc = gene_assoc_vcf.loc[matrix_temp[str(sample)] != 0]
    sample_gene_list = [] 
    sample_gene_assoc.map(lambda x: sample_gene_list.extend(x.split(',')))
    u_sample_gene_assoc = np.unique([i for i in sample_gene_list if i])
    u_sample_gene_assoc_set = set(u_sample_gene_assoc)
    gene_in_sample = len(u_sample_gene_assoc)
    
    go_gene_table = pd.concat([u_list_GO,u_list_GO_gene],axis=1)
    go_gene_table.columns = ['GO','gene']
    
    go_logical = go_gene_table['gene'].map(lambda x: True if u_sample_gene_assoc_set.intersection(set(x)) else False)
    go_term_sample = go_gene_table.loc[go_logical]["GO"]
    sample_index += 1
    print('This is sample number {}'.format(sample_index))
    term_index = 0
    for go_term in go_term_sample:
        term_index += 1
        print(term_index)
        go_term_gene = (go_gene_table.loc[go_gene_table['GO'] == go_term])['gene'].iloc[0]
        ## Getting the last two entries for the contingency table
        gene_in_GO_human = len(go_term_gene)
        gene_in_GO_sample = len(set(u_sample_gene_assoc).intersection(set(go_term_gene)))
        ##CAlling the fisher test
        short_variant_GO.loc[[go_term],[sample + '_or']], short_variant_GO.loc[[go_term],[sample + '_pval']] = stats.fisher_exact([[gene_in_GO_human, gene_in_human], [gene_in_GO_sample, gene_in_sample]])
        
## Put together the column names
new_list = list(sum(zip(pval_id,oddratio_id),()))

##Name the columns
short_variant_GO.columns = new_list

## Merge to annotate the results with the GO functions
short_variant_GO_merged = short_variant_GO.merge(go_def,left_index=True,right_on='GO')

## Make connection to the SQLite3 DataBase

conn = sqlite3.connect(args.database)

## Insert the computed dataframe as a table in the sqlite3 database
short_variant_GO_merged.to_sql('short_variant_sample_based_GO', conn, if_exists="replace", index=False)

## close the connection to save data
conn.close()

##Output a file so Snakemake can determine termination of the script when it's done running
short_variant_GO_merged.to_csv(sys.stdout, sep='\t', index=False)