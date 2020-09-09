import numpy as np
import pandas as pd
import sys
import scipy.stats as stats
import sqlite3
import argparse

## Define a function that get used often
## Adapted online 'https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744'
def read_vcf(path):
    import pandas as pd
    import io
    ## Obtain the lines without ## as the beginning
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    ## Obtain the lines with ## as the beginning - This gets export as header later
    with open(path, 'r') as f:
        header_lines = [l for l in f if l.startswith('##')]
    return header_lines, pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t')


## Define arg parser to take in inputs
parser = argparse.ArgumentParser(description="""MP term analysis: This script
                                                performs the GO term analysis
                                                by taking a look at mutations 
                                                called for each MP term and 
                                                compare them between the
                                                Experimental vs Control groups""") 

parser.add_argument('--database', required=True,
                    help='path to the sqlite3 database file that holds the sample names of the two cohorts')
parser.add_argument('--mpterms', required=True,
                    help='path to the file containing mp term and associated gene\
                    1st column is gene name, 2nd column is mp_term')
parser.add_argument('--mpfunctions', required=True,
                    help='path to the file containing mp term and its functions')

args = parser.parse_args()


## Define the final holding variable for the FINAL PRODUCT
fin_case_df = pd.DataFrame()
fin_case_and_control_df = pd.DataFrame()

################ Generate two DataFrame- WithHeart Disease DataFrame and WithOUT heart disease DataFrame #########################


#Make a connection to the database
conn = sqlite3.connect(args.database)
#conn = sqlite3.connect('/Users/duongn/WorkFolder/WorkFolder/Mike_code/WEBSITE_JAWN/new.sqlite')

## Read in the "id_table" that holds the list of all sample IDs
df = pd.read_sql_query("SELECT * from id_table", conn)
matrix_temp = pd.read_sql_query("SELECT * from short_variant_mutation_matrix", conn)
#Make a copy of the gene column
gene_assoc = matrix_temp['GENE'].copy()
#Then drop the gene column from the table so the table returns to normal GATK vcf format
matrix_temp = matrix_temp.drop(columns='GENE')

conn.close()

mutation_info = matrix_temp.iloc[:,:9]
mutation_matrix = matrix_temp.iloc[:,9:]

filtered_col_names = pd.DataFrame(mutation_matrix.columns, columns=['ID'], dtype='int64')

merged_id = df.merge(filtered_col_names, on='ID' )

## Extract DISEASED and NON-DISEASED samples
## Assumption is that the 2nd column of the id_table is the boolean of 1 and 0
## and that the 1st column is the ID numbers/names
yes_heart = merged_id.loc[merged_id.iloc[:,1] == 1]
no_heart = merged_id.loc[merged_id.iloc[:,1] == 0]

## Get unique list
yes_sample_list = list(pd.unique(yes_heart.iloc[:,0]))
yes_sample_list_str = [str(i) for i in yes_sample_list]
no_sample_list = list(pd.unique(no_heart.iloc[:,0]))
no_sample_list_str = [str(i) for i in no_sample_list]



################## Separate main file into 2 cohorts ##################

yes_heart_matrix4 = mutation_matrix[yes_sample_list_str]
no_heart_matrix4 = mutation_matrix[no_sample_list_str]


high_moderate_mutation_info = mutation_info.copy()
high_moderate_matrix_yes = yes_heart_matrix4.copy()
high_moderate_matrix_no = no_heart_matrix4.copy()


## Read in list of genes associated with each mutation from the Database
## Put these genes into a list to itterate over later for mutation load analysis
small_gene_list = [] 
gene_assoc.map(lambda x: small_gene_list.extend(x.split(',')))

#Put these genes from string into a list of string.
gene_assoc = gene_assoc.map(lambda x: x.split(','))
## Make a unique list of genes to iterate over later on
## Get rid of '' entries using list comprehension
u_small_gene_list = np.unique([i for i in small_gene_list if i])

####################
human_MP = pd.read_csv(args.mpterms, delimiter='\t')
#human_MP = pd.read_csv('/Volumes/duongn/DBHi/WorkFolder/Mike_code/Heart_Defect/human_MP_ws.txt', delimiter='\t')
mp_def = pd.read_csv(args.mpfunctions, delimiter='\t')
#mp_def = pd.read_csv('/Volumes/duongn/DBHi/WorkFolder/Mike_code/Heart_Defect/mp_def.txt', delimiter='\t')


u_list_MP = pd.Series(np.unique(human_MP['MP']))
u_list_MP_gene = pd.Series(u_list_MP.copy())


for u_MP_index,u_MP in enumerate(u_list_MP):
    bob = list((human_MP.loc[human_MP['MP'] == u_MP])['SYMBOL'])
    u_list_MP_gene.at[u_MP_index] = bob
    #print(u_MP_index)

####################################################################################################
######### This chunk is ######################## MP term analysis ##################################
####################################################################################################

## Making holder matrixes for pval, stats,
## and for the #of mutaiton in DISEASED and NON-DISEASED
gene_pval = u_list_MP.copy()
gene_stat = u_list_MP.copy()
gene_yes_sum = u_list_MP.copy()
gene_no_sum = u_list_MP.copy()

for gene_eval_index,gene_eval in enumerate(u_list_MP_gene):
    ## Find, for each MP term, if the gene association with that MP term is associated with any mutations
    ## in the mutation matrix
    gene_logical = gene_assoc.map(lambda x : True if set(gene_eval).intersection(set(x)) else False)
    ## Pull those rows out from DISEASED and NON-DISEASED matrix
    gene_yes = high_moderate_matrix_yes.loc[gene_logical]
    gene_no = high_moderate_matrix_no.loc[gene_logical]
    ## Create a matrix to store the number of mutations per gene
    ## for DISEASED and NON-DISEASED
    gene_yes_sum[gene_eval_index] = gene_yes.sum(axis=0).sum()
    gene_no_sum[gene_eval_index] = gene_no.sum(axis=0).sum()
    #print(gene_eval_index)
    try:
        gene_stat[gene_eval_index],gene_pval[gene_eval_index] = stats.mannwhitneyu(gene_yes.sum(axis=0),gene_no.sum(axis=0))
    except:
        ## If dealing with zero, mannwhitneyU throws errors. So we put it as an exception
        gene_stat[gene_eval_index] = 0
        gene_pval[gene_eval_index] = 0



## Turn gene, pval, stats, yes, and no count (DISEASED & NON-DISEASED) to DF
u_list_MP_df = pd.DataFrame(u_list_MP)
gene_pval_df = pd.DataFrame(gene_pval)
gene_stat_df = pd.DataFrame(gene_stat)
gene_yes_sum_df = pd.DataFrame(gene_yes_sum)
gene_no_sum_df = pd.DataFrame(gene_no_sum)

## Put em all in a single DataFrame
add_df = pd.concat([u_list_MP_df,gene_pval_df,gene_stat_df,gene_yes_sum_df,gene_no_sum_df],axis=1)
add_df.columns = ['MP','pval','stat','yes_count','no_count']

## Remove all of the zero values, those are the ones that met the try-EXCEPT loop above
fin_add_df = (add_df.loc[add_df['pval'] != 0]).sort_values(by='pval')

## convert muatation count columns to numerics
fin_add_df['yes_count'] = pd.to_numeric(fin_add_df['yes_count'])
fin_add_df['no_count'] = pd.to_numeric(fin_add_df['no_count'])

## Calculate ratio to know which ones are OVERREPRESENTED in DISEASED or NON-DISEASED
ratio_df = (fin_add_df['yes_count']/yes_heart_matrix4.shape[1])/(fin_add_df['no_count']/no_heart_matrix4.shape[1])

## Add the ratio column to the DataFrame
fin_add_df = pd.concat([fin_add_df,ratio_df],axis=1)

## Name the columns
fin_add_df.columns = ['MP','pval','stat','yes_count','no_count','ratio']

## Convert the pval columns to nubmers to be sorted
fin_add_df['pval'] = pd.to_numeric(fin_add_df['pval'])

## Produce the full database to export to sqlite3 database
fin_add_df = fin_add_df.merge(mp_def,on=['MP'])

## Extract ENRICHED in DISEASED cohort
cases_enriched = (fin_add_df.loc[fin_add_df['ratio'] > 1]).sort_values(by='pval')


#$# Count how many genes are enriched in DISEASED cohort
(cases_enriched.loc[(cases_enriched['pval'] < 0.05)])['MP'].shape
## Extract ENRICHED in NON-DISEASED cohort
controls_enriched = (fin_add_df.loc[fin_add_df['ratio'] < 1]).sort_values(by='pval')

## Put the ENRICHED in DISEASED cohort in the FINAL df
fin_case_df = pd.concat([fin_case_df,cases_enriched])
## Same thing for the full matrix with enrichment in DISEASED and NON-DISEASED
fin_case_and_control_df = pd.concat([fin_case_and_control_df,fin_add_df])

## Reset the index so it looks nice
fin_case_df = fin_case_df.reset_index(drop=True)
fin_case_and_control_df = (fin_case_and_control_df.reset_index(drop=True)).sort_values(by='MP')

## Extract those over-represented in Cases and are lower than 0.05 
variant_cases_pval = fin_case_df[(fin_case_df['ratio']>1) & (fin_case_df['pval']<0.05) ]




############################ Export Data #########################################
## Make connection to the SQLite3 DataBase
conn = sqlite3.connect(args.database)

## Insert the computed dataframe as a table in the sqlite3 database
fin_add_df.to_sql('short_variant_MP', conn, if_exists="replace", index=False)

## close the connection to save data
conn.close()

variant_cases_pval.to_csv(sys.stdout, sep='\t', index=False)


