################## PURPOSE #################

# This is a script to perform the Affy gene-gene analysis to compare between CASES vs CONTROLS groups
# this looks at all inputted mutations and locate which gene has which mutation.
# Then every single gene is looked at to see if there is an over represented number of mutations
# in the CASES group vs the CONTROLS group

################# ASSUMPTION ###############

# That the 0th column to and including the 8th column are the info columns. They are formatted as:
# '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO','FORMAT,
# After that, the 9th column would contain info on the 1st sample, 10th column on the 2nd sample and so on

# That there is a sqlite3 database available that holds a table named 'id_table'
# This table has two columns: 
# 1st column = sample names/IDs
# 2nd column = either 0 or 1 to indicate if sample belongs in 0 (contorls) or 1 (cases) group.


############################################

## Import necessary packages
import numpy as np
import pandas as pd
import sys
import scipy.stats as stats
import sqlite3
import argparse



## Define arg parser to take in inputs
parser = argparse.ArgumentParser(description="""Gene-Gene analysis: This script
                                                performs the gene to gene analysis
                                                by taking a look at mutations called
                                                in each gene and compare the 
                                                Experimental vs Control groups""") 

parser.add_argument('--bed_file', required=True,
                    help='path to the bed file with columns as chr,start,end,sampleid')
parser.add_argument('--database', required=True,
                    help='path to the sqlite3 database file that holds the sample names of the two cohorts')
parser.add_argument('--genelist', required=True,
                    help='path to the file containing all genes name, start, and end location in the human genome.\
                    make sure its version matches with your human reference genome \(i.e hg19/hg38\)')



args = parser.parse_args()


## Define the final holding variable for the FINAL PRODUCT
fin_case_df = pd.DataFrame()
fin_case_and_control_df = pd.DataFrame()

################ Generate two DataFrame- WithHeart Disease DataFrame and WithOUT heart disease DataFrame #########################


#Make a connection to the database
conn = sqlite3.connect(args.database)
#conn = sqlite3.connect('/Users/duongn/WorkFolder/WorkFolder/WES/pipe/data/endpoints/new.sqlite')

## Read in the "id_table" that holds the list of all sample IDs
df = pd.read_sql_query("SELECT * from id_table", conn)
matrix_temp = pd.read_sql_query("SELECT * from short_variant_mutation_matrix", conn)

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

##Generate a list of sample name stored in a DataFrame to be merged with CNV list later on
yes_sample = pd.DataFrame(pd.unique(yes_heart.iloc[:,0]))
yes_sample.columns = ['WES from U. Washington file ID']
no_sample = pd.DataFrame(pd.unique(no_heart.iloc[:,0]))
no_sample.columns = ['WES from U. Washington file ID']

bed_file = pd.read_csv(args.bed_file, delimiter='\t', dtype={'chr': str, 'start': int, 'end': int})
#bed_file = pd.read_csv('/Users/duongn/WorkFolder/WorkFolder/WES/pipe/data/bev_affy.bed', delimiter='\t', dtype={'chr': str, 'start': int, 'end': int})
bed_file['chr'] = bed_file['chr'].apply( lambda x : x.strip('chr'))

gene_assoc = (bed_file.iloc[:,0].copy()).astype(str)

## Read in list of genes in the hg19 genome
#gene_list = pd.read_csv(args.genelist, delimiter='\t')
gene_list = pd.read_csv(args.genelist, delimiter='\t', dtype={'chr': str, 'start': int, 'end': int})
#gene_list = pd.read_csv('/Users/duongn/WorkFolder/WorkFolder/WES/pipe/data/ref-data/hg19_genes.bed', delimiter='\t')
gene_list['chr'] = gene_list['chr'].apply( lambda x : x.strip('chr'))





############################################ Gene pval test LOOP ######################################
## For each mutation, find which gene that mutation is in
## Create a matrix of the same size as the mutation matrix
## For each row, get the genes in "hg19 gene DATABASE" that overlap that mutation

small_gene_list = []
for row in bed_file.itertuples():
    gene_df = gene_list.loc[(gene_list['chr'] == row.chr) & (gene_list['start'] <= row.end) & (gene_list['end'] >= row.start)]
    #$# This line doesn't save anything.... may need to delete or change
    #np.unique(gene_df['name2'])
    gene_assoc.at[row.Index] = list(gene_df['gene'])
    small_gene_list += list(gene_df['gene'])
    
    
## Make a unique list of genes to iterate over later on
u_small_gene_list = np.unique(small_gene_list)


####################################################################################################
######### This chunk is ######################## Affymetrix gene2gene analysis #####################
####################################################################################################


LOF_gene_df = pd.DataFrame(u_small_gene_list)
gene_pval = LOF_gene_df[0].copy()
gene_stat = gene_pval.copy()
gene_yes_sum = gene_pval.copy()
gene_no_sum = gene_pval.copy()




for gene_eval_index,gene_eval in enumerate(u_small_gene_list):

    # Find affy_res that has a certain gene
    gene_logical = gene_assoc.map(lambda x : True if gene_eval in x else False)
    affy_gene_res = bed_file.loc[gene_logical]

    # Merge and take only those patients who has the same WES ID to ensure comparability

    # split into those YES ( with CHD) and NO (without CHD)
    yes_affy_wes_res = affy_gene_res.merge(yes_sample,on='WES from U. Washington file ID')
    no_affy_wes_res = affy_gene_res.merge(no_sample,on='WES from U. Washington file ID')

    # Turn the counts into a series, count how many mutations does each sample have
    yes_affy_mutations = yes_affy_wes_res['WES from U. Washington file ID'].value_counts()
    no_affy_mutations = no_affy_wes_res['WES from U. Washington file ID'].value_counts()

    # Create the zeros values for those samples who don't have a mutation. Concat em, then input into mannwhiteney
#    stat_test_yes = pd.concat([pd.Series(np.zeros([yes_sample.shape[0] - yes_affy_mutations.shape[0],])),yes_affy_mutations])
#    stat_test_no = pd.concat([pd.Series(np.zeros([no_sample.shape[0] - no_affy_mutations.shape[0],])),no_affy_mutations])
    stat_test_yes = pd.concat([yes_affy_mutations,pd.Series(np.zeros([yes_sample.shape[0] - yes_affy_mutations.shape[0],]))])
    stat_test_no = pd.concat([no_affy_mutations,pd.Series(np.zeros([no_sample.shape[0] - no_affy_mutations.shape[0],]))])
    gene_yes_sum[gene_eval_index] = stat_test_yes.sum()
    gene_no_sum[gene_eval_index] = stat_test_no.sum()
    try:
        gene_stat[gene_eval_index],gene_pval[gene_eval_index] = stats.mannwhitneyu(stat_test_yes,stat_test_no)
    except:
        gene_stat[gene_eval_index] = 0
        gene_pval[gene_eval_index] = 0


gene_pval_df = pd.DataFrame(gene_pval)
gene_stat_df = pd.DataFrame(gene_stat)
gene_yes_sum_df = pd.DataFrame(gene_yes_sum)
gene_no_sum_df = pd.DataFrame(gene_no_sum)


add_df = pd.concat([LOF_gene_df,gene_pval_df,gene_stat_df,gene_yes_sum_df,gene_no_sum_df],axis=1)
add_df.columns = ['gene','pval','stat','yes_count','no_count']

fin_add_df = (add_df.loc[add_df['pval'] != 0]).sort_values(by='pval')

## convert muatation count columns to numerics
fin_add_df['yes_count'] = pd.to_numeric(fin_add_df['yes_count'])
fin_add_df['no_count'] = pd.to_numeric(fin_add_df['no_count'])

## Calculate ratio to know which ones are OVERREPRESENTED in DISEASED or NON-DISEASED
ratio_df = (fin_add_df['yes_count']/yes_sample.shape[0])/(fin_add_df['no_count']/no_sample.shape[0])

## Add the ratio column to the DataFrame
fin_add_df = pd.concat([fin_add_df,ratio_df],axis=1)

## Name the columns
fin_add_df.columns = ['gene','pval','stat','yes_count','no_count','ratio']


## Convert the pval columns to nubmers to be sorted
fin_add_df['pval'] = pd.to_numeric(fin_add_df['pval'])
## Extract ENRICHED in DISEASED cohort
cases_enriched = (fin_add_df.loc[fin_add_df['ratio'] > 1]).sort_values(by='pval')

#$# Count how many genes are enriched in DISEASED cohort
(cases_enriched.loc[(cases_enriched['pval'] < 0.05)])['gene'].shape
## Extract ENRICHED in NON-DISEASED cohort
controls_enriched = (fin_add_df.loc[fin_add_df['ratio'] < 1]).sort_values(by='pval')

## Put the ENRICHED in DISEASED cohort in the FINAL df
fin_case_df = pd.concat([fin_case_df,cases_enriched])
## Same thing for the full matrix with enrichment in DISEASED and NON-DISEASED
fin_case_and_control_df = pd.concat([fin_case_and_control_df,fin_add_df])

## Reset the index so it looks nice
fin_case_df = fin_case_df.reset_index(drop=True)
fin_case_and_control_df = (fin_case_and_control_df.reset_index(drop=True)).sort_values(by='gene')

## Extract those over-represented in Cases and are lower than 0.05 
variant_cases_pval = fin_case_df[(fin_case_df['ratio']>1) & (fin_case_df['pval']<0.05) ]



############################ Export Data #########################################
## Make connection to the SQLite3 DataBase
conn = sqlite3.connect(args.database)

## Insert the computed dataframe as a table in the sqlite3 database
fin_add_df.to_sql('CNV_gene2gene', conn, if_exists="replace", index=False)

## close the connection to save data
conn.close()

variant_cases_pval.to_csv(sys.stdout, sep='\t', index=False)
