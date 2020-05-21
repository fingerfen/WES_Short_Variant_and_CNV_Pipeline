import numpy as np
import pandas as pd
import scipy.stats as stats
import sqlite3
import sys
import argparse

## Define arg parser to take in inputs
parser = argparse.ArgumentParser(description="""This script performs MetaP fisher test for Affy and WES data""") 

parser.add_argument('--database', required=True,
                    help='path to the sqlite3 database file that holds the sample names of the two cohorts')


args = parser.parse_args()


################### Make connection to Database ################
#conn = sqlite3.connect('/Users/duongn/WorkFolder/WorkFolder/WES/pipe/data/endpoints/new.sqlite')
conn = sqlite3.connect(args.database)

########################################################################
############### Combine pval of Gene-gene analysis #####################
########################################################################

#Import Short variant and CNV tables
gene_result_full = pd.read_sql_query("SELECT * from short_variant_gene2gene", conn)
gene_Affy_result_full = pd.read_sql_query("SELECT * from CNV_gene2gene", conn)

## Merge the gene between the two
gene_merged = gene_result_full.merge(gene_Affy_result_full, on='gene')

#Make a new place holder dataframe to hold the stat and combined PVALUES
new_gene_pval = gene_merged['gene'].copy().rename('combined_pval')
new_gene_stat = gene_merged['gene'].copy().rename('combined_stat')

#combine PVALUES
for row in gene_merged.itertuples():
    [new_gene_stat.at[row.Index],new_gene_pval.at[row.Index]] = stats.combine_pvalues([row.pval_x,row.pval_y],method='fisher')

#Concat the combined pval to the merged table
gene_merged_new_gene_pval = pd.concat([gene_merged,new_gene_pval], axis=1)

#Get rid of any pval = zeros, which at this point shouldn't be any.
#This step is most likely not needed here
gene_merged_new_gene_pval_no_zeros = gene_merged_new_gene_pval.loc[gene_merged_new_gene_pval['combined_pval'] != 0]

#Get all the ones with pval lower than 0.05
gene_merged_new_gene_pval_no_zeros_significant = gene_merged_new_gene_pval_no_zeros.loc[gene_merged_new_gene_pval_no_zeros['combined_pval'] < 0.05]

#Get all of the one with ratio > 1 AKA over-represented in EXPERIEMNTAL group
#in both Short Variants and CNVs
gene_both_ratio_big = gene_merged_new_gene_pval_no_zeros_significant.loc[(gene_merged_new_gene_pval_no_zeros_significant['ratio_x'] >= 1) & (gene_merged_new_gene_pval_no_zeros_significant['ratio_y'] >= 1)]


########################################################################
############### Combine pval of GO term analysis #######################
########################################################################

GO_result_full = pd.read_sql_query("SELECT * from short_variant_GO", conn)
GO_Affy_result_full = pd.read_sql_query("SELECT * from CNV_GO", conn)

#Merge
GO_merged = GO_result_full.merge(GO_Affy_result_full, on='GO')

#Create new place holder
new_GO_pval = GO_merged['GO'].copy().rename('combined_pval')
new_GO_stat = GO_merged['GO'].copy().rename('combined_stat')
#combine the pvals
for row in GO_merged.itertuples():
    [new_GO_stat.at[row.Index],new_GO_pval.at[row.Index]] = stats.combine_pvalues([row.pval_x,row.pval_y],method='fisher')

#Concat the pval to the merged table
GO_merged_new_GO_pval = pd.concat([GO_merged,new_GO_pval], axis=1)
#Get rid of any pval = zero
GO_merged_new_GO_pval_no_zeros = GO_merged_new_GO_pval.loc[GO_merged_new_GO_pval['combined_pval'] != 0]
#Take only those with pval < 0.05
GO_merged_new_GO_pval_no_zeros_significant = GO_merged_new_GO_pval_no_zeros.loc[GO_merged_new_GO_pval_no_zeros['combined_pval'] < 0.05]
#Take only those over represented in both Short Varianta nd CNVs
GO_both_ratio_big = GO_merged_new_GO_pval_no_zeros_significant.loc[(GO_merged_new_GO_pval_no_zeros_significant['ratio_x'] >= 1) & (GO_merged_new_GO_pval_no_zeros_significant['ratio_y'] >= 1)]



########################################################################
############### Combine pval of MP term analysis #######################
########################################################################



MP_result_full = pd.read_sql_query("SELECT * from short_variant_MP", conn)
MP_Affy_result_full = pd.read_sql_query("SELECT * from CNV_MP", conn)

#Merge the data
MP_merged = MP_result_full.merge(MP_Affy_result_full, on='MP')

#Create new place holder for stats and combined pvals
new_MP_pval = MP_merged['MP'].copy().rename('combined_pval')
new_MP_stat = MP_merged['MP'].copy().rename('combined_stat')
#Combine the pvals
for row in MP_merged.itertuples():
    [new_MP_stat.at[row.Index],new_MP_pval.at[row.Index]] = stats.combine_pvalues([row.pval_x,row.pval_y],method='fisher')

#Concat the pval to the merged table
MP_merged_new_MP_pval = pd.concat([MP_merged,new_MP_pval], axis=1)
#Get rid of any pval= zero
MP_merged_new_MP_pval_no_zeros = MP_merged_new_MP_pval.loc[MP_merged_new_MP_pval['combined_pval'] != 0]
#Take only those with pval < 0.05
MP_merged_new_MP_pval_no_zeros_significant = MP_merged_new_MP_pval_no_zeros.loc[MP_merged_new_MP_pval_no_zeros['combined_pval'] < 0.05]
#Take only those over represented in both Short Varainta and CNVs
MP_both_ratio_big = MP_merged_new_MP_pval_no_zeros_significant.loc[(MP_merged_new_MP_pval_no_zeros_significant['ratio_x'] >= 1) & (MP_merged_new_MP_pval_no_zeros_significant['ratio_y'] >= 1)]

## Write the results to the database
gene_both_ratio_big.to_sql('MetaP_gene2gene', conn, if_exists="replace", index=False)
GO_both_ratio_big.to_sql('MetaP_GO', conn, if_exists="replace", index=False)
MP_both_ratio_big.to_sql('MetaP_MP', conn, if_exists="replace", index=False)

#Close the connection to save the data. 
conn.close()

sys.stdout.write('Check your sqlite3 database file for the resulted tables named: SNP_Affy_gene_cases, SNP_Affy_go_cases, and SNP_Affy_mp_cases')