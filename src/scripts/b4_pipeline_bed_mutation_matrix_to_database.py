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

#bed_file = pd.read_csv('/Users/duongn/WorkFolder/WorkFolder/WES/pipe/data/bevCombind_cnv.bed', delimiter='\t', dtype={'CHR': str, 'START': int, 'END': int})
bed_file = pd.read_csv(args.bed_file, delimiter='\t', dtype={'CHR': str, 'START': int, 'END': int})
bed_file['CHR'] = bed_file['CHR'].apply( lambda x : x.strip('chr'))

gene_assoc = (bed_file.iloc[:,0].copy()).astype(str)

## Read in list of genes in the hg19 genome
#gene_list = pd.read_csv('/Users/duongn/WorkFolder/WorkFolder/WES/pipe/data/ref-data/hg19_genes.bed', delimiter='\t')
gene_list = pd.read_csv(args.genelist, delimiter='\t')
gene_list['chr'] = gene_list['chr'].apply( lambda x : x.strip('chr'))

############################################ Gene pval test LOOP ######################################
## For each mutation, find which gene that mutation is in
## Create a matrix of the same size as the mutation matrix
## For each row, get the genes in "hg19 gene DATABASE" that overlap that mutation

small_gene_list = []
for row in bed_file.itertuples():
    gene_df = gene_list.loc[(gene_list['chr'] == row.CHR) & (gene_list['start'] <= row.END) & (gene_list['end'] >= row.START)]
    #$# This line doesn't save anything.... may need to delete or change
    #np.unique(gene_df['name2'])
    gene_assoc.at[row.Index] = list(gene_df['gene'])


gene_assoc_df = pd.DataFrame(gene_assoc.map(lambda x: ','.join(x)))
gene_assoc_df.columns = ['GENE']

bed_file_w_gene_assoc = pd.concat([gene_assoc_df,bed_file], axis=1)

conn = sqlite3.connect(args.database)

## Insert the computed dataframe as a table in the sqlite3 database
bed_file_w_gene_assoc.to_sql('CNV_mutation_matrix', conn, if_exists="replace", index=False)

## close the connection to save data
conn.close()

bed_file_w_gene_assoc.to_csv(sys.stdout, sep='\t', index=False)