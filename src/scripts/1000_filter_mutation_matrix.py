import numpy as np
import pandas as pd
import sys
import sqlite3
import argparse

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

parser = argparse.ArgumentParser(description="""Code to generate a matrix that has mutation counts for
                                                later analyses""") 
parser.add_argument('--input', required=True,
                    help='path to the input .vcf file, the file filtered for high-moderate impact')
parser.add_argument('--database', required=True,
                    help='path to the sqlite3 database file that holds the sample names of the two cohorts')

args = parser.parse_args()


_,matrix_temp = read_vcf(args.input)

#### If the input file already has the mutation info
## Check to see if file has the correct header for column 1. If not, gotta check the file again
if matrix_temp.columns[0] == '#CHROM':
    mutation_info = matrix_temp.iloc[:,:9]
    matrix = matrix_temp.iloc[:,9:]
else:
    ## Throw error and exit the code
    sys.stderr.write('Please check the format of the input file. The headers shoulder be')
    sys.stderr.write('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO','FORMAT, then individual samples')
    sys.exit()


################################# Yes heart ####################################

# Generate masks
# Find out which cell has 0,1,2 mutation counts
heart_ones = matrix.applymap(lambda x: x.startswith('1/0') | x.startswith('0/1') | x.startswith('2/0') | x.startswith('0/2'))
heart_twos = matrix.applymap(lambda x: x.startswith('1/1') | x.startswith('2/2') | x.startswith('2/1') | x.startswith('1/2'))
heart_zeros = ~(np.logical_or(heart_ones,heart_twos))


# Apply masks to have either 0 or 1 or 2 in each cell.
# This then will be used to 
heart_matrix2 = matrix.mask(heart_zeros,0)
heart_matrix3 = heart_matrix2.mask(heart_ones,1)
heart_matrix4 = heart_matrix3.mask(heart_twos,2)


################################# Filter outliers >1000 mutations #######################

## Make a new matrix for DISEASED samples with <1000 mutations
sum_mat = heart_matrix4.sum(axis=0)
sum_mat_res = (sum_mat.loc[sum_mat <= 1000]).index
less1000_heart_matrix4 = heart_matrix4[sum_mat_res]

mutation_info_and_heart_matrix4 = pd.concat([mutation_info,less1000_heart_matrix4],axis=1)


conn = sqlite3.connect(args.database)
## Insert the computed dataframe as a table in the sqlite3 database
mutation_info_and_heart_matrix4.to_sql('short_variant_mutation_matrix', conn, if_exists="replace", index=False)
## close the connection to save data
conn.close()

mutation_info_and_heart_matrix4.to_csv(sys.stdout, sep='\t', index=False)