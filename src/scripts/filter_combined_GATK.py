import numpy as np
import pandas as pd
import sys
import re
import os
import io


######################## NOTE ############################

'''
This file only filters and classify the result into certain 
class and change how the results are displayed. 

There is no cutting or adding of informations
'''

########################################################

#### Read_vcf addapted from online from dceoy/read_vcf.py

def read_vcf(path):
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


########################### Code body ###########################


## Read in file input
gatk_name = sys.argv[1]

sys.stderr.write('Importing gatk file\n')

## Call the "read_vcf" function and output the header and the actual content
header_lines, gatk_file = read_vcf(gatk_name)


###################### 5 cols format processing ##############

## Split up the file. Since there are two ways of answwer format
gatk_file_5cols = gatk_file.loc[gatk_file['FORMAT'] == 'GT:AD:DP:GQ:PL']


#split the file at the 9th column with the left side being the info
#and the right side being the actual information regarding the samples
gatk_matrix_5cols = gatk_file_5cols.iloc[:,9:]
gatk_mutation_info_5cols = gatk_file_5cols.iloc[:,:9]


## Report to stdout the status of the code
sys.stderr.write('Fitlering gatk file\n')


## Make a logical matrix for cells that begin with "./." and change them to 0/0
no_call_mask_5cols = gatk_matrix_5cols.applymap(lambda x : x.startswith('./.'))
gatk_matrix2_5cols = gatk_matrix_5cols.mask(no_call_mask_5cols,'0/0:0,0:0')

## Filter out the AD, if it has single like a dot"." dot entry, then remove it. 
AD_DP_mask_5cols = gatk_matrix2_5cols.applymap(lambda x: True if len(x.split(':')[1].split(',')) == 1 else False)
gatk_matrix2_v2_5cols = gatk_matrix2_5cols.mask(AD_DP_mask_5cols,'0/0:0,0:0')


### Filter for DP, isolate it first and put it into a variable
gatk_tot_5cols = gatk_matrix2_v2_5cols.applymap(lambda x : x.split(':')[2])


### Some has DP = .    so replace it with 0

gatk_tot_5cols = gatk_tot_5cols.mask(gatk_tot_5cols == '.','0')

### Get the mask for DP, if bigger than 8 == True, then all the false one gets replaced
tot_mask_5cols = gatk_tot_5cols.applymap(lambda x : True if int(x) >= 8 else False )

gatk_matrix3_5cols = gatk_matrix2_v2_5cols.where(tot_mask_5cols,'0/0:0,0:0')


### Filter for ALT, grab the second element of the Allel Depth column

gatk_alt_5cols = gatk_matrix3_5cols.applymap(lambda x : x.split(':')[1].split(',')[1])


## Report the situation to stdout
sys.stderr.write('Converting gatk file to int\n')

## Change the file to numeric
gatk_tot_5cols = gatk_tot_5cols.apply(pd.to_numeric)
gatk_alt_5cols = gatk_alt_5cols.apply(pd.to_numeric)

## Get rid of anything lesser than 0.25
alt_mask_5cols = gatk_alt_5cols/gatk_tot_5cols >= .25
gatk_matrix4_5cols = gatk_matrix3_5cols.where(alt_mask_5cols,'0/0:0,0:0')


gatk_fin_file_5cols = pd.concat([gatk_mutation_info_5cols,gatk_matrix4_5cols],axis=1)



###################### 3 cols format processing ##############

## Split up the file. Since there are two ways of answwer format
gatk_file_3cols = gatk_file.loc[gatk_file['FORMAT'] == 'GT:DP:PL']


#split the file up
gatk_matrix_3cols = gatk_file_3cols.iloc[:,9:]
gatk_mutation_info_3cols = gatk_file_3cols.iloc[:,:9]


## Report to stdout the status of the code
sys.stderr.write('Fitlering gatk file\n')


## Grab all mask that are of certain requirements and change them to 0/0
no_call_mask_3cols = gatk_matrix_3cols.applymap(lambda x : x.startswith('./.'))
gatk_matrix2_3cols = gatk_matrix_3cols.mask(no_call_mask_3cols,'0/0:0:0')

## Filter out the AD, if it has single like a dot"." dot entry, then remove it. 
# AD_DP_mask_3cols = gatk_matrix2_3cols.applymap(lambda x: True if len(x.split(':')[1].split(',')) == 1 else False)
# gatk_matrix2_v2_3cols = gatk_matrix2_3cols.mask(AD_DP_mask_3cols,'0/0:0:0')


### Filter for DP
### DP in this case is in the 2nd column, unlike the 5cols format from above
gatk_tot_3cols = gatk_matrix2_3cols.applymap(lambda x : x.split(':')[1])


### Some has DP = .    so replace it with 0

gatk_tot_3cols = gatk_tot_3cols.mask(gatk_tot_3cols == '.','0')

### Get the mask for DP, if bigger than 8 == True
tot_mask_3cols = gatk_tot_3cols.applymap(lambda x : True if int(x) >= 8 else False )

gatk_matrix3_3cols = gatk_matrix2_3cols.where(tot_mask_3cols,'0/0:0:0')

### Filter for ALT

# gatk_alt_3cols = gatk_matrix3_3cols.applymap(lambda x : x.split(':')[1].split(',')[1])


## Report the situation to stdout
sys.stderr.write('Converting gatk file to int\n')

## Change the file to numeric
# gatk_tot_3cols = gatk_tot_3cols.apply(pd.to_numeric)
# gatk_alt_3cols = gatk_alt_3cols.apply(pd.to_numeric)

# alt_mask_3cols = gatk_alt_3cols/gatk_tot_3cols >= .25
# gatk_matrix4_3cols = gatk_matrix3_3cols.where(alt_mask_3cols,'0/0:0:0')


gatk_fin_file_3cols = pd.concat([gatk_mutation_info_3cols,gatk_matrix3_3cols],axis=1)


################# Combine 5 and 3 cols results together ###############

combined_fin_file = (pd.concat([gatk_fin_file_5cols, gatk_fin_file_3cols],axis=0, ignore_index=False)).sort_index()



##################### Output to results ##############################


#Print out to stdout the location of the new file
sys.stderr.write('Writing new gatk file')

## Print out the header first
sys.stdout.write(''.join(header_lines))

## Print the content of the filtered file. 
combined_fin_file.to_csv(sys.stdout, sep='\t', index=False)