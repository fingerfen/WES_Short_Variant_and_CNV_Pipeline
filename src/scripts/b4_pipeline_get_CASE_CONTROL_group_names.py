################## PURPOSE #################

# This is a script to take in data from CHOP to generate a table that has 
# sample name on column 1 and whether or not the sample is in the CASES or CONTROLS
# group in the 2nd column. 
# this script then put that data as a table named 'id_table' in a sqlite3 database

# You can REPLACE this script by manualy making a sqlite3 database
# and put in it a table named 'id_table'. The table should have 2 columns
# The 1st column is the sample names/IDs
# The 2nd column is either the number 0 or 1.
# 0 is for samples belonging in CONTROLS group
# 1 is for samples belonging in CASES group

################# ASSUMPTION ###############

## Very BIG assumptions.
## Is that you are working with the same 22q11.2 DS samples as me. 
## Most likely, you would need to write your own script for this step. 

############################################

import pandas as pd
import sqlite3



############# Merge heart anno, ethnic tables , and chr22 list together #################
# Merge heart_anno and ethnic

heart_anno = pd.read_csv('/Users/duongn/WorkFolder/WorkFolder/Mike_code/Heart_Defect/wes_chd_anno.txt', delimiter='\t')
ethnic_info = pd.read_csv('/Users/duongn/WorkFolder/WorkFolder/Mike_code/Heart_Defect/ethnicGroup.txt', delimiter='\t')
list_full = pd.read_csv('/Users/duongn/WorkFolder/WorkFolder/Mike_code/Heart_Defect/list_full', delimiter='\t', header=None)


a = heart_anno.merge(ethnic_info,left_on='famID', right_on='famid')


# Extract only the caucasian from the table
b = a.loc[a['cluster'] == 'Cau']


# Merge the table with the Exome samples
c = list_full.merge(b,left_on=0, right_on='WES from U. Washington file ID')

## Establish a connection to a sqlite3 database
conn = sqlite3.connect(args.database)

## Create a table to be put into a database
id_table = c[['WES from U. Washington file ID','heart6']]
id_table.columns = ['ID','Cohort']

## Export the DataFrame to a table in the database
id_table.to_sql('id_table', conn, if_exists="replace")

conn.close()