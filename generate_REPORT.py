from flask import Flask
from flask import render_template
import sqlite3
import argparse
app = Flask(__name__)
 
parser = argparse.ArgumentParser(description="""This script performs MetaP fisher test for Affy and WES data""") 

parser.add_argument('--database', required=True,
                    help='path to the sqlite3 database file that holds the sample names of the two cohorts')


args = parser.parse_args() 

@app.route('/')
def hello_world():
    # sqlalchemy commands to get table data:
    #snp_gene_cases = db.execute("SELECT * FROM snp_gene_cases)
    #return render_template("homepage.html", snp_gene_cases=snp_gene_cases)
    conn = sqlite3.connect(args.database)
    cur = conn.cursor()
    id_table = cur.execute("SELECT * FROM id_table").fetchall()
    SNP_gene_cases_table = cur.execute("SELECT * FROM short_variant_gene2gene WHERE ratio>1 AND pval<0.05").fetchall()
    SNP_go_cases_table = cur.execute("SELECT * FROM short_variant_GO WHERE ratio>1 AND pval<0.05").fetchall()
    SNP_mp_cases_table = cur.execute("SELECT * FROM short_variant_MP WHERE ratio>1 AND pval<0.05").fetchall()
    Affy_gene_cases_table = cur.execute("SELECT * FROM CNV_gene2gene WHERE ratio>1 AND pval<0.05").fetchall()
    Affy_go_cases_table = cur.execute("SELECT * FROM CNV_GO WHERE ratio>1 AND pval<0.05").fetchall()
    Affy_mp_cases_table = cur.execute("SELECT * FROM CNV_MP WHERE ratio>1 AND pval<0.05").fetchall()
    SNP_Affy_gene_cases_table = cur.execute("SELECT * FROM MetaP_gene2gene  WHERE ratio_x>1 AND ratio_y>1 AND combined_pval<0.05").fetchall()
    SNP_Affy_go_cases_table = cur.execute("SELECT * FROM MetaP_GO WHERE ratio_x>1 AND ratio_y>1 AND combined_pval<0.05").fetchall()
    SNP_Affy_mp_cases_table = cur.execute("SELECT * FROM MetaP_MP WHERE ratio_x>1 AND ratio_y>1 AND combined_pval<0.05").fetchall()
    id = cur.execute("SELECT DISTINCT ID FROM id_table").fetchall()
    new_id = [entry[0] for entry in id]
    return render_template("homepage.html", id_table=id_table, SNP_gene_cases_table=SNP_gene_cases_table, SNP_go_cases_table=SNP_go_cases_table, SNP_mp_cases_table=SNP_mp_cases_table,\
    Affy_gene_cases_table=Affy_gene_cases_table, Affy_go_cases_table=Affy_go_cases_table, Affy_mp_cases_table=Affy_mp_cases_table, SNP_Affy_gene_cases_table=SNP_Affy_gene_cases_table,\
    SNP_Affy_go_cases_table=SNP_Affy_go_cases_table, SNP_Affy_mp_cases_table=SNP_Affy_mp_cases_table, new_id=new_id)


#@app.route('/')
#def bob():
#    conn = sqlite3.connect('/Volumes/duongn/DBHi/WorkFolder/Mike_code/WEBSITE_JAWN/chd.sqlite')
#    cur = conn.cursor()
#    SNP_gene_cases_table = cur.execute("SELECT * FROM SNP_gene_cases")
#    return render_template("homepage.html", SNP_gene_cases_table=SNP_gene_cases_table) 

 
@app.route('/<int:id>/hello_link')
def hello_link(id):
    # connect to the sqlite db
    # get stuff from the table
    # put into a list or something you can iterate over
    # this_is_a_list_of_sample_names = db.execute("SELECT sample_name FROM samples")
    #return render_template("link.html", this_is_a_list_of_sample_names=this_is_a_list_of_sample_names)
    # in html, you iterate over bob to display the elements.
    
    conn = sqlite3.connect(args.database)
    cur = conn.cursor()
    short_variant_statement = '''SELECT "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT" FROM short_variant_mutation_matrix WHERE "{}" != 0'''.format(id)
    short_variant_of_id = cur.execute(short_variant_statement).fetchall()
    CNV_statement = '''SELECT * FROM CNV_mutation_matrix WHERE "WES from U. Washington file ID" = {} '''.format(id)
    CNV_of_id = cur.execute(CNV_statement).fetchall()
    some_variable = 2
    bob = [1,2,3,4,5]
    return render_template("link.html", some_variable=some_variable, bob=bob, id=id, short_variant_of_id=short_variant_of_id, CNV_of_id=CNV_of_id)


if __name__ == "__main__":
    app.run(debug=True, port=5000)
