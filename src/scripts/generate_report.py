import sys
import os
from pathlib import Path
import sqlite3
import pandas as pd
import argparse


## Define arg parser to take in inputs
parser = argparse.ArgumentParser(description="""Script to generate a bunch of html files that together
												make up the report HTML page.""") 


parser.add_argument('--database', required=True,
                    help='path to the sqlite3 database file that holds the sample names of the two cohorts')
parser.add_argument('--project_summary', required=True,
					help='path to a simple text file with a few lines describing what the project is about')
parser.add_argument('--homepage_path', required=True,
					help='path where the output homepage html is to be at. Usually recommended to be in data/endpoints folder')

args = parser.parse_args()


###################################################################################
############################## Defining Functions #################################
###################################################################################

##Define a function to read in table query and turn it into HTML code
def get_html_table(sql_command_str, db_name):    
    import apsw
    import io
    output=io.StringIO()
    conn = apsw.Connection(db_name)
    shell=apsw.Shell(stdout=output, db=conn)
    # How to execute a dot command
    shell.process_command(".mode html")
    shell.process_sql(str(sql_command_str))
    return (output.getvalue())

## Define a function that create html linked entries in a datatable
def get_html_linked_table(sql_command_str,db_name):
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    ## Read in the "id_table" that holds the list of all sample IDs
    query = c.execute(sql_command_str).fetchall()
    conn.close()
    html_code = ''
    for i in query:
        html_code += '<tr>\n'
        html_code += '<td><a href="index/{}.html">{}</a></td>\n'.format(i[0],i[0])
        if i[1] == 0:
            group = 'Controls'
        elif i[1] == 1:
            group = "Cases"
        else:
            group = "Not applicable"
            
        html_code += '<td>{}</td>\n'.format(group)
        html_code += '</tr>\n'
    return html_code

##Define function to generate the subpages of the main html page.
def get_html_subpage(sample,db_name,path_to_index):
    html_code = '''(<!doctype html>
    <html lang="en">
    <head>
      <!--
      <link rel="stylesheet" href="http://cdn.datatables.net/1.10.20/css/jquery.dataTables.min.css">
    -->
      <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.4.1/jquery.min.js"></script>
      <script src="https://cdn.datatables.net/1.10.20/js/jquery.dataTables.min.js"></script>
      <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/4.1.3/css/bootstrap.css">
      <link rel="stylesheet" href="https://cdn.datatables.net/1.10.21/css/dataTables.bootstrap4.min.css">
      <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.4.1/jquery.min.js"></script>
      <script src="https://cdn.datatables.net/1.10.20/js/jquery.dataTables.min.js"></script>
      <script src="https://cdn.datatables.net/1.10.21/js/dataTables.bootstrap4.min.js"></script>

      <style>
        section {
        background: #181818;
        color: #ffffff;
        padding: 50px 0;
    }
      </style>
      <script>
        $(document).ready( function () {
          $('table.display').DataTable();
        } );
      </script>
    </head>


    <br>
    <br>
    <body>

    <div class="container">'''
    html_code += '''<h2><center>Sample {}</center></h2>'''.format(sample)


    html_code += '''<div>
    <table table border="1" class="display table table-hover">
      <thead class=thead-dark>
        <tr>
          <th class="tg-0lax">Category</th>
          <th class="tg-0lax">Details</th>
        </tr>
      </thead>
    <tbody>
      <tr>
        <td class="tg-0lax">Age</td>
        <td class="tg-0lax">2</td>
      </tr>
      <tr>
        <td class="tg-0pky">Ethnicity</td>
        <td class="tg-0lax">Caucasian</td>
      </tr>
      <tr>
        <td class="tg-0lax">Admission Date</td>
        <td class="tg-0lax">July 4th 2019</td>
      </tr>
      <tr>
        <td class="tg-0lax">Department</td>
        <td class="tg-0lax">Oncology</td>
      </tr>
      <tr>
        <td class="tg-0lax">Diagnosis</td>
        <td class="tg-0lax">22q11.2 DS</td>
      </tr>
    </tbody>
    </table>
    </div>
    <br>
    <br>


    <section>
    <h2><center>Short Variant Information</center></h2>
    </section>
    <br>
    <div>
    <table table border="1" class="display table table-hover">
      <thead class=thead-dark>
        <tr style="text-align: center;">
          <th>CHR</th>
          <th>POS</th>
          <th>REF</th>
          <th>ALT</th>
        </tr>
      </thead>
      <tbody align=center >'''
        
    html_code += get_html_table('''SELECT "#CHROM","POS","REF","ALT" FROM short_variant_mutation_matrix WHERE "{}" != 0'''.format(sample), db_name)
    html_code += '''  </tbody>
    </table>
    </div>

    <h4>Fisher Exact on GO term</h4>
    <div>
    <table table border="1" class="display table table-hover">
      <thead class=thead-dark>
        <tr style="text-align: center;">
          <th>GO</th>
          <th>ODD_RATIO</th>
          <th>PVALUE</th>
          <th>FUNCTION</th>
    <!--     
     <th>INFO</th>
          <th>FORMAT</th>
    -->
        </tr>
      </thead>
      <tbody align=center >'''
    html_code += get_html_table('''SELECT "GO","or_{}","pval_{}","name" FROM short_variant_sample_based_GO WHERE "pval_{}" < 0.05 AND "pval_{}" > 0'''.format(sample,sample,sample,sample), db_name)
    html_code += '''  </tbody>
    </table>
    </div>


    <h4>Fisher Exact on MP term</h4>
    <div>
    <table table border="1" class="display table table-hover">
      <thead class=thead-dark>
        <tr style="text-align: center;">
          <th>MP</th>
          <th>ODD_RATIO</th>
          <th>PVALUE</th>
          <th>FUNCTION</th>
    <!--     
     <th>INFO</th>
          <th>FORMAT</th>
    -->
        </tr>
      </thead>
      <tbody align=center >'''
    html_code += get_html_table('''SELECT "MP","or_{}","pval_{}","name" FROM short_variant_sample_based_MP WHERE "pval_{}" < 0.05 AND "pval_{}" > 0'''.format(sample,sample,sample,sample), db_name)
    html_code += '''  </tbody>
    </table>
    </div>

    <br>
    <br>
    <section>
    <h2><center>Copy Number Variation Information</center></h2>
    </section>

    <div>
    <table table border="1" class="display table table-hover">
      <thead class=thead-dark>
        <tr style="text-align: center;">
          <th>CHR</th>
          <th>START</th>
          <th>END</th>
          <th>CNV_TYPE</th>
        </tr>
      </thead>
      <tbody>'''
    html_code += get_html_table('''SELECT "CHR","START","END","CNV_TYPE" FROM CNV_mutation_matrix WHERE "ID" = {} '''.format(sample), db_name)
    html_code += '''  </tbody>
    </table>
    </div>

    <h4>Fisher Exact on GO term</h4>
    <div>
    <table table border="1" class="display table table-hover">
      <thead class=thead-dark>
        <tr style="text-align: center;">
          <th>GO</th>
          <th>ODD_RATIO</th>
          <th>PVALUE</th>
          <th>FUNCTION</th>
    <!--     
     <th>INFO</th>
          <th>FORMAT</th>
    -->
        </tr>
      </thead>
      <tbody align=center >'''
    html_code += get_html_table('''SELECT "GO","or_{}","pval_{}","name" FROM CNV_sample_based_GO WHERE "pval_{}" < 0.05 AND "pval_{}" > 0'''.format(sample,sample,sample,sample), db_name)
    html_code += '''  </tbody>
    </table>
    </div>


    <h4>Fisher Exact on MP term</h4>
    <div>
    <table table border="1" class="display table table-hover">
      <thead class=thead-dark>
        <tr style="text-align: center;">
          <th>MP</th>
          <th>ODD_RATIO</th>
          <th>PVALUE</th>
          <th>FUNCTION</th>
    <!--     
     <th>INFO</th>
          <th>FORMAT</th>
    -->
        </tr>
      </thead>
      <tbody align=center >'''
    html_code += get_html_table('''SELECT "MP","or_{}","pval_{}","name" FROM CNV_sample_based_MP WHERE "pval_{}" < 0.05 AND "pval_{}" > 0'''.format(sample,sample,sample,sample), db_name)
    html_code += '''  </tbody>
    </table>
      </div>
    </div>
    </body>
    </html>'''
    text_file = open("{}/{}.html".format(path_to_index,sample), "w")
    n = text_file.write(html_code)
    text_file.close()
    return n



###################################################################################
############################## BEGINNING OF SCRIPT ################################
###################################################################################



## Check if path exists, if not, create it. Hypothetically, the "data/endpoints/index" folder should already exist
## because the pipeline comes with the endpoints/index folder

## This should be the default location. 
path_for_index_file = "data/endpoints/index"


#Check to see if this path exists, if not, make it and its parent directories
Path(path_for_index_file).mkdir(parents=True, exist_ok=True)


## Defind the beginning of the HTML file.
head = '''<!doctype html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
  <!--
  <link rel="stylesheet" href="http://cdn.datatables.net/1.10.20/css/jquery.dataTables.min.css">
-->
  <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/4.1.3/css/bootstrap.css">
  <link rel="stylesheet" href="https://cdn.datatables.net/1.10.21/css/dataTables.bootstrap4.min.css">
  <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.4.1/jquery.min.js"></script>
  <script src="https://cdn.datatables.net/1.10.20/js/jquery.dataTables.min.js"></script>
  <script src="https://cdn.datatables.net/1.10.21/js/dataTables.bootstrap4.min.js"></script>


  <!-- the side bar links and scripts --> 
  <!-- Scrollbar Custom CSS -->
  <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/malihu-custom-scrollbar-plugin/3.1.5/jquery.mCustomScrollbar.min.css">
    <!-- jQuery Custom Scroller CDN -->
    <script src="https://cdnjs.cloudflare.com/ajax/libs/malihu-custom-scrollbar-plugin/3.1.5/jquery.mCustomScrollbar.concat.min.js"></script>
    <!-- Popper.JS -->
    <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.0/umd/popper.min.js" integrity="sha384-cs/chFZiN24E4KMATLdqdvsezGxaGsi4hLGOzlXwp5UZB1LY//20VyM2taTB4QvJ" crossorigin="anonymous"></script>
    <!-- Font Awesome JS -->
    <script defer src="https://use.fontawesome.com/releases/v5.0.13/js/solid.js" integrity="sha384-tzzSw1/Vo+0N5UhStP3bvwWPq+uvzCMfrN1fEFe+xBmv1C/AtVX5K0uZtmcHitFZ" crossorigin="anonymous"></script>
    <script defer src="https://use.fontawesome.com/releases/v5.0.13/js/fontawesome.js" integrity="sha384-6OIrr52G08NpOFSZdxxz1xdNSndlD4vdcf/q2myIUVO0VsqaGHJsB0RaBE01VTOY" crossorigin="anonymous"></script>



  <script>
    $(document).ready( function () {
      $('table.display').DataTable();
    } );
  </script>
  <style>
    @media (min-width: 1200px) {
    .container{
        max-width: 1400px;
    }
}

@import "https://fonts.googleapis.com/css?family=Poppins:300,400,500,600,700";
body {
    font-family: 'Poppins', sans-serif;
    background: #fafafa;
}

p {
    font-family: 'Poppins', sans-serif;
    font-size: 1.1em;
    font-weight: 300;
    line-height: 1.7em;
    color: #999;
}

section {
    background: #181818;
    color: #ffffff;
    padding: 50px 0;
}

a,
a:hover,
a:focus {
    color: inherit;
    text-decoration: none;
    transition: all 0.3s;
}

.navbar {
    padding: 15px 10px;
    background: #fff;
    border: none;
    border-radius: 0;
    margin-bottom: 40px;
    box-shadow: 1px 1px 3px rgba(0, 0, 0, 0.1);
}

.navbar-btn {
    box-shadow: none;
    outline: none !important;
    border: none;
}

.line {
    width: 100%;
    height: 1px;
    border-bottom: 1px dashed #ddd;
    margin: 40px 0;
}

/* ---------------------------------------------------
    SIDEBAR STYLE
----------------------------------------------------- */

#sidebar {
    width: 250px;
    position: fixed;
    top: 0;
    left: -250px;
    height: 100vh;
    z-index: 999;
    background: #7386D5;
    color: #fff;
    transition: all 0.3s;
    overflow-y: scroll;
    box-shadow: 3px 3px 3px rgba(0, 0, 0, 0.2);
}

#sidebar.active {
    left: 0;
}

#dismiss {
    width: 35px;
    height: 35px;
    line-height: 35px;
    text-align: center;
    background: #7386D5;
    position: absolute;
    top: 10px;
    right: 10px;
    cursor: pointer;
    -webkit-transition: all 0.3s;
    -o-transition: all 0.3s;
    transition: all 0.3s;
}

#dismiss:hover {
    background: #fff;
    color: #7386D5;
}

.overlay {
    display: none;
    position: fixed;
    width: 100vw;
    height: 100vh;
    background: rgba(0, 0, 0, 0.7);
    z-index: 998;
    opacity: 0;
    transition: all 0.5s ease-in-out;
}
.overlay.active {
    display: block;
    opacity: 1;
}

#sidebar .sidebar-header {
    padding: 20px;
    background: #6d7fcc;
}

#sidebar ul.components {
    padding: 20px 0;
    border-bottom: 1px solid #47748b;
}

#sidebar ul p {
    color: #fff;
    padding: 10px;
}

#sidebar ul li a {
    padding: 10px;
    font-size: 1.1em;
    display: block;
}

#sidebar ul li a:hover {
    color: #7386D5;
    background: #fff;
}

#sidebar ul li.active>a,
a[aria-expanded="true"] {
    color: #fff;
    background: #6d7fcc;
}

a[data-toggle="collapse"] {
    position: relative;
}

.dropdown-toggle::after {
    display: block;
    position: absolute;
    top: 50%;
    right: 20px;
    transform: translateY(-50%);
}

ul ul a {
    font-size: 0.9em !important;
    padding-left: 30px !important;
    background: #6d7fcc;
}

ul.CTAs {
    padding: 20px;
}

ul.CTAs a {
    text-align: center;
    font-size: 0.9em !important;
    display: block;
    border-radius: 5px;
    margin-bottom: 5px;
}

a.download {
    background: #fff;
    color: #7386D5;
}

a.article,
a.article:hover {
    background: #6d7fcc !important;
    color: #fff !important;
}

  </style>

<script type="text/javascript">
    $(document).ready(function () {
        $("#sidebar").mCustomScrollbar({
            theme: "minimal"
        });

        $('#dismiss, .overlay').on('click', function () {
            // hide sidebar
            $('#sidebar').removeClass('active');
            // hide overlay
            $('.overlay').removeClass('active');
        });

        $('#sidebarCollapse').on('click', function () {
            // open sidebar
            $('#sidebar').addClass('active');
            // fade in the overlay
            $('.overlay').addClass('active');
            $('.collapse.in').toggleClass('in');
            $('a[aria-expanded=true]').attr('aria-expanded', 'false');
        });
    });
</script>
</head>



<body>
<!-- container class is for bootstrap's layout. Wrapper is for the side bar navigation -->
<div class="container wrapper">
      <nav id="sidebar">

        <div id="dismiss">
            <i class="fas fa-arrow-left"></i>
        </div>

        <div class="sidebar-header">
            <h3>Navigation</h3>
        </div>

        <ul class="list-unstyled components">
            <p>Sections</p>
            <li class="active">
                <a href="#">Top of Page</a>
            </li>
            <li>
                <a href="#short_variants_gene_section">Short Variants</a>
            </li>
            <li>
                <a href="#cnv_gene_section">Copy Number Variations</a>
            </li>
            <li>
                <a href="#metap_gene_section">MetaP Fisher Method</a>
            </li>
        </ul>
    </nav>

<div id="content">

        <nav class="navbar navbar-expand-lg navbar-light bg-light">
            <div class="container-fluid">

                <button type="button" id="sidebarCollapse" class="btn btn-info">
                    <i class="fas fa-align-left"></i>
                    <span>Toggle Sidebar</span>
                </button>
            </div>
        </nav>
<!-- Begins the first text in the page -->
<h1><center><a id="top">Report</a></center></h1>
<h2><center>Structural Variants Data</center></h2>
<hr></hr>



<hr></hr>
<br>
<br>
<div class="row">
  <div class="col-md-6">
    <h2>Project Summaries</h2>
      <p> '''


## Read the summary in the input
summary = open(args.project_summary, "r").read()


## Add into the HEAD variable
## Add the project summary.
head += summary


##continue to add information to HEAD after the summary was added
head += '''      </p>
  </div>
</div>


<hr></hr>
<br>
<br>

<h2>Samples included in the study</h2>
<div>
<table border="1" class="display table table-hover">
  <thead class=thead-dark>
    <tr style="text-align: center;" >
      <th>Sample_ID</th>
      <th>Congenital_Heart_Defects</th>
    </tr>
  </thead>
  <tbody align=center >'''


## Call on the "get_html_linked_table" function to make a table out of the Database
## This table is special because it has links in it. Linking entries to another webpage
head += get_html_linked_table('select * from id_table', args.database)


## Add onto the HEAD variable after the table above was made
head += '''</tbody>
</table>
</div>



<br>
<br>
<section>
<h2><a id="short_variants_gene_section"><center>Short Variants Results</center></a></h2>
</section>
<br>
<h4>Genetic Load Analysis Results</h4>
<div>
<table border="1" class="display table table-hover">
  <thead class=thead-dark>
    <tr style="text-align: center;">
      <th>LOF_gene</th>
      <th>pval</th>
      <th>stat</th>
      <th>cases_count</th>
      <th>control_count</th>
      <th>ratio</th>
    </tr>
  </thead>
  <tbody align=center>'''

#Make an NORMAL, no links, table and add it to the file
head += get_html_table('SELECT * FROM short_variant_gene2gene WHERE ratio>1 AND pval<0.05', args.database)

## Add more to the HEAD variable
head += '''  </tbody>
</table>
</div>


<br>
<br>
<h4>GO Term Functional Analysis Results</h4>

<div>
<table border="1" class="display table table-hover">
  <thead class=thead-dark>
    <tr style="text-align: center;">
      <th>GO</th>
      <th>pval</th>
      <th>stat</th>
      <th>cases_count</th>
      <th>control_count</th>
      <th>ratio</th>
      <th>function</th>
    </tr>
  </thead>
  <tbody align=center>'''

#Make another NORMAL table
head += get_html_table('SELECT "GO","pval","stat","yes_count","no_count","ratio","name" FROM short_variant_GO WHERE ratio>1 AND pval<0.05', args.database)


## Add onto the head
head += '''  </tbody>
</table>
</div>

<br>
<br>
<h4>MP Term Functional Analysis Results</h4>

<div>
<table border="1" class="display table table-hover">
  <thead class=thead-dark>
    <tr style="text-align: center;">
      <th>MP</th>
      <th>pval</th>
      <th>stat</th>
      <th>cases_count</th>
      <th>control_count</th>
      <th>ratio</th>
      <th>function</th>
    </tr>
  </thead>
  <tbody align=center>'''


## Make table
head += get_html_table('SELECT "MP","pval","stat","yes_count","no_count","ratio","name" FROM short_variant_MP WHERE ratio>1 AND pval<0.05', args.database)


## Contionue to add more into the HEAD variable
head += '''  </tbody>
</table>
</div>

<hr></hr>







<br>
<br>
<section>
<h2><a id="cnv_gene_section"><center>Copy Number Variation Results</center></a></h2>
</section>
<br>
<h4>Genetic Load Analysis Results</h4>
<a href="#top">Back to top</a>
<div>
<table border="1" class="display table table-hover">
  <thead class=thead-dark>
    <tr style="text-align: center;">
      <th>LOF_gene</th>
      <th>pval</th>
      <th>stat</th>
      <th>cases_count</th>
      <th>control_count</th>
      <th>ratio</th>
    </tr>
  </thead>
  <tbody align=center>'''

#Add a NORMAL table
head += get_html_table('SELECT * FROM CNV_gene2gene WHERE ratio>1 AND pval<0.05', args.database)

#Add more to the HEAD variable
head += '''  </tbody>
</table>
</div>


<br>
<br>
<h4>GO Term Functional Analysis Results</h4>

<div>
<table border="1" class="display table table-hover">
  <thead class=thead-dark>
    <tr style="text-align: center;">
      <th>GO</th>
      <th>pval</th>
      <th>stat</th>
      <th>cases_count</th>
      <th>control_count</th>
      <th>ratio</th>
      <th>function</th>
    </tr>
  </thead>
  <tbody align=center>'''

## Add another NORMAL table
head += get_html_table('SELECT "GO","pval","stat","yes_count","no_count","ratio","name" FROM CNV_GO WHERE ratio>1 AND pval<0.05', args.database)


## Continue to add onto the HEAD variable
head += '''  </tbody>
</table>
</div>


<br>
<br>
<h4>MP Term Functional Analysis Results</h4>

<div>
<table border="1" class="display table table-hover">
  <thead class=thead-dark>
    <tr style="text-align: center;">
      <th>MP</th>
      <th>pval</th>
      <th>stat</th>
      <th>cases_count</th>
      <th>control_count</th>
      <th>ratio</th>
      <th>function</th>
    </tr>
  </thead>
  <tbody align=center>'''

## Add another NORMAL table
head += get_html_table('SELECT "MP","pval","stat","yes_count","no_count","ratio","name" FROM CNV_MP WHERE ratio>1 AND pval<0.05', args.database)

## Continue to add onto the HEAD variable
head += '''  </tbody>
</table>
</div>





<br>
<br>
<section>
<h2><a id="metap_gene_section"><center>MetaP Fisher Results (Short Variants +  Copy Number Variation)</center></a></h2>
</section>
<br>
<h4>Genetic Load short variants + CNVs Results</h4>

<a href="#top">Back to top</a>
<div>
<table border="1" class="display table table-hover">
  <thead class=thead-dark>
    <tr style="text-align: center;">
      <th>LOF_gene</th>
      <th>pval</th>
      <th>stat</th>
      <th>cases_count</th>
      <th>control_count</th>
      <th>ratio</th>
    </tr>
  </thead>
  <tbody align=center>'''

## Add another table
head += get_html_table('SELECT "gene","pval_x","stat_x","yes_count_x","no_count_x","ratio_x" FROM MetaP_gene2gene  WHERE ratio_x>1 AND ratio_y>1 AND combined_pval<0.05', args.database)


## Add more to the HEAD variable
head += ''' </tbody>
</table>
</div>

<br>
<br>
<h4>GO Term Functional Analysis Results</h4>

<div>
<table border="1" class="display table table-hover">
  <thead class=thead-dark>
    <tr style="text-align: center;">
      <th>GO</th>
      <th>pval</th>
      <th>stat</th>
      <th>cases_count</th>
      <th>control_count</th>
      <th>ratio</th>
      <th>function</th>
    </tr>
  </thead>
  <tbody align=center>'''


## Add another NORMAL table
head += get_html_table('SELECT "GO","pval_x","stat_x","yes_count_x","no_count_x","ratio_x","name_x" FROM MetaP_GO WHERE ratio_x>1 AND ratio_y>1 AND combined_pval<0.05', args.database)

## Add to the HEAD variable
head += '''  </tbody>
</table>
</div>



<br>
<br>
<h4>MP Term Functional Analysis Results</h4>

<div>
<table border="1" class="display table table-hover">
  <thead class=thead-dark>
    <tr style="text-align: center;">
      <th>MP</th>
      <th>pval</th>
      <th>stat</th>
      <th>cases_count</th>
      <th>control_count</th>
      <th>ratio</th>
      <th>function</th>
    </tr>
  </thead>
  <tbody align=center>'''

## Add another NORMAL table
head += get_html_table('SELECT "MP","pval_x","stat_x","yes_count_x","no_count_x","ratio_x","name_x" FROM MetaP_MP WHERE ratio_x>1 AND ratio_y>1 AND combined_pval<0.05', args.database)


## Finally, last bit of the HEAD variable
head += '''  </tbody>
</table>
</div>

</div>
    <!-- Dark Overlay element -->
    <div class="overlay"></div>
</div>
</body>
</html>'''

##########################################################################################
######################## Printing out the main HOMEPAGE.HTML #############################
##########################################################################################


##Print out the homepage, write it to a file
homepage_file = open(args.homepage_path, "w")
n = homepage_file.write(head)
homepage_file.close()




##########################################################################################
####### Printing out the main individual sample pages, each same gets a page #############
##########################################################################################

##Connect to the database to pull out the sample
conn = sqlite3.connect(args.database)

## Read in the "id_table" that holds the list of all sample IDs
id_table = pd.read_sql_query("SELECT * from id_table", conn)

#Make a copy of the gene column
conn.close()

##For every ID in the ID table, generate a suub html page for it. 
for i in id_table['ID']:
    get_html_subpage(i,args.database,path_for_index_file)

