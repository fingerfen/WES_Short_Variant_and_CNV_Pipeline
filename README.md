# WES_Short_Variant_and_CNV_Pipeline
Bioinformatics pipeline on Conda Virtual Environment for Short Variant (WES) and Copy Number Variation Analyses



# DAG - Pipeline 's Flow

[<img src="dag.svg">]()

# Program Requirements:
### Preinstalled on the system
- Anaconda - or just Conda (for virtual environment setup)
- Snakemake v5.4.4 (Or anything above that)
- Python3
- xlrd (for dealing with excel files)
- Flask
- Vt
- Vcfanno
- SnpEff
- Gatk 4.0 (For filtering of Short Variants)
- Scientific packages i.e (Pandas, Spicy, Numpy)

# Input Files Requirements:
- **.vcf**
  * Variant Call Format file
  * Contains Short Variant called of all samples
  * Columns should be: ` #CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, then individual samples`

- **.bed**
  * File that stores CNV discovery results from your desired method (SNPs array, WGS, WES)
  * SHOULD ALREADY BE FILTERED
  * Columns should be: ` chr, start, end, WES from U. Washington file ID`
- **.sqlite** 
  * A database file that has a table named `id_table` in it already
  * The table is user generated
  * Columns of the `id_table` should be: `ID, Heart_Phenotype`
  * ID is the sample id
  * Heart_Phenotype is either 0 or 1 (int). 0 is controlled group, 1 is cases group
  * This table will get accessed to separate the Controlled from Cases group
- **.fasta/.fa**
  * Reference file that you used for Short Variant Discovery
- **.fasta.fai**
  * Indexed file for the reference file


# Operation Instruction

### Assumptions
- All .vcf, .bed .fa, .fai, and .sqlite files are in subdirectory of the Snakemake and config.yaml files

### Quick To-Do list to run the pipeline

- Quick overview of running the pipeline:
  * Clone this github repository
  * Move all required file into this repository, preferably make a folder and put them all in there
  * Make sure that all files above are downstream from the Snakemake file
  * Edit the config.yaml
  * Execute the Snakemake command to run the pipeline
  * The output is in data/endpoints – relative to where the Snakefile is
  * In the main directory of the repo, or where the `generate_REPORT.py` is. 
	  * Run `python3 generate_REPORT.py --database put_sqlite_database_path_here`
	  * Go to google chrome
	  * Type in `localhost:5000` to see the webpage
