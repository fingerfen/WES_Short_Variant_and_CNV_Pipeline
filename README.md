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

  * Create a new environment in Conda
	  * Install all of the above necessary packages
	  * **Note** - Get Scipy from conda-forge channel 
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

## DON'T USE ABSOLUTE PATH IN THE CONFIG.YAML FILE
## EVERYFILES NEEDED FOR THE PIPELINE MUST BE AT LEAST AT THE SAME OR 1 DIRECTORY BELOW THE SNAKEFILE


Once all files are in the repo and the config file is edited, then in the folder where the snakemake and config.yaml files are located, run:
```
snakemake --use-singularity -s Snakefile_wes --configfile config.yaml -d . -j -p --max-jobs-per-second 10 --latency-wait 30
```

**The options/flags above are:**
- -s : Point to the location of the Snakemake file (in this case, it is in the current directory)
- --configfile : Point to the location of the config file
- -d : Specifying working directory. The "." after -d is to show the working dir is at the current folder
- -j : Set available cores
- -p : Print shell command that will be executed
- --latency-wait: Define the number of seconds to wait for a file to show up after that file has been created

-------------------------

- After running the command, sit back, relax and wait for the pipeline to finish. Then run Flask as below:
	  * Run `python3 generate_REPORT.py --database put_sqlite_database_path_here`
	  * Go to google chrome
	  * Type in `localhost:5000` to see the webpage


# Output

- The final webpage will be hosted on local host port 5000
- Relative to where the Snakemake file is, the pipeline will make `data/interim/` folder to store the intermediate files.

# Notes

- The pipeline should take a few hours to run for roughly 300 samples
- Speed may vary depending on hardware and size of your data

  
# Discussion

Please direct any questions towards      |                           |
---------------------------------------- |:-------------------------:|
Mike Xie                                 | xiem1@email.chop.edu      |
Nhat Duong                               |     |
 
 
# Citation
- Adrian Tan, Gonçalo R. Abecasis and Hyun Min Kang. Unified Representation of Genetic Variants. Bioinformatics (2015) 31(13): 2202-2204
- Pedersen, B.S., Layer, R.M. & Quinlan, A.R. _Vcfanno_: fast, flexible annotation of genetic variants. _Genome Biol_  **17,** 118 (2016). https://doi.org/10.1186/s13059-016-0973-5
- Cingolani P, Platts A, Wang le L, et al. A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3. _Fly (Austin)_. 2012;6(2):80‐92. doi:10.4161/fly.19695
- Van der Auwera GA, Carneiro M, Hartl C, Poplin R, del Angel G, Levy-Moonshine A, Jordan T, Shakir K, Roazen D, Thibault J, Banks E, Garimella K, Altshuler D, Gabriel S, DePristo M, 2013 _CURRENT PROTOCOLS IN BIOINFORMATICS 43:11.10.1-11.10.33_

