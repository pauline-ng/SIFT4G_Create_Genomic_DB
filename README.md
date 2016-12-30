# SIFT4G_Create_Genomic_DB
Create genomic databases with SIFT predictions. Input is an organism's genomic DNA (.fa) file and the gene annotation file (.gtf). Output will be a database that can be used with SIFT4G_Annotator.jar to annotate VCF files.

## Requirements

1. [SIFT 4G Algorithm](../../../../rvaser/sift4g)
2. Perl  
  *DBI  
  *[Bioperl](http://www.bioperl.org/wiki/Installing_BioPerl_on_Ubuntu_Server) for running DB::Fasta  
  *LWP  

## Installation
    git clone https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB.git scripts_to_build_SIFT_db

## Test Installation

### C.ruddii example
Test that it works on C. ruddii, on of the smallest known genomes. Files will be automatically downloaded from Ensembl.

In __test_files/candidatus_carsonella_ruddii_pv_config.txt__ , set *\<PARENT_DIR\>, \<SIFT4G_PATH\>, \<PROTEIN_DB\>*     

Then make the database:

    perl make-SIFT-db-all.pl -config test_files/candidatus_carsonella_ruddii_pv_config.txt --ensembl_download 

It takes ~30 minutes for this database to be generated.

### Partial Homo sapiens example (~ 1 hour for MT and chr 21)

This example uses local files to build the database. The structure and placement of files is important.

In __test_files/homo_sapiens-test.txt__ , set *\<SIFT4G_PATH\>* and *\<PROTEIN_DB\>* 

Partial genomic DNA, gene annotation, and dbSNP annotations are provided in test_files/homo_sapiens_small

Testing human (chr21 and MT only, ~2 hrs):

    perl make-SIFT-db-all.pl -config test_files/homo-sapiens-test.txt
    
## Usage

    perl make-SIFT-db-all.pl -config <config_file> [--ensembl_download] 
    
### Making a SIFT database from local genomic and gene annotation file (.gtf)

1. Make a config file. Use __test_files/homo_sapiens-test.txt__ for a template

 a. Set *\<PARENT_DIR\>, \<ORG\>, \<ORG_VERSION\>, \<SIFT4G_PATH\>, \<PROTEIN_DB\>*
    Optional: *\<DBSNP_VCF_FILE\>*
   
    See config for details.

2. Put the genomic fasta files and the gene annotation files in their proper place:

   a. Make the folders:
   
    `mkdir <PARENT_DIR>
    mkdir <PARENT_DIR>/gene-annotation-src
    mkdir <PARENT_DIR>/chr-src
    mkdir <PARENT_DIR>/dbSNP`

   b. Put files in folders:
   
    Put genomic .fa.gz files in *\<PARENT_DIR\>chr-src*
    Put compressed gene annotation file (gtf.gz) in *\<PARENT_DIR\>/gene-annotation-src*

    Optional: Put compressed dbSNP VCF file in *\<PARENT_DIR\>/dbSNP*
    Optional: Put compressed protein file (.pep.all.fa.gz) in *\<PARENT_DIR\>/gene-annotation-src*  (used for checking)

    Example of the file structure can be found in *test_files/homo_sapiens_small*
    
3. Run command: 

    perl make-SIFT-db-all.pl -config <config_file>
    
    
    
### Creating a SIFT 4G Database based on Ensembl gene annotations  

1. Set parameters in the config file
   Use __test_files/candidatus_carsonella_ruddii_pv_config.txt__ as a template.

  a. Links to Ensembl genome and gene annotation files:*GENE_DOWNLOAD_SITE, PEP_FILE, CHR_DOWNLOAD_SITE, DBSNP_ORGANISM_DOWNLOAD_SITE (optional)*    
  
  b. Output  & SIFT4G Path: *SIFT4G_PATH, PROTEIN_DB, PARENT_DIR, ORG, ORG_VERSION*  
  
  c. Set genetic codes in *_GENETIC_CODE*  
  
2. Create the database:  

    `perl make-SIFT-db-all.pl -config <config file> --ensembl_download`
    
3. Check and use the SIFT database (directions below.)


## Check the Database

The database is stored in \<PARENT_DIR\>/\<ORG_VERSION\>

    cd <PARENT_DIR>/<ORG_VERSION>
SIFT4G_DB_STRUCTURE
The SIFT database structure is decribed [here](http://sift-dna.org/sift4g/AnnotateVariants.html#SIFT4G_DB_STRUCTURE)

## Annotate VCF files with the SIFT Database

    java -jar <Path to SIFT4G_Annotator> -c -i <Path to input vcf file> -d <Path to SIFT4G database directory> -r <Path to your results folder> -t 

Complete instructions [here](http://sift-dna.org/sift4g/AnnotateVariants.html)

## Configuration File

| Parameter | Description  |
|--- | --- |
| SIFT4G_PATH | Path to the executable sift4g (installed in step 1 of Requirements) |
| PROTEIN_DB | Path to the protein database (.fa or .fasta) Recommend UniProt 90 |
| PARENT_DIR    | Path where all the output will be generated. User must have write access. |
| ORG           | Organism name, e.g. homo_sapiens, rat, etc. Should be one word and no special characters. |
| ORG_VERSION   | The final prediction database will be in *\<PARENT_DIR\>/\<ORG_VERSION\>* |
|  GENE_DOWNLOAD_SITE (Ensembl only) | URL to download gene annotation (gtf) file |
|  PEP_FILE (Ensembl only) | 	URL to download the protein sequence file. This file is optional and used for quality control. |
| CHR_DOWNLOAD_SITE (Ensembl only) | URL to download the genome sequence (.fa) |
| DBSNP_ORGANISM_DOWNLOAD_SITE (Ensembl only) |URL folder to download dbSNP annotations (optional) |  
| DBSNP_VCF_FILE *(optional)* | A \*.vcf.gz file that will be used to annotate variants with dbSNP rs id's |
| GENETIC_CODE_TABLE | Genetic code to be used to translate the DNA sequence into proteins ([integer value based on NCBI](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi))* |
| GENETIC_CODE_TABLENAME *(not used)* | String to remind user what GENETIC_CODE_TABLE is being used | 
| MITO_GENETIC_CODE_TABLE |  	Fasta sequences named Mt, chrM, Mito will use this genetic code. (This is for mitochondrial sequences). To disable or edit this feature, edit the function chr_is_mito in dna_protein_subs.pl  | 
| MITO_GENETIC_CODE_TABLENAME *(not used)* | String to remind user what MITO_GENETIC_CODE_TABLE is being used |

## Monitoring the Database Creation Process

Because the database can take hours/days to complete, here is what to look for:

| If the Terminal says .... |  Check the following: |
| --- | --- |
| `making single records file` |  `ls -lt <PARENT_DIR>/singleRecords/*` is being updated |
| `make the fasta sequences` | `ls fasta/* | wc -l` number of files should be increasing |
| \*processing database | SIFT 4G Algorithm is running |
| `populating databases` | `ls -lt <PARENT_DIR>/singleRecords_with_scores/*` is being updated |



    
    



