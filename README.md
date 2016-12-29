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
Test that it works on C. ruddii, on of the smallest known genomes.

    cd test_files

Edit the config file for C. ruddii in test_files/candidatus_carsonella_ruddii_pv_config.txt    

Set the following parameters: 

| Parameter | Description  |
|--- | --- |
| SIFT4G_PATH | Path to the executable sift4g (installed in step 1 of Requirements) |
| PROTEIN_DB | Path to the protein database (.fa or .fasta) |
| PARENT_DIR    | Path where all the output will be generated. User must have write access. |
| ORG           | Organism name, e.g. homo_sapiens, rat, etc. Should be one word and no special characters. |
| ORG_VERSION   | The final prediction database will be in PARENT_DIR/ORG_VERSION |

After the above variables have been set, make the database:

    perl make-SIFT-db-all.pl -config test_files/candidatus_carsonella_ruddii_pv_config.txt --ensembl_download 

It takes ~30 minutes for the database to be generated.

## Usage

    perl make-SIFT-db-all.pl -config <config_file> [--ensembl_download] 
    
### Creating a SIFT 4G Database based on Ensembl gene annotations  

1. Set parameters in the config file  


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
