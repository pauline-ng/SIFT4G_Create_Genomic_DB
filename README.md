# SIFT4G_Create_Genomic_DB
Create genomic databases with SIFT predictions. Input is an organism's genomic DNA (.fa) file and the gene annotation file (.gtf). Output will be a database that can be used with SIFT4G_Annotator.jar to annotate VCF files.

*We can also build your database for you. However, we ask that the person who builds your SIFT database to be listed as a  middle author on your paper. <A HREF="http://siftdna.org/sift-bin/contact.pl">You can email us the details.</A>*

## Requirements

1. [SIFT 4G Algorithm](../../../../rvaser/sift4g)
2. Perl  
  *DBI  
  *[Bioperl](http://www.bioperl.org/) for running DB::Fasta  
  *LWP  
  *Switch.pm   (`sudo apt-get install libswitch-perl`)

## Installation
    git clone https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB.git scripts_to_build_SIFT_db

## Test Installation

### *C.ruddii* example
*C. ruddii* is a small genome and can quickly test if everything is working. The genome and gene files are automatically downloaded from Ensembl.

1. Set variables in configuration file.  
   `cd test_files/candidatus_carsonella_ruddii_pv_config.txt`
   
   Edit the config file to set *\<PARENT_DIR\>, \<SIFT4G_PATH\>, \<PROTEIN_DB\>*     
   [(See config details)](#configFile)
   
   `mkdir <PARENT_DIR>`
   
   Go back to the __scripts_to_build_SIFT_db__ directory
   
   `cd ..`
   


2. Make the database:

    `perl make-SIFT-db-all.pl -config test_files/candidatus_carsonella_ruddii_pv_config.txt --ensembl_download` 

   It takes ~30 minutes for this database to be generated in *\<PARENT_DIR\>/\<ORG_VERSION\>*.  

3. [Check the database](#checkDB)

   The database should be in a folder named something like: candidatus_carsonella_ruddii_pv/ASM1036v1.34
   

### Partial *Homo sapiens* example 

This example uses local files to build a database of human chr21 and mitochondrial genes. Do this exercise if you are building a SIFT database with a genome that is on your local computer.

Files are already provided in [scripts_to_build_SIFT_db/test_files/homo_sapiens_small](./test_files/homo_sapiens_small)
   - Genomic DNA (.fa.gz)  
   - Gene annotation (.gtf.gz)  
   - dbSNP annotations (.vcf.gz) 


1.  Set the variables in the config file.

    `cd scripts_to_build_SIFT_db/test_files/`
    
    Set variables in the config file __homo_sapiens-test.txt__:   *\<SIFT4G_PATH\>* and *\<PROTEIN_DB\>*
    
    Note that \<PARENT_DIR\> is __already__ set to ./test_files/homo_sapiens_small  
    SIFT scripts will look for the genome and gene annotation files in that folder (which are provided in this example).

2.  Make the database:

    `cd ..`  
    `perl make-SIFT-db-all.pl -config test_files/homo_sapiens-test.txt`
    
    It takes ~2 hours for human chr21 and mitochondria predictions to be generated in *\<PARENT_DIR\>/\<ORG_VERSION\>*.

3.  [Check the database](#checkDB)

## Usage

    perl make-SIFT-db-all.pl -config <config_file> [--ensembl_download] 
    
__Directions for making a SIFT database from:__
* [Ensembl download of genomic and gene annotation files](#DBfromEnsembl)  
* [local genomic and gene annotation file (.gtf)](#DBfromGTF)  
* [local genomic and gene annotation file (.gff)](#DBfromGFF)  
    
### <a name="DBfromGTF"></a>Making a SIFT database from local genomic and gene annotation file (.gtf)

1. Create a config file. 

   a.  Use __test_files/homo_sapiens-test.txt__ as a template

   `cd test_files/`  
   `cp homo_sapiens-test.txt <my_org_config.txt>`
   
   b. In \<my_org_config.txt\>, set *\<PARENT_DIR\>, \<ORG\>, \<ORG_VERSION\>, \<SIFT4G_PATH\>, \<PROTEIN_DB\>*  
      Check *GENETIC_CODE_TABLE* and *MITO_GENETIC_CODE_TABLE* is set correctly.  
      Optional to set: *\<DBSNP_VCF_FILE\>*
   
      [See config details](#configFile).

2. Put the genomic fasta files and the gene annotation files in their proper place:

   a. Make the folders:
   
    `mkdir <PARENT_DIR>`  
    `mkdir <PARENT_DIR>/gene-annotation-src`  
    `mkdir <PARENT_DIR>/chr-src`  
    `mkdir <PARENT_DIR>/dbSNP`  

   b. Move your files into the appropriate folders:
   
    Put compressed genomic fasta files (.fa.gz) in *\<PARENT_DIR\>/chr-src*  
    
    `mv *.fa.gz <PARENT_DIR>/chr-src`
    
    Put compressed gene annotation file (gtf.gz) in *\<PARENT_DIR\>/gene-annotation-src*  

    `mv *.gtf.gz <PARENT_DIR>/gene-annotation-src`
    
    Optional: Put compressed dbSNP VCF file in *\<PARENT_DIR\>/dbSNP*  
    `mv *.vcf.gz <PARENT_DIR>/dbSNP`
    
    Optional: Put compressed protein file (.pep.all.fa.gz) in *\<PARENT_DIR\>/gene-annotation-src*  
    
    `mv *.pep.all.fa.gz <PARENT_DIR>/gene-annotation-src`
    
    This file is used for checking.  

    Example of the file structure can be found in [test_files/homo_sapiens_small](./test_files/homo_sapiens_small)
    
3. Run command: 

    `perl make-SIFT-db-all.pl -config <config_file>`
    
4. [Check the database](#checkDB)

5. [Annotate a VCF file with your database](#annotate)
    
### <a name="DBfromGFF"></a>Making a SIFT database from genomic DNA (.fa.gz) and gene annotation file (.gff)

Use this if you have a gff file (like that supplied from Phytozyme)

1. Download and install [gffread](https://github.com/gpertea/gffread)

2. Convert the gene annotation .gff file to a .gtf file (because SIFT processes gtf files).
   In the gtf file, make sure the 9th column (attribute column) says **gene_biotype "protein_coding;"** for rows which are labelled as exon, CDS, stop_codon, and start_codon. 

   `mv <gff3.gz> <PARENT_DIR>/gene-annotation-src`  
   `gunzip <*.gff3.gz>`  
   `gffread <*.gff3> -T -o [FILENAME].gene.gtf` 
   
   `#below command is for Phytozyme, to ensure 'protein_coding' is in 2nd column`  
   `perl -pe 's/phytozomev10/protein_coding/g' [FILENAME].gene.gtf > FILENAME.mod.gtf`  
   
   `mv FILENAME.mod.gtf [FILENAME].gtf`  
   `gzip [FILENAME].gtf`   

3. Then follow instructions for [building a database using a gtf file](#DBfromGTF)
    
### <a name="DBfromEnsembl"></a>Creating a SIFT 4G Database based on Ensembl gene annotations  

This will download genome and gene files directly from Ensembl with the option __--ensembl_download__

1. Set parameters in the config file.  

   Use __test_files/candidatus_carsonella_ruddii_pv_config.txt__ as a template.

  a. Set weblinks to Ensembl genome and gene annotation files: *GENE_DOWNLOAD_SITE, PEP_FILE, CHR_DOWNLOAD_SITE*  
     Optional: *DBSNP_ORGANISM_DOWNLOAD_SITE*    
     [Config file details](#configFile)
     
  b. Set output paths: *PARENT_DIR, ORG, ORG_VERSION*  
  
  c. Set [genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) in *GENETIC_CODE_TABLE, MITO_GENETIC_CODE_TABLE* 
  
  If you're working with a vertebrate, it's
  ```GENETIC_CODE_TABLE=1
GENETIC_CODE_TABLENAME=Standard
MITO_GENETIC_CODE_TABLE=2
MITO_GENETIC_CODE_TABLENAME=Vertebrate Mitochondrial

```
  
  d. Set SIFT4G paths: *SIFT4G_PATH, PROTEIN_DB*
  
2. Create the database:  

    `perl make-SIFT-db-all.pl -config <config file> --ensembl_download`
    
3. [Check the database](#checkDB)

4. [Annotate a VCF file with your database](#annotate)

## <a name="checkDB"></a>Check the Database

The database is stored in *\<PARENT_DIR\>/\<ORG_VERSION\>* which was set in the config file.

    cd <PARENT_DIR>/<ORG_VERSION>
    more <PARENT_DIR>/<ORG_VERSION>/CHECK_GENES.LOG
    
The file CHECK_GENES.LOG contains a summary of SIFT predictions by chromosome. There are 4 columns:
  - Chr	 
  - Genes with SIFT Scores  
  - Pos with SIFT scores	 
  - Pos with Confident Scores
  
The last line summarizes predictions for the entire genome:

__ALL   # (#/#)  # (#/#) # (#/#)__  

Your database is __done__ if the percentages are high for the last 3 different columns. Woohoo!

A SIFT database is made for each chromosome in the file \<chr\>.gz  The SIFT database structure is described [here](http://sift-dna.org/sift4g/AnnotateVariants.html#SIFT4G_DB_STRUCTURE)

    zcat <chr>.gz | more   # does it look all right?
    zcat <chr>.gz | grep CDS | more # SIFT numeric scores will be in columns 10-12. If too many rows say "NA", that's a problem
    
## <a name="annotate"></a> Annotate VCF files with the SIFT Database

1. Download the SIFT 4G Annotator (a Java executable) [here](http://sift-dna.org/sift4g/AnnotateVariants.html)  

2. Commandline:
   `java -jar <Path to SIFT4G_Annotator> -c -i <Path to input vcf file> -d <Path to SIFT4G database directory> -r <Path to your results folder> -t`

Complete instructions [here](http://sift-dna.org/sift4g/AnnotateVariants.html)

---

## <a name="configFile"></a>Configuration File

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
| MITO_GENETIC_CODE_TABLE |  	Fasta sequences named Mt, chrM, or Mito will use this genetic code. (This is for mitochondrial sequences). To disable or edit this feature, edit the function chr_is_mito in dna_protein_subs.pl  | 
| MITO_GENETIC_CODE_TABLENAME *(not used)* | String to remind user what MITO_GENETIC_CODE_TABLE is being used |
| PLASTID_GENETIC_CODE_TABLE | Fasta sequences named Pt, PT, chrPt, chrPT, or chloroplast will use this genetic code. (This is for chlorplasts). To disable or edit this feature, edit the function chr_is_plastid in dna_protein_subs.pl  |


---

## Monitoring the Database Creation Process

Because the database can take hours/days to complete, here is what to look for:

| If the Terminal says .... |  Check the following: |
| --- | --- |
| `making single records file` |  `ls -lt <PARENT_DIR>/singleRecords/*` is being updated |
| `make the fasta sequences` | `ls fasta/* | wc -l` number of files should be increasing |
| \*processing database | SIFT 4G Algorithm is running |
| `populating databases` | `ls SIFT_predictions/*` Inspect the SIFT prediction files  |
|                        |   `ls -lt <PARENT_DIR>/singleRecords_with_scores/*` is being updated |



    
    



