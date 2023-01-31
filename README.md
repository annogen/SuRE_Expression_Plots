# SuRE Expression Plots
This repository has scirpts that will help you plot expression plots for given processed SuRE data.

## Sub-directories
* envs - Directory to store all conda environments required to run the pipeline. 
* scripts - Directory to store all associated scripts.
* data - Directory to store meta-data tsv sheets, sample names as text files etc.  

## Specifying positions of interest
Specify the SNPs of interest in expression_plots.yml as chr_pos or give a SNP_FILE of interest with posiitons of of the format chr_pos in new line.

## Limitations:
Only works for one cell line at a time

## Reuired Entries in the YAML
* OUTDIR - output directory
* NORM_DIR - directory with normalised SuRE Counts file
* NORM_CDNA - names of the ipcr normalised and scaled cDNA columns for each librray
* SHOWNAME - names for the libraries that you would like to be displayed on the plots
* SNPS - input position of interests as:
    * SNP_FILE - a file 
    * SNP_LIST - a list in the yaml
* COLNAMES - names of column in your dataframe corresponding to :
    * SNPabspos - absolute position of the mutation
    * start - start of the fragment
    * end - end of the fragment
    * lib - library name
    * SNPVAR - SNP VAR value assigned by the pipeline. eg: 0 - REF, 1 - ALT etc.
