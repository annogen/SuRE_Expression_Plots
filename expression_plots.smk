# AUTHOR / DATE
# Vartika Bisht; January 31, 2023

## Plot expression plots and statistics for the given region of interest using processed SuRE data.
## plot.R : plots expression plots for a given mutation and saves it as a PNG
## html.Rmd : Combines all PNG plots from plot.R and puts it into a HTML

import os
import git
import pandas as pd


#Establish snakefile paths.
SNAKEFILE = workflow.snakefile
REPO_DIR = os.path.dirname(SNAKEFILE)

# get the Repo object corresponding to the pipeline git repo
REPO = git.Repo(REPO_DIR)

# Set up script and environment directories which are relative to the repo 
SCRIPT_DIR = os.path.join(REPO_DIR,"scripts")
ENV_DIR = os.path.join(REPO_DIR,"envs")

# Set up plot and html script path , relative to the repo
PLOT_SCRIPT = os.path.join(SCRIPT_DIR,"plot.R")
HTML_SCRIPT = os.path.join(SCRIPT_DIR,"html.Rmd")

# Set up variables from config file 
OUTDIR = config["OUTDIR"]
MUT_FILE = config["MUTS"]["MUT_FILE"]
COLNAMES = config["COLNAMES"]
NORM_DIR = config["NORM_DIR"]
NORM_CDNA = config["NORM_CDNA"]
SHOWNAME = config["SHOWNAME"]

# Set up list of libraries
LIBS = list(NORM_CDNA.keys())

# Set up positions of interest using mutataions position list mentioned in config file or using the path of the file with position list.
if MUT_FILE is None :
    MUTS = config["MUTS"]["MUT_LIST"]
else:
    MUTS = pd.read_table(MUT_FILE, header = None)[0].tolist()

# Set up input normalised SuRE files to be loaded
def what_file_to_load(wildcards):
    chr = wildcards.mut.split("_")[0]
    return( os.path.join(NORM_DIR,wildcards.lib,"normalise."+ wildcards.lib + "." + chr +".txt.gz") )


# Output PNGS
def get_pngs():
    return(expand(os.path.join(OUTDIR,"{mut}.plot.png"), mut = MUTS) + expand(os.path.join(OUTDIR,"{mut}.stats.png"), mut = MUTS))

# Output HTMLS
def get_combines_html():
    return(os.path.join(OUTDIR,"Final.html"))

rule all:
    input:
        get_pngs(),
        get_combines_html()

# We know the muation, subset the data
rule subset_relavent_info_for_mut:
    input:
        what_file_to_load
    output:
        temp(os.path.join(OUTDIR,"{mut}.{lib}.txt"))
    conda: os.path.join(ENV_DIR,"csvtool.yml")
    params:
        cols = lambda wildcards: ",".join(list(COLNAMES.values()) + NORM_CDNA[wildcards.lib]),
        pos = lambda wildcards: wildcards.mut.split("_")[1]
    shell:
        '''
            csvcut -t -c {params.cols} <(zcat {input} )  | csvformat -T | awk -v OFS=\"\\t\" -v FS=\"\\t\" \'NR==1; NR > 1 {{ if($1 == {params.pos} ){{print $0}} }}\' > {output}
        '''

rule plot:
    input:
        lambda wildcards: expand(os.path.join(OUTDIR,wildcards.mut+".{lib}.txt"),lib = LIBS)
    output:
        plotout = os.path.join(OUTDIR,"{mut}.plot.png"),
        statsout = os.path.join(OUTDIR,"{mut}.stats.png")
    params:
        cols = COLNAMES,
        cdna = NORM_CDNA,
        script = PLOT_SCRIPT,
        libs = LIBS,
        mut = "{mut}",
        jitter = 5,
        zero_exp_offset = -50,
        libname = SHOWNAME
    conda: os.path.join(ENV_DIR,"expression_plot.yml")
    script: '{params.script}'

rule combine_pngs:
    input:
        plotout = expand(os.path.join(OUTDIR,"{mut}.plot.png"), mut = MUTS),
        statsout = expand(os.path.join(OUTDIR,"{mut}.stats.png"), mut = MUTS)
    output:
        os.path.join(OUTDIR,"Final.html")        
    conda: os.path.join(ENV_DIR,"expression_plot.yml")
    params:
        script = HTML_SCRIPT
    script: '{params.script}'