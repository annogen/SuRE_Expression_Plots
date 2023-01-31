# AUTHOR / DATE
# Vartika Bisht; January 31, 2023

## Plot expression plots and statistics for the given region of interest using processed SuRE data.
## plot.R : plots expression plots for a given mutation and saves it as a PNG
## html.Rmd : Combines all PNG plots from plot.R and puts it into a HTML


import os
import pandas as pd

PLOT_SCRIPT = "/scratch/vartika/GitHub/HLHS_Reports/scripts/Expression_Plots/plot.R"
HTML_SCRIPT = "/scratch/vartika/GitHub/HLHS_Reports/scripts/Expression_Plots/html.Rmd"
ENV = "/scratch/vartika/GitHub/HLHS_Reports/envs"

OUTDIR = config["OUTDIR"]
SNP_FILE = config["SNPS"]["SNP_FILE"]
COLNAMES = config["COLNAMES"]
NORM_DIR = config["NORM_DIR"]
NORM_CDNA = config["NORM_CDNA"]
SHOWNAME = config["SHOWNAME"]
LIBS = list(NORM_CDNA.keys())

if SNP_FILE is None :
    SNPS = config["SNPS"]["SNP_LIST"]
else:
    SNPS = pd.read_table(SNP_FILE, header = None)[0].tolist()

def what_file_to_load(wildcards):
    chr = wildcards.snp.split("_")[0]
    return( os.path.join(NORM_DIR,wildcards.lib,"normalise."+ wildcards.lib + "." + chr +".txt.gz") )


def get_pngs():
    return(expand(os.path.join(OUTDIR,"{snp}.plot.png"), snp = SNPS) + expand(os.path.join(OUTDIR,"{snp}.stats.png"), snp = SNPS))

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
        temp(os.path.join(OUTDIR,"{snp}.{lib}.txt"))
    conda: os.path.join(ENV,"csvtool.yml")
    params:
        cols = lambda wildcards: ",".join(list(COLNAMES.values()) + NORM_CDNA[wildcards.lib]),
        pos = lambda wildcards: wildcards.snp.split("_")[1]
    shell:
        '''
            csvcut -t -c {params.cols} <(zcat {input} )  | csvformat -T | awk -v OFS=\"\\t\" -v FS=\"\\t\" \'NR==1; NR > 1 {{ if($1 == {params.pos} ){{print $0}} }}\' > {output}
        '''

rule plot:
    input:
        lambda wildcards: expand(os.path.join(OUTDIR,wildcards.snp+".{lib}.txt"),lib = LIBS)
    output:
        plotout = os.path.join(OUTDIR,"{snp}.plot.png"),
        statsout = os.path.join(OUTDIR,"{snp}.stats.png")
    params:
        cols = COLNAMES,
        cdna = NORM_CDNA,
        script = PLOT_SCRIPT,
        libs = LIBS,
        snp = "{snp}",
        jitter = 5,
        zero_exp_offset = -50,
        libname = SHOWNAME
    conda: os.path.join(ENV,"expression_plot.yml")
    script: '{params.script}'

rule combine_pngs:
    input:
        plotout = expand(os.path.join(OUTDIR,"{snp}.plot.png"), snp = SNPS),
        statsout = expand(os.path.join(OUTDIR,"{snp}.stats.png"), snp = SNPS)
    output:
        os.path.join(OUTDIR,"Final.html")        
    conda: os.path.join(ENV,"expression_plot.yml")
    params:
        script = HTML_SCRIPT
    script: '{params.script}'