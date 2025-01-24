## SuRE Expression Plots
This repository has scirpts that will help you plot expression plots for given processed SuRE data.

### Sub-directories
* envs - Directory to store all conda environments required to run the pipeline. 
* scripts - Directory to store all associated scripts.

### Using bash script for plotting all fragments overlapping the mutation of interest.
```
bash plot.bash [OPTIONS]
OPTIONS:
  -s: column name for position of SNPs in the normlised file [required]
  -v: column name for variable defining the mutation type (REF=0, ALT=1, ALT=2 etc) in the normlised file [required]
  -st: column name for start of a fragment in the normlised file [required]
  -en: column name for end of a fragment in the normlised file [required]
  -l: column name for library name in the normlised file
  -t: temp directory
  -C: Chromosome of interest [required]
  -p: Position of interest [required]
  -n: gzipped normalised file [required]
  -o: Output directory [required]
  -R: R script for plotting [required]
  -h: print this message
```
#### Using plot.R
```
Usage: Rscript plot.R [options]

Options:
        -f NORMFILE, --normfile=NORMFILE
                Path to the combined and subsetted normalised file [required]

        -o OUTDIR, --outdir=OUTDIR
                Output directory [required]

        -p POSITION, --position=POSITION
                Position of the mutation of interest [required]

        -c CHROMOSOME, --chromosome=CHROMOSOME
                Chromosome of interest [required]

        -j JITTER, --jitter=JITTER
                itter between the fragments

        -z ZEROOFFSET, --zerooffset=ZEROOFFSET
                Ofset for fragments with zero expression, so that we can look at them seperately.

        -h, --help
                Show this help message and exit

```

```
### bash plot.bash -C chr4 -p 57119312 -n "/HDD/data2/scratch/vartika/projects/VB230109_HLHS_combined/Combination_of_6/X68_X59_X67_X57_X38_X86/INDELpipeline_Functional_Output/NORMALISE/SuRE_X38/normalise.SuRE_X38.chr4.txt.gz /HDD/data2/scratch/vartika/projects/VB230109_HLHS_combined/Combination_of_6/X68_X59_X67_X57_X38_X86/INDELpipeline_Functional_Output/NORMALISE/SuRE_X57/normalise.SuRE_X57.chr4.txt.gz /HDD/data2/scratch/vartika/projects/VB230109_HLHS_combined/Combination_of_6/X68_X59_X67_X57_X38_X86/INDELpipeline_Functional_Output/NORMALISE/SuRE_X59/normalise.SuRE_X59.chr4.txt.gz /HDD/data2/scratch/vartika/projects/VB230109_HLHS_combined/Combination_of_6/X68_X59_X67_X57_X38_X86/INDELpipeline_Functional_Output/NORMALISE/SuRE_X67/normalise.SuRE_X67.chr4.txt.gz /HDD/data2/scratch/vartika/projects/VB230109_HLHS_combined/Combination_of_6/X68_X59_X67_X57_X38_X86/INDELpipeline_Functional_Output/NORMALISE/SuRE_X68/normalise.SuRE_X68.chr4.txt.gz /HDD/data2/scratch/vartika/projects/VB230109_HLHS_combined/Combination_of_6/X68_X59_X67_X57_X38_X86/INDELpipeline_Functional_Output/NORMALISE/SuRE_X86/normalise.SuRE_X86.chr4.txt.gz" -o "/HDD/data2/scratch/vartika/projects/VB240620_JA3_OFT_figures/fig1/SuREOverview" -R "/HDD/data2/scratch/vartika/GitHub/SuRE_Expression_Plots/scripts/plot.R"
```
