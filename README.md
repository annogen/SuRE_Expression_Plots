# SuRE Expression Plots
This repository has scirpts that will help you plot expression plots for given processed SuRE data.

![chr12_128797634_plot](https://user-images.githubusercontent.com/53393505/220905951-da3884bf-e9f8-4090-a2a3-529b590948dd.png)

## Sub-directories
* envs - Directory to store all conda environments required to run the pipeline. 
* scripts - Directory to store all associated scripts.

## Using bash script for plotting all fragments overlapping the mutation of interest.
```
bash plot.bash [OPTIONS]
OPTIONS:
  -s: column name for position of SNPs in the normlised file
  -v: column name for variable defining the mutation type (haplotype 1, haplotype 2, unreas etc) in the normlised file
  -st: column name for start of a fragment in the normlised file
  -en: column name for end of a fragment in the normlised file
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
