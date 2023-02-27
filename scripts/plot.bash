#!/bin/bash
# AUTHOR / DATE
# Vartika Bisht; January 31, 2023

## Bash scipt for plotting fragment overlapping mutation of interest.

# Example Bash command
# bash plot.bash -C chr12 -p 128797634 -n "/HDD/data2/scratch/vartika/projects/VB230109_HLHS_combined/Combination_of_1/SuRE_X68/NORMALISE/SuRE_X68/normalise.SuRE_X68.chr12.txt.gz /HDD/data2/scratch/vartika/projects/VB230109_HLHS_combined/Combination_of_1/SuRE_X59/NORMALISE/SuRE_X59/normalise.SuRE_X59.chr12.txt.gz /HDD/data2/scratch/vartika/projects/VB230109_HLHS_combined/Combination_of_1/SuRE_X67/NORMALISE/SuRE_X67/normalise.SuRE_X67.chr12.txt.gz" -o "/HDD/data2/scratch/vartika/projects/VB230109_HLHS_combined/Combination_of_3/X68_X59_X67/Expression_plot" -R "/HDD/data2/scratch/vartika/GitHub/SuRE_Expression_Plots/scripts/plot.R"

# PARSE OPTIONS
# Reset in case getopts has been used previously in the shell.
OPTIND=1         

SCRIPTNAME="plot.bash"


# Defaults
snpposcol="SNP_ABS_POS_hg19"
snpvarcol="SNP_VAR"
startcol="start_hg19"
endcol="end_hg19"
libcol="Lib"
bccol="BC"
TMP=$TMP

# Usage 
USAGE=
usage() {
  echo >&2 "usage: ${SCRIPTNAME} -s:v:st:en:l:t:C:p:n:o:R:"
  echo >&2 "OPTIONS:"
  echo >&2 "  -s: column name for position of SNPs in the normlised file"
  echo >&2 "  -v: column name for variable defining the mutation type (haplotype 1, haplotype 2, unreas etc) in the normlised file"
  echo >&2 "  -st: column name for start of a fragment in the normlised file"
  echo >&2 "  -en: column name for end of a fragment in the normlised file"  
  echo >&2 "  -l: column name for library name in the normlised file"  
  echo >&2 "  -t: temp directory"  
  echo >&2 "  -C: Chromosome of interest [required]"
  echo >&2 "  -p: Position of interest [required]"
  echo >&2 "  -n: gzipped normalised file [required]"
  echo >&2 "  -o: Output directory [required]"
  echo >&2 "  -R: R script for plotting [required]"  
  echo >&2 "  -h: print this message"
  echo >&2 ""
  exit 1;
}
# Opts
while getopts "?:h:s:v:st:en:l:t:C:p:n:o:R:" opt; do
  case $opt in
    s)
      snpposcol=$OPTARG;
      ;;
    v)
      snpvarcol=$OPTARG;
      ;;
    st)
      startcol=$OPTARG;
      ;;
    en)
      endcol=$OPTARG;
      ;;
    l)
      libcol=$OPTARG;
      ;;
    t)
      TMP=$OPTARG;
      ;;    
    C)
      chromosome=$OPTARG;
      ;;
    p)
      position=$OPTARG;
      ;;
    n)
      normfiles=$OPTARG;
      ;;
    o)
      outdir=$OPTARG;
      ;;
    R)
      RSCRIPTNAME=$OPTARG;
      ;;
    h)
      usage;
      ;;
    \?)
      echo "option not recognized: "$opt
      usage
      ;;
  esac
done
shift $(( OPTIND - 1 ))


# Check if output dir exists, if not, then create it
if [ -d $outdir ] 
then
    echo "$outdir directory exists." 
else
    mkdir $outdir
    echo "$outdir directory does not exists, now it is created."
fi


echo "Subset normlised count files to only contain the SNP of interest $chromosome $position"

for norm in $normfiles
do
  echo "Reading and subsetting $norm ..."
  # Select the columns needed
  # Column refering to the position of the mutation
  snpposcolnum=$(head -1  <(zcat $norm) | sed 's/\t/\n/g' | nl | grep $snpposcol | awk '{print $1}')
  snpvarcolnum=$(head -1  <(zcat $norm) | sed 's/\t/\n/g' | nl | grep $snpvarcol | awk '{print $1}')
  startcolnum=$(head -1  <(zcat $norm) | sed 's/\t/\n/g' | nl | grep $startcol | awk '{print $1}')
  endcolnum=$(head -1  <(zcat $norm) | sed 's/\t/\n/g' | nl | grep $endcol | awk '{print $1}')
  libcolnum=$(head -1  <(zcat $norm) | sed 's/\t/\n/g' | nl | grep $libcol | awk '{print $1}')
  bccolnum=$(head -1  <(zcat $norm) | sed 's/\t/\n/g' | nl | grep $bccol | awk '{print $1}')
  normcdnacol=$(head -1  <(zcat $norm) | sed 's/\t/\n/g' | nl | grep "ipcr.norm.sum" | awk '{print $1}' | paste -s -d, -)
  
  # Add 
  paste -d"\t" <(zcat $norm | tail -n +2) <(cut -d$'\t' -f $normcdnacol <(zcat $norm | tail -n +2) | awk -v OFS="\t" -v FS="\t" '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print sum}' - ) > $TMP/tmp.norm.txt
  
  awk -v OFS="\t" -v FS="\t" -v position=$position -v snpposcolnum=$snpposcolnum -v snpvarcolnum=$snpvarcolnum -v startcolnum=$startcolnum -v endcolnum=$endcolnum -v libcolnum=$libcolnum -v bccolnum=$bccolnum '{if($snpposcolnum == position ){print $bccolnum"\t"$snpposcolnum"\t"$startcolnum"\t"$endcolnum"\t"$snpvarcolnum"\t"$libcolnum"\t"$NF }}' $TMP/tmp.norm.txt | sed '1ibc\tposition\tstart\tend\tsnpvar\tlibrary\texpression' - | gzip -c > $outdir/${position}.$(basename $norm ) ;

done

# Concatenate subsetted files ( with only fragments overlapping with the mutation )
cat <(zcat $outdir/${position}.$(basename $norm)  | head -n 1) <(zcat $outdir/${position}.*.txt.gz | grep -v position) > $outdir/$chromosome.$position.txt

echo "Combined normalised count file with all fragments overlapping the mutation of interest throughout all libraries written to $outdir/${chromosome}.${position}.txt"

# Plot using the R plot scipt
echo "Plotting.."
Rscript $RSCRIPTNAME -f ${outdir}/${chromosome}.${position}.txt -o $outdir -p $position -c $chromosome
echo "Done."