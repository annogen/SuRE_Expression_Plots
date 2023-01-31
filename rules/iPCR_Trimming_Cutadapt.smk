#Module using Cutadapt to trim adapters out of raw FASTQ reads and filter for a minimum read length. Also produces INFO files of the forward reads containing barcode sequences.
#Input is paired end FASTQ files.
#Output is trimmed FASTQs and INFO + Stats file.
#Requires Cutadapt v3.5 (Development branch on github right now, not available on bioconda yet)


def get_iPCR_fastq(wildcards):
    fwd_reads = os.path.join(config["iPCR"]["FASTQ_DIR"], config["iPCR"]["SAMPLES"][wildcards.s]["R1"])
    rev_reads = os.path.join(config["iPCR"]["FASTQ_DIR"], config["iPCR"]["SAMPLES"][wildcards.s]["R2"])
    return [fwd_reads, rev_reads]

rule iPCR_trim:
  input:
    get_iPCR_fastq
  output:
    fwd_trimmed = os.path.join(config["OUTDIR"], config["iPCR"]["OUTDIR"], "{s}","fastq1","{s}_forw.fastq.gz"),
    rev_trimmed = os.path.join(config["OUTDIR"], config["iPCR"]["OUTDIR"], "{s}","fastq1","{s}_rev.fastq.gz"),
    info = os.path.join(config["OUTDIR"], config["iPCR"]["OUTDIR"], "{s}","fastq1","{s}_forw.info"),
    stats = os.path.join(config["OUTDIR"], config["iPCR"]["OUTDIR"], "{s}","fastq1","{s}.stats")
  params:
    cutmotif=lambda wildcards: config["iPCR"]["SAMPLES"][wildcards.s]["CUTMOTIF"],
    forwAdaptr=config["ADPTR_IPCR_FORW_SEQ"],
    revAdaptr=config["ADPTR_IPCR_REV_SEQ"],
    min_read_length = 5 #Hard coded in ludos iPCR trim.
  # benchmark:
  #   "/scratch/james/projects/General/JB210114_Multi_core_trimming/iPCR_trim_benchmark_results/modular/multi_core/JB210420_{s}_multi-core-ipcr-trim.txt"
  conda: os.path.join("envs", "iPCR_trimming.yml")
  threads:
    20
  shell:
     "cutadapt -j {threads} -g {params.forwAdaptr} -G {params.revAdaptr} --pair-filter=any -m {params.min_read_length} --info-file={output.info} "
     "-o {output.fwd_trimmed} -p {output.rev_trimmed} -O4 --discard-untrimmed {input} >> {output.stats}"
