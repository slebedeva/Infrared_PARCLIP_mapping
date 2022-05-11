shell.executable("/bin/bash")
from os import path
from glob import glob
import sys

""" Snakemake pipeline for PARpipe, ported from Neel Mukherjee """ 

configfile: "config_parpipe.yaml"

PHIX=config["PHIX"]
BOWTIE_INDEX=config["BOWTIE_INDEX"]
SEGE_IDX=config["SEGE_IDX"]
FA=config["FA"]
THREE_PRIME_ADAPTER_SEQUENCE=config["THREE_PRIME_ADAPTER_SEQUENCE"]
EASY_CLIP_FIVE_PRIME_ADAPTER_SEQUENCE=config["EASY_CLIP_FIVE_PRIME_ADAPTER_SEQUENCE"]
NINETEEN_NUCLEOTIDE_MARKER_SEQUENCE=config["NINETEEN_NUCLEOTIDE_MARKER_SEQUENCE"]
TWENTY_FOUR_NUCLEOTIDE_MARKER_SEQUENCE=config["TWENTY_FOUR_NUCLEOTIDE_MARKER_SEQUENCE"]
FORWARD_PRIMER_SEQUENCE=config["FORWARD_PRIMER_SEQUENCE"]
REVERSE_PRIMER_SEQUENCE=config["REVERSE_PRIMER_SEQUENCE"]
CHROM_SIZES=config["CHROM_SIZES"]
STAR_GENOME=config["STAR_GENOME"]
GTF=config["GTF"]

THREADS = config["THR"] ##my addition


## extract sample names 
SAMPLE, = glob_wildcards("fastq/{sample}_R1_001.fastq.gz")

rule all:
  input:
    #### aligned reads ####
    #expand("mapped/star/{sample}_star_uniq.bam", sample=SAMPLE),
    #expand("mapped/star/{sample}_star_nosoft.bam", sample=SAMPLE),
    #expand("mapped/star/{sample}_star_soft.bam", sample=SAMPLE),
    #expand("mapped/segemehl/{sample}_sege2.sorted.bam",sample=SAMPLE),
    expand("merged/{sample}.merged_uniq.bam",sample=SAMPLE),
    expand("merged/{sample}.merged_uniq.pile",sample=SAMPLE),
    expand("mapped/bigwigs/{sample}_fwd.bw", sample=SAMPLE),
    "merged/feature_counts.txt",

## only CLIP
rule cutadapt:
    input:
      fastq="fastq/{sample2}_CLIP_R1_001.fastq.gz",
    output:
      fastq="filtered/{sample2}_CLIP_trim.fastq.gz",
    shell:
      """
     cutadapt -a {THREE_PRIME_ADAPTER_SEQUENCE} -g {EASY_CLIP_FIVE_PRIME_ADAPTER_SEQUENCE} -b {NINETEEN_NUCLEOTIDE_MARKER_SEQUENCE} -b {TWENTY_FOUR_NUCLEOTIDE_MARKER_SEQUENCE} -b {REVERSE_PRIMER_SEQUENCE} -b {FORWARD_PRIMER_SEQUENCE} -u 6 -m 21 -j {THREADS} -n 2 -o {output.fastq} {input.fastq} 
      """

## modify for removal of also the 5'adapter barcode if using multiple ones! right now implemented in the step above (-u 6) because same barcode
rule extract_umis:
    input:
      fq="filtered/{sample2}_CLIP_trim.fastq.gz",
    output:
      fq="filtered/{sample2}_CLIP_trim_umi.fastq.gz",
    shell:
      """
      umi_tools extract -I {input.fq} --extract-method=regex --bc-pattern='(?P<umi_1>.{{7}}).*(?P<umi_2>.{{4}})$' -S {output.fq} 
      """

## only input
rule cutadapt_truseq:
    input:
      R1 = "fastq/{sample2}_input_R1_001.fastq.gz",
    output:
      R1 = "filtered/{sample2}_input_trim.fastq.gz",
    shell:
      """
      cutadapt -q 20 -m 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -j {THREADS} -o {output.R1} {input.R1} 
      """

rule extract_UMI:
    input:
      R1 = "filtered/{sample2}_input_trim.fastq.gz",
    output:
      R1 = "filtered/{sample2}_input_trim_umi.fastq.gz",
    shell:
      """
      umi_tools extract -I {input.R1} --bc-pattern=NNNNNNNN --stdout={output.R1}
      """



## remove contaminants (currently phiX, maybe expand to miRNAs/snoRNAs/tRNAs/rRNA/repeats with stats)
idx_nums = [1, 2, 3, 4, "rev.1", "rev.2"]
bt_idx_path = BOWTIE_INDEX
bt_idxs = [(bt_idx_path + ".{}.ebwt".format(x)) for x in idx_nums]

rule bowtie_idx:
  input:
    phix=PHIX,
  output:
    bt_idxs,
  params:
    outdir = BOWTIE_INDEX,
  shell:
    """
    bowtie-build {input.phix} {params.outdir}
    """

rule bowtie:
  input:
    bt_idxs,
    fq="filtered/{sample}_trim_umi.fastq.gz",
  output:
    fq="filtered/{sample}_trim_umi_phix.fastq.gz",
    phix="filtered/{sample}.phix.bam",
  params:
    idx_path = bt_idx_path,
    fq_suffix = "filtered/{sample}_trim_umi_phix.fastq",
  shell:
    """
    gunzip -c {input.fq} | bowtie -p {threads} --sam --un {params.fq_suffix} {params.idx_path} - | samtools view -bhSF4 - | samtools sort - -o {output.phix} &&
    gzip {params.fq_suffix} &&
    samtools index {output.phix}
    """

rule star_idx:
  input:
    FA
  output:
    path.join(STAR_GENOME, "Genome")
  message:
    "building star index "
  shell:
    """
    STAR --runMode genomeGenerate \
      --genomeDir {STAR_GENOME}  \
      --genomeFastaFiles {FA} \
      --runThreadN {THREADS} \
      --outFileNamePrefix {STAR_GENOME} \
      --sjdbGTFfile {GTF} --sjdbOverhang 50  \
    """

rule star_align:
    input:
      fq = "filtered/{sample}_trim_umi_phix.fastq.gz",
      genome = STAR_GENOME + "/Genome" 
    output:
      bam = "mapped/star/{sample}_Aligned.sortedByCoord.out.bam",
      bai = "mapped/star/{sample}_Aligned.sortedByCoord.out.bam.bai",
    params:
      out = "mapped/star/{sample}_",
    shell:
      """
      STAR --runMode alignReads  \
        --scoreDelOpen 0 --scoreDelBase 0  scoreInsOpen -5   \
        --outFilterMismatchNoverReadLmax 0.1  \
        --seedSearchStartLmax 25  \
        --genomeDir {STAR_GENOME}  \
        --runThreadN {THREADS}  \
        --outBAMsortingThreadN {THREADS}  \
        --readFilesIn {input.fq}  \
        --outFileNamePrefix {params.out}  \
        --readFilesCommand zcat  \
        --outFilterMultimapNmax 10  \
        --outSAMattributes All  \
        --outMultimapperOrder Random  --outSAMprimaryFlag AllBestScore  \
        --outFilterType BySJout --outSJfilterReads Unique --outFilterIntronMotifs RemoveNoncanonicalUnannotated  \
        --outSAMtype BAM SortedByCoordinate  \
        --outWigType None  --quantMode GeneCounts  &&

      samtools index {output.bam} 

      """

rule dedup_star_bam:
    input:
      bam = "mapped/star/{sample}_Aligned.sortedByCoord.out.bam",
      bai = "mapped/star/{sample}_Aligned.sortedByCoord.out.bam.bai",
    output:
      bam = "mapped/star/{sample}_star_dedup.bam",
    shell:
      """
      umi_tools dedup --method="directional" -I {input.bam} -S {output.bam} &&
      samtools index {output.bam}
      """

## extract only uniquely mapping reads
rule unique_star:
    input:
      bam="mapped/star/{sample}_star_dedup.bam",
    output:
      bam="mapped/star/{sample}_star_uniq.bam",
      bai="mapped/star/{sample}_star_uniq.bam.bai",
    params:
      header="mapped/star/{sample}_header",
      uniq="mapped/star/{sample}_uniq",
      tmp="mapped/star/{sample}_uniq.sam",
    shell:
      """
      samtools view -H {input.bam} > {params.header} &&
      samtools view {input.bam} | grep 'NH:i:1\s' > {params.uniq} &&
      cat {params.header} {params.uniq} > {params.tmp} &&
      samtools view -b {params.tmp} > {output.bam} &&
      samtools index {output.bam} &&
      rm {params.header} {params.uniq} {params.tmp}
      """

## extract soft-clipped alignments from unique reads and remap those with segemehl to recover reads with diagnostic deletions
rule extract_soft:
    input:
      bam="mapped/star/{sample}_star_uniq.bam",
    output:
      soft="mapped/star/{sample}_star_soft.bam",
      nosoft="mapped/star/{sample}_star_nosoft.bam",
    params:
      header="mapped/star/{sample}_header",
      nosoft="mapped/star/{sample}_nosoft",
      soft="mapped/star/{sample}_soft",
    shell:
      """
      samtools view -H {input.bam}  > {params.header} &&
      samtools view {input.bam} | awk '$6 ~ /S/{{print}}' > {params.soft} &&
      samtools view {input.bam} | awk '$6 !~ /S/{{print}}' > {params.nosoft} &&
      cat {params.header} {params.soft} | samtools view -b > {output.soft} &&
      cat {params.header} {params.nosoft} | samtools view -b > {output.nosoft} &&
      samtools index {output.soft} &&
      samtools index {output.nosoft} &&
      rm {params.header} {params.nosoft} {params.soft}
      """

rule bam2fq:
    input:
      soft="mapped/star/{sample}_star_soft.bam",
    output:
      fq="mapped/star/{sample}_star_soft.fastq",
    shell:
      """
      samtools fastq {input.soft} > {output.fq} 
      """

rule segemehl2:
    input:
      fq="mapped/star/{sample}_star_soft.fastq",
      sege_idx=SEGE_IDX,
    output:
      sam="mapped/segemehl/{sample}_sege2.sam",
    shell:
      """
      segemehl.x -D 2 -M 3 -H 1 --briefcigar -t 4 -i {input.sege_idx} -d {FA} -q {input.fq}  > {output.sam}
      """

rule sege_idx:
    input: FA,
    output: SEGE_IDX,
    shell: "segemehl.x -x {SEGE_IDX} -d {FA}"

rule sam2bam2:
    input:
      sam="mapped/segemehl/{sample}_sege2.sam",
    output:
      bam="mapped/segemehl/{sample}_sege2.sorted.bam",
      bai = "mapped/segemehl/{sample}_sege2.sorted.bam.bai",
    params:
      tmp="mapped/segemehl/{sample}_sege2.bam",
    shell:
      """
      samtools view -b {input.sam} > {params.tmp} &&
      samtools sort {params.tmp} -o {output.bam} &&
      samtools index {output.bam} &&
      rm {input.sam} {params.tmp}
      """

## merge back non soft clipped STAR alignments and segemehl remapped soft clipped alignments
rule merge:
    input:
      sege="mapped/segemehl/{sample}_sege2.sorted.bam",
      star="mapped/star/{sample}_star_nosoft.bam",
    output:
      bam="merged/{sample}.merged.bam",
    shell:
      """
      samtools merge {output.bam} {input.star} {input.sege} &&
      samtools index {output.bam}
      """

## take unique reads again from the merged bam
rule unique:
    input:
      bam="merged/{sample}.merged.bam",
    output:
       bam="merged/{sample}.merged_uniq.bam",
    params:
      header="merged/{sample}_header",
      uniq="merged/{sample}_uniq",
      tmp="merged/{sample}_uniq.sam",
    shell:
      """
      samtools view -H {input.bam} > {params.header} &&
      samtools view {input.bam} | grep 'NH:i:1\s' > {params.uniq} &&
      cat {params.header} {params.uniq} > {params.tmp} &&
      samtools view -b {params.tmp} > {output.bam} &&
      samtools index {output.bam} &&
      rm {params.header} {params.uniq} {params.tmp}
      """

rule pileup:
  input:
    bam="merged/{sample}.merged_uniq.bam",
  output:
    pile="merged/{sample}.merged_uniq.pile",
  shell:
    """
    samtools mpileup -f {FA} {input.bam} > {output.pile}
    """


rule make_bigwigs:
  input:
    {CHROM_SIZES},
    bam="merged/{sample}.merged_uniq.bam",
  output:
    fwd = "mapped/bigwigs/{sample}_fwd.bw",
    rev = "mapped/bigwigs/{sample}_rev.bw",
    revneg = "mapped/bigwigs/{sample}_rev_neg.bw",
  params:
    tmp_for="mapped/bigwigs/{sample}_for.bedgraph",
    tmp_rev="mapped/bigwigs/{sample}_rev.bedgraph",
  shell:
    """
    genomeCoverageBed -bg -split -ibam {input.bam} -strand + -g {CHROM_SIZES} > {params.tmp_for} &&
    genomeCoverageBed -bg -split -ibam {input.bam} -strand - -g {CHROM_SIZES} > {params.tmp_rev} &&

    bedSort {params.tmp_for} {params.tmp_for} &&
    bedSort {params.tmp_rev} {params.tmp_rev} &&

    bedGraphToBigWig {params.tmp_for} {CHROM_SIZES} {output.fwd} &&
    bedGraphToBigWig {params.tmp_rev} {CHROM_SIZES} {output.rev} &&

    awk '{{$4=$4*-1; print}}' {params.tmp_rev} > {params.tmp_rev}.tmp &&
    bedGraphToBigWig {params.tmp_rev}.tmp {CHROM_SIZES} {output.revneg} &&
    rm {params.tmp_for} {params.tmp_rev} {params.tmp_rev}.tmp
    """


rule make_chrom_sizes:
  input:
    {FA},
  output:
    {CHROM_SIZES},
  shell:
    """
    samtools faidx {FA}
    """

rule featurecounts:
    input:
      bam = expand("merged/{sample}.merged_uniq.bam", sample=SAMPLE),
      gtf = {GTF},
    output:
      "merged/feature_counts.txt",
    shell:
      """
            featureCounts \
        -a {input.gtf} \
        -T {THREADS} \
        -o {output} \
        -s 1 \
        -Q 10 \
        {input.bam}
      """ 

