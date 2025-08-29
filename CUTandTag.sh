#!/bin/bash

##QC
mkdir -p fastqc
for dir in /path/to/CUTandTAG/fastq/01.RawData/*; do
  sample=$(basename "$dir")
  mkdir -p "fastqc/$sample"
  fastqc -t "$(nproc)" "$dir"/*.fq.gz -o "fastqc/$sample"
done


set -euo pipefail


# Directory for your MultiQC output
outdir="multiqc_report"
mkdir -p "$outdir"


# Run MultiQC over all FastQC result folders
multiqc fastqc/ -o "$outdir"


##Trim

mkdir -p trim_and_filter
for dir in /path/to/CUTandTAG/fastq/01.RawData/*; do
  sample=$(basename "$dir")
  mkdir -p "trim_and_filter/$sample"
  /home/bappa.ghosh/bappa/projects/36_UND_Motoki/5_pipelines/scripts/trim_and_filter_PE.pl \
    -1 "$dir"/*_1.fq.gz \
    -2 "$dir"/*_2.fq.gz \
    -a 1 -b 151 -c 1 -d 151 -m 20 -q sanger \
    -o "trim_and_filter/$sample"
done


##CUTADAPT
mkdir -p cutadapt
for fq1 in trim_and_filter/*.1.trim_1_151.minQS_20.fastq; do
  # sample = filename before first dot
  sample=$(basename "$fq1" | cut -d. -f1)
  # corresponding R2 file
  fq2="trim_and_filter/${sample}.2.trim_1_151.minQS_20.fastq"


  # run paired‐end cutadapt
  cutadapt \
    -a CTGTCTCTTATACACATCT -O 5 -q 0 \
    -o cutadapt/"${sample}"_R1_cutadapt.fastq \
    -p cutadapt/"${sample}"_R2_cutadapt.fastq \
    "$fq1" "$fq2"
done

# 1. Make a directory to hold the fixed reads & reports
mkdir -p fixed


# 2. Loop over every R1‐cutadapt file and run the fixer
for fq1 in cutadapt/*_R1_cutadapt.fastq; do
  sample=$(basename "$fq1" _R1_cutadapt.fastq)
  fq2="cutadapt/${sample}_R2_cutadapt.fastq"


  perl /home/bappa.ghosh/bappa/projects/36_UND_Motoki/5_pipelines/scripts/fix_adapterTrimmed_PE.pl \
       "$fq1" \
       "$fq2" \
       "fixed/${sample}_R1_fixed.fastq" \
       "fixed/${sample}_R2_fixed.fastq" \
       151 \
       "fixed/${sample}.fix.report"
done


## Fastqc and multiqc_post trim:

mkdir -p post_qc && fastqc -t "$(nproc)" fixed/*_fixed.fastq -o post_qc
multiqc post_qc -o post_qc

## Rm temporary files

rm *minQS_20.fastq *_cutadapt.fastq *.FilterStats.txt *.fix.report


## mapping
mkdir -p aligned && \
for f in fixed/*_R1_fixed.fastq; do \
  s=$(basename "$f" _R1_fixed.fastq); \
  bowtie2 -p "$(nproc)" \
    --local --very-sensitive --no-mixed --no-discordant \
    -I 25 -X 700 \
    -x /data/shared_genomics_data/ref_genomes/mm10/bowtie2/mm10_index \
    -1 "$f" \
    -2 "fixed/${s}_R2_fixed.fastq" \
  | samtools view -@ "$(nproc)" -bS - > "aligned/${s}.bam"; \
Done


##Sorting:


 mkdir -p sorted && for bam in aligned/*.bam; do sample=$(basename "$bam" .bam); samtools sort -@ "$(nproc)" -o sorted/${sample}_sorted.bam "$bam" && samtools index sorted/${sample}_sorted.bam; done


###Deduplicate removal:


mkdir -p dedup && \
for bam in sorted/*_sorted.bam; do \
  sample=$(basename "$bam" _sorted.bam) && \
  picard AddOrReplaceReadGroups I="$bam" O="dedup/${sample}_rg.bam" RGID="$sample" RGLB="$sample" RGPL=illumina RGPU="$sample" RGSM="$sample" && \
  picard MarkDuplicates I="dedup/${sample}_rg.bam" O="dedup/${sample}_dedup.bam" M="dedup/${sample}_dup_metrics.txt" REMOVE_DUPLICATES=true; \
done


##Bedgraph files:


mkdir -p bedgraph && \
for bam in dedup/*_dedup.bam; do \
  bedtools genomecov -ibam "$bam" -bg > "bedgraph/$(basename "$bam" _dedup.bam).bedgraph"; \
done

##Bedgraph to bigwig


mkdir -p bigwig && \
for bg in bedgraph/*.bedgraph; do \
  bedGraphToBigWig "$bg" /data/shared_genomics_data/ref_genomes/mm10/mm10.chrom.sizes \
    "bigwig/$(basename "$bg" .bedgraph).bigWig"; \
done



##

ssh ftpUser@***.***.***.*** ## Ask Motoki for IP and details
Path: /var/ftp/Takaku/



UCSC Genome Browser
track type=bigWig name="bigwitg Gianna 3" db=mm10 visibility=full alwaysZero=on maxHeightPixels=40 color=0,0,153 bigDataUrl=ftp://***.***.***.***/Takaku/file_name.bigWig



