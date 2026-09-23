#!/bin/bash

# Align trimmed paired-end reads to the reference genome using BWA-MEM
#
# Input:
#   - Trimmed paired-end FASTQ files
#   - Reference genome
#
# Output:
#   - Coordinate-sorted BAM file
#
# The read group information is added during alignment and is used
# by downstream variant-calling tools.

# Input files
REFERENCE=/Reference.fasta
R1=~/Sample_Trimmed_1P.fastq.gz
R2=~/Sample_Trimmed_2P.fastq.gz

# Output
OUTPUT=~/Sample_sorted.bam

# Align reads, remove unmapped reads, and sort
bwa mem \
    -M \
    -t 24 \
    -R "@RG\tID:Sample\tPL:ILLUMINA\tSM:Sample" \
    "$REFERENCE" \
    "$R1" \
    "$R2" |
    samtools view -Sb -F4 - |
    samtools sort -o "$OUTPUT"

# Index sorted BAM
samtools index "$OUTPUT"

###########################################################################

# Remove PCR duplicates from aligned reads using Picard
#
# Input:
#   - Coordinate-sorted BAM file
#
# Output:
#   - BAM file with duplicate reads removed
#   - Picard duplicate metrics file
#
# Duplicate reads are removed before downstream variant calling.

INPUT=~/Sample_sorted.bam
OUTPUT=~/Sample_sorted_unique.bam
METRICS=~/Sample_sorted_Metrics.txt

PICARD=~/picard.jar

java -jar "$PICARD" MarkDuplicates \
    I="$INPUT" \
    O="$OUTPUT" \
    REMOVE_DUPLICATES=TRUE \
    METRICS_FILE="$METRICS" \
    MAX_FILE_HANDLES=1000
