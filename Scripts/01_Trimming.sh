#!/bin/bash

# Trimming paired-end Illumina reads using Trimmomatic
#
# This script removes adapter contamination and low-quality bases
# from paired-end sequencing reads.
#
# Input:
#   - Raw paired-end FASTQ files
#
# Output:
#   - Paired trimmed reads (*_1P and *_2P)
#   - Unpaired trimmed reads (*_1U and *_2U)
#
# Main filtering:
#   - TruSeq paired-end adapter removal
#   - Leading/trailing bases with quality < 3 removed
#   - Reads with average quality < 30 removed
#   - Reads shorter than 30 bp removed

# Input files
R1=/Sample_R1_001.fastq.gz
R2=/Sample_R2_001.fastq.gz

# Output files
R1_P=/Sample_Trimmed_1P.fastq.gz
R1_U=/Sample_Trimmed_1U.fastq.gz
R2_P=/Sample_Trimmed_2P.fastq.gz
R2_U=/Sample_Trimmed_2U.fastq.gz

# Trimmomatic
TRIMMOMATIC=/softwares/trimmomatic/trimmomatic.jar
ADAPTERS=/softwares/trimmomatic/adapters/TruSeq3-PE.fa

# Run Trimmomatic
java -jar "$TRIMMOMATIC" PE \
    -threads 16 \
    -phred33 \
    "$R1" \
    "$R2" \
    "$R1_P" \
    "$R1_U" \
    "$R2_U" \
    "$R2_P" \
    "ILLUMINACLIP:${ADAPTERS}:2:30:10:2" \
    KEEP_BOTH_READS \
    AVGQUAL:30 \
    LEADING:3 \
    TRAILING:3 \
    MINLEN:30
