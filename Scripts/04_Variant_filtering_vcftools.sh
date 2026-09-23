#!/bin/bash

########   IMP    ################
# I usually do the filtering as different steps. This is just amodel
##################################


# Filter FreeBayes VCF files using VCFtools
#
# Input:
#   - Raw VCF produced by FreeBayes
#
# Output:
#   - Filtered VCF
#
# Filters documented in the PhD analysis include:
#   - Removal of indels
#   - Allele-count filtering
#   - Minimum/maximum read depth
#   - Minimum genotype quality
#   - Minimum variant quality
#   - Minor allele count
#   - Dataset-specific missing-data threshold
#
# NOTE:
# The exact --max-missing value varied between datasets and should
# be specified explicitly for each analysis.

INPUT=input.vcf
OUTPUT=output

MAX_MISSING=0.XX

vcftools \
    --vcf "$INPUT" \
    --remove-indels \
    --min-alleles 10 \
    --max-alleles 500 \
    --minDP 2 \
    --maxDP 2 \
    --minGQ 30 \
    --minQ 30 \
    --mac 3 \
    --max-missing "$MAX_MISSING" \
    --recode \
    --out "$OUTPUT"
