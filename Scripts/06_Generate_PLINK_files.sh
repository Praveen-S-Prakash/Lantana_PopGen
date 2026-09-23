#!/bin/bash

# Convert a VCF file to PLINK format and perform PCA
#
# Input:
#   - bgzip-compressed VCF with contig names formatted for PLINK
#
# Output:
#   - PLINK binary files
#   - PLINK text-format files
#   - PCA results
#
# --double-id assigns the same value to both family and individual IDs.

INPUT=input.vcf.gz
OUTPUT=out_file

plink \
    --vcf "$INPUT" \
    --double-id \
    --aec \
    --pca \
    --make-bed \
    --recode \
    --out "$OUTPUT"
