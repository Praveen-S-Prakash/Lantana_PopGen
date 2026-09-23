#!/bin/bash

# Rename VCF contigs for downstream PLINK analysis
#
# PLINK requires chromosome names in the format used by the
# downstream analysis. In this workflow, numeric contig names
# are changed to the format:
#
#     1   -> contig1
#     2   -> contig2
#     ...
#
# Input:
#   - VCF file
#
# Output:
#   - Contig-name conversion table
#   - bgzip-compressed and indexed VCF
#   - VCF with renamed contigs

INPUT=input.vcf

CONTIGS=contig.txt
RENAMED_CONTIGS=convContig.txt
RENAME_TABLE=newContig.txt

COMPRESSED=input.vcf.gz
OUTPUT=output.vcf.gz

# Extract contig names from the VCF
grep -v "^##" "$INPUT" |
    cut -f1 |
    sort -u > "$CONTIGS"

# Add the "contig" prefix
sed 's/^/contig/' "$CONTIGS" > "$RENAMED_CONTIGS"

# Create old-name / new-name conversion table
paste "$CONTIGS" "$RENAMED_CONTIGS" > "$RENAME_TABLE"

# IMPORTANT:
# Check the first line of the conversion table before continuing.
# The original workflow included a manual edit of this file.

# Compress the VCF
bgzip "$INPUT"

# Index the compressed VCF
bcftools index "$COMPRESSED"

# Rename chromosomes
bcftools annotate \
    --rename-chrs "$RENAME_TABLE" \
    "$COMPRESSED" \
    -Oz \
    -o "$OUTPUT"

# Check the resulting chromosome names
zcat "$OUTPUT" |
    grep -v "^##" |
    less -S
