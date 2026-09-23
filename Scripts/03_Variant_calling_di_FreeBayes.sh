#!/bin/bash

# Variant calling with FreeBayes under a diploid genotype model
#
# Input:
#   - Reference genome
#   - List of BAM files
#
# Output:
#   - VCF containing called variants
#
# The FreeBayes ploidy (-p) is explicitly set to 4 in the original
# analysis workflow documented in the PhD notes. This script preserves
# that command exactly.

REFERENCE=lc_cwa_10x.fasta
BAM_LIST=bam_tetra_list.txt
OUTPUT=ddRAD_freebayse_all.vcf

freebayes \
    -p 4 \
    -f "$REFERENCE" \
    --genotype-qualities \
    -L "$BAM_LIST" \
    > "$OUTPUT"
