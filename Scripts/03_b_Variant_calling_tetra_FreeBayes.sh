#!/bin/bash

# Variant calling with FreeBayes under the tetraploid workflow
#
# Input:
#   - Reference genome
#   - List of BAM files
#
# Output:
#   - VCF containing called variants
#
# This command follows the FreeBayes configuration documented
# in the original PhD analysis.

REFERENCE=lc_cwa_10x.fasta
BAM_LIST=bam_tetra_list.txt
OUTPUT=ddRAD_freebayse_all.vcf

freebayes \
    -f "$REFERENCE" \
    --genotype-qualities \
    -L "$BAM_LIST" \
    > "$OUTPUT"
