# Heterozygosity analysis using adegenet
#
# This script converts genotype data into an adegenet genind object
# and calculates heterozygosity-related statistics.
#
# Population assignments are taken from the population column in
# groupswise$population.
#
# Input:
#   - datx: formatted genotype data prepared in a previous step
#   - groupswise: population information
#   - whole_genome_sample_name.csv
#
# Output:
#   - Heterozygosity results printed to the R console


# ============================================================
# 1. Load packages
# ============================================================

library(adegenet)


# ============================================================
# 2. Convert genotype data to a genind object
# ============================================================

# Convert the formatted genotype data to an adegenet genind object.
#
# sep = ":"       genotype allele separator
# NA.char = "9"   missing-data character
# ploidy = 2      diploid data
# type = "codom"  codominant markers
# pop             population assignments

file_genepop <- df2genind(
    data.frame(datx),
    sep = ":",
    NA.char = "9",
    ploidy = 2,
    type = "codom",
    pop = groupswise$population
)


# ============================================================
# 3. Calculate heterozygosity
# ============================================================

# Read whole-genome sample/population information.
populations <- read.csv(
    "whole_genome_sample_name.csv"
)

# Calculate heterozygosity using Hs().
#
# NOTE:
# The original script calculated Hs() twice. The first calculation
# supplied the 'populations' object, while the second did not.
# The second calculation overwrites the first result.

het_results <- Hs(
    file_genepop,
    populations,
    diploid = TRUE
)

het_results <- Hs(
    file_genepop,
    diploid = TRUE
)


# ============================================================
# 4. Print heterozygosity results
# ============================================================

print(het_results)


# ============================================================
# 5. Calculate heterozygosity scores
# ============================================================

het_results2 <- hetscore(
    file_genepop
)
