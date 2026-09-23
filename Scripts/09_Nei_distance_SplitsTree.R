# Nei's genetic distance and SplitsTree input
#
# This script imports a VCF file, converts the genotype data to an
# adegenet genlight object using a custom function for tetraploid data,
# assigns population information from individual names, and calculates
# Nei's (1972) genetic distance.
#
# Genetic distances are calculated:
#   1. Between individuals
#   2. Between populations
#
# The resulting distance matrices are exported in PHYLIP distance
# format for downstream analysis and visualization in SplitsTree.
#
# Input:
#   - input.vcf
#   - adegenet_functions.R
#
# Output:
#   - indiv_Neis_distance_mix_ploidy.phy.dst
#   - pops_Neis_distance_mix_ploidy.phy.dst


# ============================================================
# 1. Load custom functions
# ============================================================

source("./adegenet_functions.R")


# ============================================================
# 2. Import VCF data
# ============================================================

library(vcfR)
library(adegenet)
library(StAMPP)

vcf <- read.vcfR("input.vcf")


# ============================================================
# 3. Convert VCF to genlight
# ============================================================

# Convert the VCF to a genlight object using the custom function
# defined in adegenet_functions.R.
#
# The custom function is specifically used for tetraploid data.

aa.genlight <- vcfR2genlight.tetra(vcf)


# ============================================================
# 4. Assign SNP names
# ============================================================

# Create SNP names from the VCF contig/chromosome and position.

locNames(aa.genlight) <- paste(
    vcf@fix[, 1],
    vcf@fix[, 2],
    sep = "_"
)


# ============================================================
# 5. Assign population information
# ============================================================

# Population names are derived from the first five characters
# of each individual name.

pop(aa.genlight) <- substr(
    indNames(aa.genlight),
    1,
    5
)


# ============================================================
# 6. Calculate Nei's (1972) genetic distance
#    between individuals
# ============================================================

# Calculate Nei's (1972) distance between individuals.

aa.D.ind <- stamppNeisD(
    aa.genlight,
    pop = FALSE
)

# Export the individual-level distance matrix in PHYLIP format
# for use in SplitsTree.

stamppPhylip(
    aa.D.ind,
    file = "indiv_Neis_distance_mix_ploidy.phy.dst"
)


# ============================================================
# 7. Calculate Nei's (1972) genetic distance
#    between populations
# ============================================================

# Calculate Nei's (1972) distance between populations.

aa.D.pop <- stamppNeisD(
    aa.genlight,
    pop = TRUE
)

# Export the population-level distance matrix in PHYLIP format
# for use in SplitsTree.

stamppPhylip(
    aa.D.pop,
    file = "pops_Neis_distance_mix_ploidy.phy.dst"
)
