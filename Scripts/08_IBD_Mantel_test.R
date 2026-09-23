# Isolation-by-Distance analysis using Mantel tests
#
# This script prepares genotype data for an IBD analysis, calculates
# pairwise genetic and geographic distance matrices, and tests the
# relationship between genetic and geographic distance using Mantel tests.
#
# Main steps:
#   1. Read and format genotype data
#   2. Convert genotype data to an adegenet genind object
#   3. Calculate pairwise genetic distances
#   4. Read geographic coordinates and calculate geographic distances
#   5. Perform Mantel tests using vegan and ade4/adegenet methods
#
# Input files:
#   - input.csv
#   - location.csv
#   - Mantal_UTM_location.csv
#
# Main outputs:
#   - proportion_of_shared_allele_matrix.csv
#   - Mantel test results
#   - Mantel test plot


# ============================================================
# 1. Read and prepare genotype data
# ============================================================

data <- read.csv("input.csv", header = TRUE)

# Truncate entries in all genotype columns to the first three characters.
# The first column is assumed to contain sample/individual information.
for (i in 2:length(data[1, ])) {
    data[, i] <- substr(data[, i], 1, 3)
}

# Write the modified genotype data back to input.csv
write.table(
    data.frame(data),
    "input.csv",
    append = FALSE,
    sep = ",",
    eol = "\n",
    row.names = FALSE,
    col.names = TRUE
)


# ============================================================
# 2. Replace "/" with ":" in genotype data
# ============================================================

datp <- read.csv("input.csv", header = TRUE)

# Replace "/" with ":" in genotype columns.
# This converts genotype separators to the format expected by
# df2genind() below.
for (i in 2:length(datp[1, ])) {
    datp[, i] <- gsub("/", ":", datp[, i])
}


# ============================================================
# 3. Transpose genotype data
# ============================================================

datt <- t(datp)
datt <- as.data.frame(datt)

# Remove the first row after transposing.
datt <- datt[-1, ]


# ============================================================
# 4. Format missing genotype data
# ============================================================

# Missing data represented by "." are converted to "9".
# "9" is subsequently specified as the missing-data character
# when creating the genind object.

datx <- datt

for (i in 1:length(datt[1, ])) {
    datx[, i] <- gsub("[.]", "9", datt[, i])
}


# ============================================================
# 5. Read population information
# ============================================================

groupswise <- read.csv("location.csv")

groups <- as.factor(groupswise$population)


# ============================================================
# 6. Convert genotype data to a genind object
# ============================================================

# Convert the formatted genotype data into an adegenet genind object.
#
# sep = ":"       genotype allele separator
# NA.char = "9"   missing-data character
# ploidy = 2      diploid data
# type = "codom"  codominant markers

file_genepop <- df2genind(
    data.frame(datx),
    sep = ":",
    NA.char = "9",
    ploidy = 2,
    type = "codom"
)


# ============================================================
# 7. Calculate genetic distance
# ============================================================

# Calculate the proportion of shared alleles between individuals.

ps <- propShared(file_genepop)

# Convert similarity to genetic distance.
gen_dist <- 1 - ps

# Save the genetic distance matrix.
write.csv(
    data.frame(gen_dist),
    "proportion_of_shared_allele_matrix.csv",
    eol = "\n",
    col.names = FALSE,
    row.names = FALSE
)


# ============================================================
# 8. Read geographic coordinates
# ============================================================

locations <- read.csv(
    "Mantal_UTM_location.csv",
    header = TRUE
)

longlat <- as.matrix(
    cbind(locations$long, locations$lat)
)

spatialPoints <- SpatialPointsDataFrame(
    data = locations,
    coords = longlat
)

# Plot sampling locations
plot(spatialPoints)


# ============================================================
# 9. Calculate geographic distance matrix
# ============================================================

# Calculate pairwise Euclidean distances between sampling locations.

dm <- dist(
    spatialPoints@coords,
    method = "euclidean",
    upper = TRUE,
    diag = TRUE
)

eucl.dist <- as.matrix(dm)

# eucl.dist = geographic distance matrix
# gen_dist  = genetic distance matrix


# ============================================================
# 10. Mantel test using vegan
# ============================================================

mantel(
    gen_dist,
    eucl.dist,
    method = "pearson",
    permutations = 999,
    strata = NULL,
    na.rm = FALSE,
    parallel = getOption("mc.cores")
)


# ============================================================
# 11. Mantel test using ade4/adegenet
# ============================================================

library(MASS)

# Calculate genetic distances from the genind allele-frequency table.
Dgen <- dist(file_genepop$tab)

# Calculate geographic distances.
Dgeo <- dist(dm)

# Perform Mantel test
ibd <- mantel.randtest(Dgen, dm)

# Display results
ibd

# Plot Mantel test result
plot(ibd)
