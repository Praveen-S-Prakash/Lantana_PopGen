
############################################################
# Population genetic analysis
# IBD, genetic distance, NJ tree, heterozygosity and inbreeding
############################################################


############################
# 1. Load packages
############################

library(PopGenReport)
library(lattice)
library(sp)
library(raster)
library(adegenet)
library(vegan)
library(hierfstat)
library(dplyr)
library(MASS)
library(ape)
library(boot)
library(phylogram)
library(NAM)


############################
# 2. Set working directory
############################

workingDir <- "/Users/praveenp/Desktop/ddRad_result_all_data/"
setwd(workingDir)


############################
# 3. Prepare genetic data
############################

# Import genetic data CSV

data <- read.csv(
  "ddRAD_freebayse_all_condig_concat_final_biallelic_minQ30_minD_maxD500_mac3_minGQ_hwe5_14less_onlyTetra_maxM95_chrRe_mantel_T.csv",
  header = TRUE
)

View(data)
head(data)

# Truncate genotype entries to the first three characters

for (i in 2:length(data[1, ])) {
  data[, i] <- substr(data[, i], 1, 3)
}

head(data)

# Write modified genetic data

write.table(
  data.frame(data),
  "ddRAD_freebayse_all_condig_concat_final_biallelic_minQ30_minD_maxD500_mac3_minGQ_hwe5_14less_onlyTetra_maxM95_chrRe_mantel_T_modified2.csv",
  append = FALSE,
  sep = ",",
  eol = "\n",
  row.names = FALSE,
  col.names = TRUE
)


# Read modified data

datp <- read.csv(
  "ddRAD_freebayse_all_condig_concat_final_biallelic_minQ30_minD_maxD500_mac3_minGQ_hwe5_14less_onlyTetra_maxM95_chrRe_mantel_T_modified2.csv",
  header = TRUE
)

View(datp)


# Replace "/" with ":" in genotype entries

for (i in 2:length(datp[1, ])) {
  datp[, i] <- gsub("/", ":", data[, i])
}

head(datp)
class(datp)


# Transpose genetic data

datt <- t(datp)

dim(datt)
View(datt)

datt <- as.data.frame(datt)

View(datt)


# Remove first row

datt <- datt[-1, ]

View(datt)


# Replace missing data "." with "9"

datx <- datt

for (i in 1:length(datt[1, ])) {
  datx[, i] <- gsub("[.]", "9", datt[, i])
}

head(datx)
View(datx)

str(datx)
dim(datx)
class(datx)


############################
# 4. Import population information
############################

groupswise <- read.csv(
  "Mantal_UTM_location_all_tetra.csv"
)

groups <- as.factor(
  groupswise$population
)

groups


############################
# 5. Convert genetic data to genind
############################

file_genepop <- df2genind(
  data.frame(datx),
  sep = ":",
  NA.char = "9",
  ploidy = 2,
  type = "codom"
)

str(file_genepop)
View(file_genepop)


############################################################
# 6. Genetic distance
#    Proportion of shared alleles
############################################################

ps <- propShared(file_genepop)

gen_dist <- 1 - ps

head(gen_dist)
class(gen_dist)
View(gen_dist)


# Save genetic distance matrix

write.csv(
  data.frame(gen_dist),
  "proportion_of_shared_allele_matrix_ddRAD_freebayse_14less_onlyTetra_maxM95_chrRe_Mantel_T.csv",
  eol = "\n",
  col.names = FALSE,
  row.names = FALSE
)


############################################################
# 7. Geographic distance
############################################################

locations <- read.csv(
  "Mantal_UTM_location_all_tetra.csv",
  header = TRUE
)

View(locations)


# Extract longitude and latitude

longlat <- as.matrix(
  cbind(locations$long, locations$lat)
)

View(longlat)


# Convert coordinates to spatial object

spatialPoints <- SpatialPointsDataFrame(
  data = locations,
  coords = longlat
)

spatialPoints

plot(spatialPoints)


# Calculate geographic distance

dm <- dist(
  spatialPoints@coords,
  method = "euclidean",
  upper = TRUE,
  diag = TRUE
)

eucl.dist <- as.matrix(dm)

dim(eucl.dist)
tail(eucl.dist)
View(eucl.dist)


# gen_dist = genetic distance
# eucl.dist = geographic distance

View(gen_dist)
View(eucl.dist)


############################################################
# 8. Isolation by distance
#    Mantel test using vegan
############################################################

mantel(
  gen_dist,
  eucl.dist,
  method = "pearson",
  permutations = 999,
  strata = NULL,
  na.rm = FALSE,
  parallel = getOption("mc.cores")
)


############################################################
# 9. Isolation by distance
#    Mantel test using adegenet
############################################################

Dgen <- dist(file_genepop$tab)

ibd <- mantel.randtest(
  Dgen,
  dm
)

ibd

plot(ibd)


############################################################
# 10. Isolation by distance plot
############################################################

dens <- kde2d(
  dm,
  Dgen,
  n = 124
)

myPal <- colorRampPalette(
  c("white", "blue", "gold", "orange", "red")
)

png(
  "IBR_YP.png",
  units = "in",
  width = 8,
  height = 5,
  res = 300
)

plot(
  dm,
  Dgen,
  pch = 20,
  cex = 0.8,
  xlab = "Geographic Distance (in meters)",
  ylab = "Genetic Distance"
)

image(
  dens,
  col = transp(myPal(300), 0.7),
  add = TRUE
)

abline(
  lm(Dgen ~ dm),
  col = "red"
)

lines(
  loess.smooth(Dgen, dm),
  col = "red"
)

title("Isolation by distance plot")

dev.off()


############################################################
# 11. NJ tree from genetic distance
############################################################

# Direct NJ tree from genetic distance matrix

nj_tree <- nj(gen_dist)

plot(nj_tree)


# Hierarchical clustering approach

hclust_object <- hclust(
  as.dist(gen_dist),
  method = "average"
)

nj_tree <- as.phylo(
  hclust_object
)

plot(
  nj_tree,
  cex = 0.2,
  edge.width = 0.5,
  main = "Neighbor-Joining Tree",
  lwd = 0.001
)


# Root tree using BENWhite_S341 as outgroup

rooted_tree <- root(
  nj_tree,
  outgroup = "BENWhite_S341"
)

plot(
  rooted_tree,
  cex = 0.1,
  main = "Rooted Neighbor-Joining Tree"
)

nodelabels(
  pch = 16,
  col = "red",
  frame = "circle",
  adj = c(0, 2)
)


# Circular tree

plot(
  as.phylo(nj_tree),
  type = "fan",
  label.offset = 0.001,
  no.margin = TRUE,
  edge.width = 0.3,
  cex = 0.3,
  show.tip.label = TRUE
)


############################################################
# 12. Bootstrap analysis of NJ tree
############################################################

distance_new <- dist.gene(datx)

View(distance_new)


hclust_object_2 <- hclust(
  distance_new,
  method = "average"
)

nj_tree_2 <- as.phylo(
  hclust_object_2
)

plot(nj_tree_2)


num_replicates <- 10

bootstrap_trees <- boot.phylo(
  nj_tree_2,
  datx,
  FUN = function(xx) nj(dist.gene(xx)),
  B = 10,
  trees = TRUE
)


# Plot NJ tree

plot(
  nj_tree_2,
  main = "Neighbor-Joining Tree"
)


# Add bootstrap values

nodelabels(
  bootstrap_trees$edge[, 2],
  cex = 0.7,
  col = "red"
)


############################################################
# 13. Nei genetic distance
############################################################

nei_dis <- Gdist(
  datx,
  method = 1
)

View(nei_dis)


############################################################
# 14. Heterozygosity
############################################################

# Recreate genind object with population information

file_genepop <- df2genind(
  data.frame(datx),
  sep = ":",
  NA.char = "9",
  ploidy = 2,
  type = "codom",
  pop = groupswise$population
)


# Calculate heterozygosity

populations <- read.csv(
  "whole_genome_sample_name.csv"
)

het_results <- Hs(
  file_genepop,
  populations,
  diploid = TRUE
)

# Calculate Hs without population information

het_results <- Hs(
  file_genepop,
  diploid = TRUE
)

View(file_genepop)
View(het_results)

print(het_results)


# Heterozygosity score

het_results2 <- hetscore(
  file_genepop
)


############################################################
# 15. Inbreeding coefficient
############################################################

inbreeding_coefficient_estimate_noN <- inbreeding(
  file_genepop,
  res.type = c("estimate")
)

View(inbreeding_coefficient_estimate_noN)


# Plot inbreeding coefficient

plot(
  inbreeding_coefficient_estimate_noN,
  ylab = "F",
  col = "blue",
  main = "inbreeding coefficient"
)


# Calculate mean if required

Fbar <- sapply(
  inbreeding_coefficient,
  mean
)


# Histogram of inbreeding coefficient

hist(
  inbreeding_coefficient_estimate_noN,
  col = "#A5D6A7",
  xlab = "F",
  main = "inbreeding coefficient"
)


# Save inbreeding coefficient

write.csv(
  inbreeding_coefficient_estimate_noN,
  file = "/Users/praveenp/Desktop/ddRad_result_all_data/inbreeding_coefficient_estimate_noN.csv"
)


############################################################
# 16. Summary of genetic data
############################################################

data_summary <- summary(
  file_genepop
)

View(data_summary)
