
############################################################
# DAPC and PCA from SNP data
############################################################


############################
# 1. Load libraries
############################

library(vcfR)
library(adegenet)
library(scales)


############################
# 2. Set working directory
############################

workingDir <- "/Users/praveenp/Desktop/ddRad_result_all_data/"
setwd(workingDir)


############################
# 3. Import VCF and create genlight object
############################

gll <- vcfR2genlight(
  read.vcfR(
    "ddRAD_freebayse_all_condig_concat_final_biallelic_minQ30_minD_maxD500_mac3_minGQ_hwe5_14less_onlyTetra_maxM90_chrRe.vcf"
  )
)

class(gll)


############################################################
# 4. DAPC analysis
#    New analysis
############################################################

# Find genetic clusters using find.clusters
# The number of clusters is evaluated up to 50

clus <- find.clusters(
  gll,
  max.n.clus = 50
)


# Run DAPC using the clusters identified above

dp <- dapc(
  gll,
  clus$grp
)


############################
# 5. DAPC plots
############################

# Basic DAPC scatter plot

scatter(dp)


# DAPC scatter plot with modified position and symbols

scatter(
  dp,
  posi.da = "topleft",
  bg = "white",
  pch = 17:22
)


# DAPC scatter plot with custom colours

myCol <- c(
  "darkblue",
  "purple",
  "green",
  "orange",
  "red",
  "blue",
  "black"
)

scatter(
  dp,
  posi.da = "bottom",
  bg = "white",
  pch = 17:22,
  cstar = 0,
  col = myCol,
  scree.pca = TRUE,
  posi.pca = "topleft"
)


# DAPC scatter plot with larger points

scatter(
  dp,
  bg = "white",
  pch = 20,
  cell = 0,
  cstar = 0,
  col = myCol,
  solid = 0.4,
  cex = 3,
  clab = 0,
  leg = TRUE
)

legend(
  "topleft",
  legend = paste("Cluster", 1:6),
  col = myCol,
  pch = 20,
  cex = 1.2
)


# Plot samples against one DAPC axis

scatter(
  dp,
  1,
  1,
  bg = "white",
  scree.da = FALSE,
  legend = TRUE,
  solid = 0.4
)


############################
# 6. DAPC summary
############################

Sum <- summary(dp)


############################
# 7. Cluster assignment plots
############################

assignplot(
  dp,
  subset = 1:100
)


# Structure-like plot

compoplot(
  dp,
  posi = "bottomright",
  txt.leg = paste("Cluster", 1:3),
  lab = "",
  ncol = 1,
  xlab = "individuals",
  col = funky(3)
)


############################
# 8. Identify individuals in a cluster
############################

# Find individuals assigned to Cluster 1

individuals_in_cluster1 <- which(
  dp$grp == 1
)

cat(
  "Individuals in Cluster 1:",
  individuals_in_cluster1,
  "\n"
)


# Get sample names for Cluster 1

sample_names_cluster1 <- gll$ind.names[
  individuals_in_cluster1
]


############################################################
# 9. Older PCA/DAPC analysis
############################################################

# Convert genlight object to matrix

x <- as.matrix(gll)

gi <- as.genind(x)


############################
# 10. Import population information
############################

pop_info <- read.table(
  "popmap_FC_14_less_tetra_dapc.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

pop(gll) <- pop_info$population_code


############################
# 11. Principal component analysis
############################

pca <- glPca(
  gll,
  nf = 2
)


# Plot PCA

quartz()

col <- scales::hue_pal()(
  length(unique(pop(gll)))
)

s.class(
  pca$scores,
  pop(gll),
  col = col,
  axesell = FALSE,
  cstar = 0,
  grid = FALSE
)


############################################################
# 12. Find clusters
############################################################

clus <- find.clusters(
  gll,
  max.n.clus = 40
)


# Reassign clusters using known population identities

clus$grp <- pop(gll)


############################################################
# 13. DAPC using known populations
############################################################

dp <- dapc(
  gi,
  pop = pop(gll),
  n.da = 1,
  perc.pca = 80
)


############################
# 14. Plot DAPC distributions
############################

quartz()

scatter(
  dp,
  bg = "white",
  scree.da = FALSE,
  legend = TRUE,
  solid = 0.4,
  col = col
)


############################################################
# 15. Compare population distributions on PC1
############################################################

scores <- pca$scores

ldens <- tapply(
  scores[, "PC1"],
  pop(gll),
  density
)

allx <- unlist(
  lapply(
    ldens,
    function(e) e$x
  )
)

ally <- unlist(
  lapply(
    ldens,
    function(e) e$y
  )
)


quartz()

plot(
  allx,
  ally,
  type = "n",
  xlab = "Discriminant function",
  ylab = "Density",
  axes = TRUE
)


for (i in 1:length(ldens)) {
  polygon(
    c(
      ldens[[i]]$x,
      rev(ldens[[i]]$x)
    ),
    c(
      ldens[[i]]$y,
      rep(
        0,
        length(ldens[[i]]$x)
      )
    ),
    col = alpha(
      col[i],
      0.6
    ),
    lwd = 1,
    border = col[i]
  )
}


############################################################
# End
############################################################

# DAPC is used to visualize genetic differentiation among
# groups after retaining variation in the selected PCs.
