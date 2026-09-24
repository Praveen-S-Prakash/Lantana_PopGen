
## Correlation between flower colour and genetic structure
## A CSV was created using the Q matrix at K = 11 and flower colour for each individual

workingDir <- "/Users/praveenp/Desktop/"
setwd(workingDir)

data <- read.csv("Flower_colour_Structure_Q.csv", header = TRUE)

colnames(data)

# MANOVA testing whether genetic structure differs among flower colours

manova_res <- manova(
  cbind(Cluster1, Cluster2, Cluster3, Cluster4, Cluster5,
        Cluster6, Cluster7, Cluster8, Cluster9, Cluster10, Cluster11)
  ~ Flower_colour,
  data = data
)

summary(manova_res, test = "Pillai")


# ANOVA testing the significance of each STRUCTURE cluster

for (i in 1:11) {
  cluster_name <- paste0("Cluster", i)
  model <- aov(data[[cluster_name]] ~ data$Flower_colour)

  cat("\n", cluster_name, "\n")
  print(summary(model))
}


# Tukey tests for each STRUCTURE cluster

for (i in 1:11) {
  cluster_name <- paste0("Cluster", i)

  model <- aov(data[[cluster_name]] ~ data$Flower_colour)
  tukey <- TukeyHSD(model)

  cat("\n", cluster_name, "\n")
  print(tukey)
}


## Alternative MANOVA code tested previously
## manova_res <- manova(as.matrix(data[, 3:13])[ , -11] ~ Flower_colour, data = data)
## summary(manova_res, test = "Pillai")
