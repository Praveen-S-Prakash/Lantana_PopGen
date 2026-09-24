
#### Selection Manhattan plot - SNP-level FST

library(qqman)
library(ggplot2)
library(dplyr)

setwd("/Users/praveenp/Desktop/")

# Read FST file
fst <- read.table("fst_output_BJ_HW.weir.fst", header = TRUE)

# Rename columns
colnames(fst) <- c("CHR", "BP", "FST")

# Remove NA / negative values
fst <- fst[!is.na(fst$FST) & fst$FST >= 0, ]

# CHR must be numeric
fst$CHR <- as.numeric(as.factor(fst$CHR))

# Add SNP column
fst$SNP <- paste(fst$CHR, fst$BP, sep = "_")

# Make sure columns are correct
colnames(fst) <- c("CHR", "BP", "FST", "SNP")

# Clean data
fst <- fst[!is.na(fst$FST) & fst$FST >= 0, ]

# Ensure correct data types
fst$CHR <- as.numeric(as.factor(fst$CHR))
fst$BP <- as.numeric(fst$BP)
fst$FST <- as.numeric(fst$FST)


# Manhattan plot
manhattan(
  fst,
  chr = "CHR",
  bp = "BP",
  p = "FST",
  snp = "SNP",
  logp = FALSE,
  ylab = "FST",
  col = c("pink", "lightgreen"),
  annotatePval = 0.9,
  annotateTop = F,
  genomewideline = 0.6,
  log1p = F,
  suggestiveline = FALSE
)


# Highlight top 0.1%
threshold <- quantile(fst$FST, 0.999)

manhattan(
  fst,
  chr = "CHR",
  bp = "BP",
  snp = "SNP",
  p = "FST",
  logp = FALSE,
  ylab = "FST",
  col = c("pink", "lightgreen"),
  annotatePval = 0.9,
  annotateTop = F,
  genomewideline = FALSE,
  suggestiveline = FALSE,
  ylim = c(0, max(fst$FST))
)

abline(h = threshold, col = "red")


# Reduce point size
manhattan(
  fst,
  chr = "CHR",
  bp = "BP",
  p = "FST",
  logp = FALSE,
  cex = 0.5,
  genomewideline = FALSE,
  suggestiveline = FALSE
)

