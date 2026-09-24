
#### Window-based FST Manhattan plot

library(qqman)

setwd("/Users/praveenp/Desktop/")

# Read windowed FST file
fst_window <- read.table(
  "fst_output_BJ_HW_window10kb.windowed.weir.fst",
  header = TRUE
)

# Clean data
fst_window <- fst_window[!is.na(fst_window$WEIGHTED_FST), ]

# Rename / create required columns
fst_window$CHR <- as.numeric(as.factor(fst_window$CHROM))
fst_window$BP <- fst_window$BIN_START
fst_window$FST <- fst_window$WEIGHTED_FST

# Add SNP column
fst_window$SNP <- paste(
  fst_window$CHR,
  fst_window$BP,
  sep = "_"
)

# Plot
manhattan(
  fst_window,
  chr = "CHR",
  bp = "BP",
  p = "FST",
  snp = "SNP",
  logp = FALSE,
  ylab = "Windowed FST",
  genomewideline = FALSE,
  suggestiveline = FALSE
)
