
#### Extract high-FST SNPs for annotation

setwd("/Users/praveenp/Desktop/")

# Read FST file
fst <- read.table(
  "fst_output_BJ_HW.weir.fst",
  header = TRUE
)

colnames(fst) <- c("CHR", "BP", "FST")

# Clean data
fst <- fst[!is.na(fst$FST) & fst$FST >= 0, ]

# Top 1%
threshold <- quantile(fst$FST, 0.99)

top <- fst[fst$FST >= threshold, ]

# Convert to BED format
top$START <- top$BP - 1
top$END <- top$BP

write.table(
  top[, c("CHR", "START", "END", "FST")],
  "top_snps.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
