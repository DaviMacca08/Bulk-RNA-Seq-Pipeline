library(DESeq2)
library(ggplot2)

# ------------------------------------------
#            Building DESeq Object
# ------------------------------------------

# 01-Read counts.txt from featureCounts

counts <- read.table("counts.txt", sep = "\t", header = T, skip = 1)
counts <- counts[, -c(2:6)]
rownames(counts) <- counts$Geneid #Ensembl Gene ID
counts <- counts[, -1]  

# 02-Set the correct samples names (siOGT --> siRNA, siNC --> Negative CTRL) 

sample <- c("siOGT2-2", "siOGT2-1", "siOGT1-2", "siOGT1-1", "siNC-2", "siNC-1")

colnames(counts) <- sample
dim(counts)

counts <-counts[rowSums(counts) > 10, ]
dim(counts)

# 03-Build a colData with samples info and check the names and orders

colData <- data.frame(condition = factor(c("siOGT2", "siOGT2", "siOGT1", "siOGT1", "siNC", "siNC")))
rownames(colData) <- colnames(counts)

all(colnames(counts) %in% rownames(colData))
all(colnames(counts) == rownames(colData))

# 04-Build a DESeq Dataset object

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ condition)

# 05-Set the factor level and run DESeq function

dds$condition <- relevel(dds$condition, ref = "siNC")
dds <- DESeq(dds)

# 06-Variance stabilization using vst()

vsdata <- vst(dds, blind = FALSE)

pdf("01_PCA_vsdata.pdf")
pca <- plotPCA(vsdata, intgroup = "condition")
pca + ggtitle("PCA - KGN cells, siOGT vs Control")
dev.off()

pdf("02_DispEsts_dds.pdf")
plotDispEsts(dds,  main = "Dispersion Estimates - KGN cells, siOGT vs Control") 
dev.off()

# 07-Save output

message("DESeq2 object built with ", nrow(dds), " genes and ", ncol(dds), " samples")

saveRDS(dds, "dds.rds")
