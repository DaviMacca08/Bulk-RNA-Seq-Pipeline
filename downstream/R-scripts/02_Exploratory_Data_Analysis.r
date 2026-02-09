library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(DESeq2)


# ------------------------------------------
#            Exploratory Data Analysis
# ------------------------------------------

# 01-Load DESeq object

dds <- readRDS("dds.rds")
colData <- as.data.frame(colData(dds))["condition"]

message("Loaded DESeq2 object with ", nrow(dds), " genes and ", ncol(dds), " samples")

vsdata <- vst(dds, blind = FALSE)


# 02-Box-Plot vsdata

vsdata_df <- as.data.frame(assay(vsdata))
vsdata_df <- pivot_longer(cbind(gene = rownames(vsdata_df), vsdata_df),
                          -gene, names_to = "sample", 
                          values_to = "expr")

pdf("03_BoxPlot_VstData.pdf")
ggplot(vsdata_df, aes(x = sample, y = expr, fill = sample)) + 
  geom_boxplot() + 
  theme_classic() +
  scale_fill_manual(
    values = c(
    "siNC-1" = "#A6A6A6",
    "siNC-2" = "#A6A6A6",
    "siOGT1-1" = "#8DD3C7",
    "siOGT1-2" = "#8DD3C7",
    "siOGT2-1" = "#FB8072",
    "siOGT2-2" = "#FB8072"
  )) +
  labs(
    title = "Boxplot of VST counts - KGN cells",
    y = "VST value",
    x = "Sample"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# 03-Take 500 most variable genes from VstData

vars <- rowVars(assay(vsdata))
vars <- sort(vars, decreasing = TRUE)
vars <- vars[1:500]
topVars <- assay(vsdata)[names(vars),]

# 04-Calculate euclidean distances and cluster 500 most variable genes

dist <- dist(t(topVars), method = "euclidean")
dist_all <- dist(t(assay(vsdata)), method = "euclidean")

pdf("04_HierarchicalClustering_AllMethods.pdf")
plot(hclust(dist, method = "ward.D2"), main = "Hierarchical Clustering (Ward.D2) - Top 500 Most Variable Genes")
plot(hclust(dist, method = "complete"), main = "Hierarchical Clustering (Complete) - Top 500 Most Variable Genes")
plot(hclust(dist, method = "average"), main = "Hierarchical Clustering (Average) - Top 500 Most Variable Genes")
plot(hclust(dist_all, method = "ward.D2"), main = "Hierarchical Clustering (Ward.D2) - All Genes")
dev.off()

# 05-Heatmap with gene and samples dendograms of 500 most variables genes

colData$condition <- factor(colData$condition)
htmap_color <- colorRampPalette(c("dodgerblue4", "white", "darkorange2"))(100)
my_colors <- list(
  condition = c(
    siNC   = "#A6A6A6",  
    siOGT1 = "#8DD3C7",  
    siOGT2 = "#FB8072"  
  )
)

pheatmap(
  topVars,
  scale = "row",
  color = htmap_color,
  annotation_col = colData,
  annotation_colors = my_colors, 
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_method = "ward.D2",
  show_rownames = FALSE,
  fontsize_col = 12,
  main = "Heatmap of top 500 most variable genes (VST)",
  filename = "05_Heatmap_Top500.pdf",
  width = 9,
  height = 9
)

# 06-Heatmap sample-sample correlation

correlation_sample <- cor(assay(vsdata)) 

pheatmap(correlation_sample, 
         main="Sample-Sample Correlation",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         filename = "06_Heatmap_SampleCor.pdf",
         width = 7,
         height = 7
           )


# 07-Principal Components (PCA) of 500 most variables genes

pca_topVars <- prcomp(t(topVars)) 
names(pca_topVars)
summary(pca_topVars)

var_exp <- (pca_topVars$sdev^2) / sum(pca_topVars$sdev^2)

pdf("07_PrincipalComponent_Top500.pdf")
barplot(var_exp * 100,
        names.arg = paste0("PC", seq_along(var_exp)),
        las = 2,
        ylab = "Variance explained (%)",
        main = "PCA: Variance Explained by Principal Components",
        col = "grey70",
        border = NA
        )
dev.off()

PC1_PC2 <- data.frame(PC1 = pca_topVars$x[,1], 
                      PC2 = pca_topVars$x[,2],
                      condition = colData$condition)

pdf("08_PCA_Top500.pdf")
ggplot(PC1_PC2, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1 (", round(var_exp[1]*100,1), "%)")) +
  ylab(paste0("PC2 (", round(var_exp[2]*100,1), "%)")) +
  theme_bw() +
  ggtitle("PCA of top 500 most variable genes")
dev.off()

