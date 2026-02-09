library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)


# ------------------------------------------------
#            Gene Onthology (GO) on commong Genes
# ------------------------------------------------ 

# 01-Load list of common DEGs 

Gene_ID_bg <- readLines("background_genes.txt")

DEG_common_up <- readLines("DEG_common_up.txt")
DEG_common_down <- readLines("DEG_common_down.txt")

# 02-Search GO terms
# up-regulated genes 

GO_up <- enrichGO(gene = DEG_common_up,
                  OrgDb = "org.Hs.eg.db",
                  keyType = "ENSEMBL",
                  ont = "BP")

GO_up <- as.data.frame(GO_up)
dim(GO_up) 

message("No GO terms for up-regulated genes")

# down-regulated genes

GO_down<- enrichGO(gene = DEG_common_down, 
                   OrgDb = "org.Hs.eg.db", 
                   keyType = "ENSEMBL", 
                   ont = "BP")

GO_down<- as.data.frame(GO_down)
dim(GO_down)

message("Found ", nrow(GO_down), " terms for down-regulated genes")

GO_top_down <- GO_down |>
  dplyr::arrange(p.adjust) |>
  dplyr::slice_head(n = 20)


# barplot of top 20 GO terms

go <- ggplot(GO_top_down, aes(
  x = reorder(Description, Count),
  y = Count,
  fill = -log10(p.adjust)))+
  geom_col() +
  coord_flip() +
  scale_fill_viridis_c(name = "-log10(adj p-value)") +
  theme_classic() +
  labs(
    x = NULL,
    y = "Gene count",
    title = "GO enrichment – Biological Process"
  )
ggsave("13_GO_20BP.pdf", plot = go, width = 8, height = 6)

# ------------------------------------------------
#            Pathway Enrichment
# ------------------------------------------------ 

# 03-Pathway Enrichment using KEGG and Reactome
# Setup

DEG_common_down_entrezID <- mapIds(org.Hs.eg.db, 
                                   keys = DEG_common_down, 
                                   keytype = "ENSEMBL", 
                                   column = "ENTREZID", 
                                   multiVals = "first"
)

DEG_common_down_entrezID <- na.omit(DEG_common_down_entrezID)

background <- mapIds(org.Hs.eg.db,
                     keys = Gene_ID_bg,
                     keytype = "ENSEMBL",
                     column = "ENTREZID",
                     multiVals = "first")

# KEGG database

kegg <- enrichKEGG(gene = DEG_common_down_entrezID,
                   organism = "hsa",
                   universe = background, 
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1
)
kegg_df <- as.data.frame(kegg)

# Reactome databse

reactome <- enrichPathway(gene = DEG_common_down_entrezID,
                          organism = "human",
                          universe = background,
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.1)

reactome_df <- as.data.frame(reactome)

# Visualization of pathway enrichment
# KEGG

pdf("14_PathwayEnrichment_KEGG.pdf", width = 10, height = 10)

dotplot(kegg,
        showCategory = 15, 
        title = "KEGG Pathway Enrichment – Top 15 pathways (Down-regulated DEGs)"
) 

kegg_readable <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
cnetplot(kegg, 
         showCategory = 8, 
         node_label = "category"
) + 
  ggtitle("KEGG Pathway Enrichment: Top 8 Pathways") +
  theme(plot.title = element_text(hjust = 0.5))

kegg_df<- kegg_df|>
  dplyr::arrange(p.adjust) |>
  dplyr::slice_head(n = 20)

ggplot(kegg_df, aes(
  x = reorder(Description, Count),
  y = Count,
  fill = -log10(p.adjust)))+
  geom_col() +
  coord_flip() +
  scale_fill_viridis_c(name = "-log10(adj p-value)") +
  theme_classic() +
  labs(
    x = NULL,
    y = "Gene count",
    title = "KEGG  Pathway Enrichment – Downregulated DEGs"
  )

dev.off()

# Reactome

pdf("15_PathwayEnrichment_Reactome.pdf", width = 10, height = 10)

dotplot(reactome,
        showCategory = 20,
        title  = "Reactome Pathway Enrichment – Top 20 pathways (Down-regulated DEGs)"
) 

reactome_readable <- setReadable(reactome, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
cnetplot(reactome_readable, 
         showCategory = 8, 
         node_label = "category") + 
  ggtitle("Reactome Pathway Enrichment: Top 8 Pathways") +
  theme(plot.title = element_text(hjust = 0.5))

reactome_df <- reactome_df |>
  dplyr::arrange(p.adjust) |>
  dplyr::slice_head(n = 20)

ggplot(reactome_df, aes(
  x = reorder(Description, Count),
  y = Count,
  fill = -log10(p.adjust)))+
  geom_col() +
  coord_flip() +
  scale_fill_viridis_c(name = "-log10(adj p-value)") +
  theme_classic() +
  labs(
    x = NULL,
    y = "Gene count",
    title = "Reactome Pathway Enrichment – Downregulated DEGs"
  )

dev.off()

