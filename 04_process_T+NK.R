


rm(list=ls())
gc()
library(bigreadr)
library(Seurat)
library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)

# Set parameters
PROJ_PATH="/public/home/biostat03/project/hivProject/"
setwd(PROJ_PATH)
CT_COLS <- c("#BB3E03", "#FFDDD2", "#F19C79", "#006D77", "#6A3569",
             "#E6194B", "#4363D8", "#FFE119", "#3CB44B", "#F58231", 
             "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE", 
             "#008080", "#E6BEFF", "#9A6324", "#FFFAC8", "#800000", 
             "#AAFFC3", "#808000", "#FFD8B1", "#000075", "#808080")

# Load object
seu_merge_obj_T_NK <- readRDS("03_result/03_cluster/all/seu_merge_obj_new_T_NK.rds")
sample_info <- fread2("01_data/sampleInfo.txt")
cluster_T_NK_s <- merge(seu_merge_obj_T_NK@meta.data, sample_info, 
                        by.x = "orig.ident", by.y = "code")

# Add labels
cluster_T_NK_s$cluster_level3 <- as.character(cluster_T_NK_s$seurat_cluster)
cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 0] <- "CD8+ T Effect_1"

cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 1] <- "NK_1"
cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 2] <- "CD4+ T Naive_1"
cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 3] <- "CD4+ T Memory"
cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 4] <- "CD8+ T Exhaust_1"
cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 5] <- "CD8+ T Naive_1"

cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 6] <- "CD8+ T Effect_2"
cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 7] <- "ISG"
cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 8] <- "CD4+ T Naive_3"
cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 9] <- "MAIT_1"
cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 10] <- "gdT"

cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 11] <- "CD8+ T Effect_3"
cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 12] <- "NKT_1"
cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 13] <- "CD8+ T Effect_4"
cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 14] <- "CD8+ T Naive_2"
cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 15] <- "NKT_2"

cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 16] <- "CD8+ T Effect_5"
cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 17] <- "CD8+ T Effect_6"
cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 18] <- "NK_2"
cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 19] <- "CD4+ T Effect"
cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 20] <- "NKT_3"

cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 21] <- "CD8+ T Naive_3"
cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 22] <- "CD8+ T Exhaust_2"
cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 23] <- "MAIT_2"
cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 24] <- "CD4+ T Naive_2"

# Define DEGs
Idents(seu_merge_obj_T_NK) <- cluster_T_NK_s$cluster_level3 
de_markers <- FindMarkers(object = seu_merge_obj_T_NK,
                          test.use = "wilcox",
                          only.pos = T,
                          min.pct = 0.1, 
                          return.thresh = 0.1,
                          logfc.threshold = 0)
head(x = markers)





prop_df = cluster_T_NK_s %>%
  group_by(group, seurat_clusters) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)*100)

cluster_T_NK_s$cluster_level3[cluster_T_NK_s$cluster_level3 == 0] <- "CD8+ T Effect"prop_plt <- ggplot(prop_df, aes(x = group, y = freq, fill = seurat_clusters)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  scale_fill_manual("Cell type", values = CT_COLS) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold")) +
  xlab("Groups") +
  ylab("Cell type proportions") +
  coord_flip()
tiff(file = "03_result/03_cluster/all/CT_T_NK_no_annot_proportion.tiff",
     width = 9, height = 7,units = "in",res = 300,compression = "lzw")
plot(prop_plt)
dev.off()


prop_df = cluster_T_NK_s %>%
  group_by(group, cluster) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)*100)

prop_plt <- ggplot(prop_df, aes(x = group, y = freq, fill = cluster)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  scale_fill_manual("Cell type", values = CT_COLS) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold")) +
  xlab("Groups") +
  ylab("Cell type proportions") +
  coord_flip()
tiff(file = "03_result/03_cluster/all/CT_T_NK_annot_proportion.tiff",
     width = 9, height = 7,units = "in",res = 300,compression = "lzw")
plot(prop_plt)
dev.off()
