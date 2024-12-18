############
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

# Load data
load("03_result/03_cluster/all/cluster_level1.RData")
load("03_result/03_cluster/all/cluster_level2_MP.RData")
cluster_level2_MP <- cluster_level2
load("03_result/03_cluster/all/cluster_level2_T_NK.RData")
cluster_level2_T_NK <- cluster_level2

# Redefine 
cluster_level1$cluster_refine <- cluster_level1$cluster %>% as.character()
cluster_level1$cluster_refine[match(cluster_level2_T_NK$cell, cluster_level1$cell)] <- 
  cluster_level2_T_NK$cluster %>% as.character()
cluster_level1$cluster_refine[match(cluster_level2_MP$cell, cluster_level1$cell)] <- 
  cluster_level2_MP$cluster %>% as.character()
unique(cluster_level1$cluster_refine)


cluster_level1$cluster_refine <- 
  factor(cluster_level1$cluster_refine,
         levels = c("CD4+ T", "CD8+ T", "gdT", "MAIT", "NKT", "NK",
                    "B", "Plasma",
                    "CD14+ Monocyte", "CD16+ Monocyte", "cDC", "pDC",
                    "PLT", "RBC"))
saveRDS(cluster_level1, 
        file = "03_result/03_cluster/all/cluster_new.rds")

# UMAP for cell types
CT_COLS <- c("#98C1D9", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51", 
             "#E0FBFC", "#A8DADC", "#3D5A80", "#81B29A", "#E07A5F", 
             "#DBC9D8", "#b388eb", "#A4277C", "#BC93B2", "#0077b6", 
             "#BB3E03", "#FFDDD2", "#F19C79", "#006D77", "#6A3569",
             "#E6194B", "#4363D8", "#FFE119", "#3CB44B", "#F58231", 
             "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE", 
             "#008080", "#E6BEFF", "#9A6324", "#FFFAC8", "#800000", 
             "#AAFFC3", "#808000", "#FFD8B1", "#000075", "#808080")
cluster_level1_sub <- subset(cluster_level1, !cluster_level1$cluster_refine %in% c("PLT", "RBC"))
umap_plt <- ggplot(cluster_level1_sub,
                   aes(x = umap_1, y = umap_2, color = cluster_refine)) + 
  geom_point(alpha = 1, size = 0.6) + 
  scale_color_manual("Cell Types", values = CT_COLS) + 
  guides(color = guide_legend(byrow = TRUE, 
                              ncol = 6, 
                              keyheight = 0.4)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())
tiff(file = paste0("03_result/03_cluster/all/UMAP_ct_new.tiff"),
     width = 9, height = 7,units = "in",res = 300,compression = "lzw")
plot(umap_plt)
dev.off()

## UMAP for sample
sample_info <- fread2("01_data/sampleInfo.txt")
cluster_level1_s <- merge(cluster_level1, sample_info, 
                          by.x = "orig.ident", by.y = "code")
cluster_level1_s_sub <- subset(cluster_level1_s, !cluster_refine %in% c("PLT", "RBC"))

for (gg in c("HIV-supress", "HIV-progress", "Control")){
  
  cluster_level1_s_tmp <- subset(cluster_level1_s_sub, group == gg)
  cat(dim(cluster_level1_s_tmp), "\n")
  umap_plt <- ggplot(cluster_level1_s_tmp,
                     aes(x = umap_1, y = umap_2, color = cluster_refine)) + 
    geom_point(alpha = 1, size = 0.6) + 
    scale_color_manual("Cell Types", values = CT_COLS) + 
    guides(color = guide_legend(byrow = TRUE, 
                                ncol = 6, 
                                keyheight = 0.4)) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank())
  tiff(file = paste0("03_result/03_cluster/all/UMAP_", gg, "_new.tiff"),
       width = 9, height = 7,units = "in",res = 300,compression = "lzw")
  plot(umap_plt)
  dev.off()
  
}

# Proportion
prop_df = cluster_level1_s_sub %>%
  group_by(group, cluster_refine) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)*100)
prop_plt <- ggplot(prop_df, aes(x = group, y = freq, fill = cluster_refine)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  scale_fill_manual("Cell type", values = CT_COLS) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold")) +
  xlab("Groups") +
  ylab("Cell type proportions") +
  coord_flip()
tiff(file = "03_result/03_cluster/all/CT_proportion.tiff",
     width = 9, height = 7,units = "in",res = 300,compression = "lzw")
plot(prop_plt)
dev.off()
