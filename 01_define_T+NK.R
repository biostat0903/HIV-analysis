#
rm(list=ls())
gc()

# Load packages
library(Seurat)
library(harmony)

library(bigreadr)
library(plyr)
library(dplyr)
library(stringr)
library(patchwork)
library(tibble)

# Set data parameters
# PROJ_PATH="D:\\hivProject\\"
PROJ_PATH="/public/home/biostat03/project/hivProject"
setwd(PROJ_PATH)
CELL_TYPE = "T_NK"

# Set scrna parameters
RESOLUTION <- 0.6
N_FEATURE <- 100
N_COUNT <- 20000
SEED <- 20240528

# Load data
load("03_result/03_cluster/all/seu_merge_obj.RData")
seu_merge_obj_sub <- subset(seu_merge_obj, 
                            subset = cluster == CELL_TYPE)


# Integrate samples
seu_split_obj <- SplitObject(seu_merge_obj_sub, split.by = "orig.ident") %>%
  llply(., function(seu_s){

    SCTransform(seu_s, variable.features.n = 2000, method = "glmGamPoi")
  })
# seu_split_objx = seu_split_obj
DefaultAssay(seu_merge_obj_sub) <- "SCT"
featureID <- SelectIntegrationFeatures(object.list = seu_split_obj,
                                       nfeatures = 3000) 
seu_split_obj <- PrepSCTIntegration(object.list = seu_split_obj,
                                    anchor.features = featureID)
seu_merge_obj_new <- merge(seu_split_obj[[1]], 
                           seu_split_obj[2:length(seu_split_obj)],
                           merge.data = TRUE)
DefaultAssay(seu_merge_obj)<-"RNA"
seu_merge_obj <- JoinLayers(seu_merge_obj)
use_assay <- "SCT"

seu_merge_obj_new <- RunPCA(seu_merge_obj_new, 
                            features = featureID, 
                            assay = use_assay, 
                            npcs = 50, 
                            verbose = FALSE) %>%
                     RunHarmony(., 
                                reduction.use = "pca",
                                reduction.save = "harmony",
                                assay.use = use_assay,
                                group.by.vars = "orig.ident")
eig_val <- (seu_merge_obj_new@reductions$harmony@stdev)^2
var_explained <- eig_val / sum(eig_val)
pc_clust <- min(which(cumsum(var_explained)>0.8))
cat(paste0("We select ", pc_clust, " PCs to explain 80% variance.\n"))

# Cluster 
DefaultAssay(seu_merge_obj_new) <- use_assay
seu_merge_obj_new <- FindNeighbors(seu_merge_obj_new, 
                               reduction = "harmony", 
                               dims = 1: pc_clust) %>% 
  FindClusters(.,
               n.iter = 200,
               algorithm = 1,
               resolution = RESOLUTION,
               random.seed = SEED) %>%
  RunUMAP(., 
          n.neighbors = 30,
          reduction = "harmony",
          dims = 1: pc_clust,
          n.epochs = 100,
          min.dist = 0.3,
          negative.sample.rate = 5L,
          seed.use = SEED)

# Output UMAP
tiff(file = "03_result/03_cluster/all/UMAP_raw_level2_T_NK.tiff",
     width = 8,height = 7,units = "in",res = 600,compression = "lzw")
DimPlot(seu_merge_obj_new, 
        reduction = "umap", 
        label = T) %>% print()
dev.off()

# Output dotplot
feature_list <- list(
  "T_features_list" = c("CD3G", "CD3D", "CD3E", "CD2"),
  "CD4_features_list" = c("CD4","FOXP3","CTLA4", "SOX4", "TNFSF13B"),
  "CD8_features_list" = c("CD8A","CD8B", "LTB", "S100B"),
  "NK_features_list" = c("KLRF1","KLRD1","NCAM1","CD160", "XCL1", "XCL2"),
  "MAIT_features_list" = c("SLC4A10","TRAV1-2"),
  "gdT_features_list" = c("TRGV9", "TRDV2", "TRGC1", "TRDC"),
  "ISG_features_list" = c("ISG15","IFI6","IFI44L","LY6E"),
  "Prolif_features_list" = c("MKI67","TYMS","PCNA"),
  "Naive_features_list" = c("CCR7","SELL","LEF1","TCF7"),
  "Effect_features_list" = c("GZMK","GZMA","GNLY","GZMB","PRF1","NKG7"),
  "Memery_features_list" = c("GPR183","S100A4"),
  "Exhaust_features_list" = c("HAVCR2", "TIGIT", "LAG3", "TOX"),
  "B_features_list" = c("CD79A", "MS4A1", "JSRP1", "MZB1", "CD38"),
  "PLT_features_list" = c("PF4", "PPBP", "GNG11")
)
tiff(file = "03_result/03_cluster/all/dotplot_level2_T_NK.tiff",
     width = 15,height = 8,units = "in",res = 200,compression = "lzw")
(DotPlot(seu_merge_obj_new, 
         features = unlist(feature_list) %>% unique) + 
    RotatedAxis()) %>% print
dev.off()

# Annotation
seurat_obj_newids <- c("CD8+ T",                                       ## 0
                       "NK", "CD4+ T", "CD4+ T", "CD8+ T", "CD8+ T",   ## 1-5
                       "CD8+ T", "CD8+ T", "CD4+ T", "MAIT", "gdT",    ## 6-10
                       "CD8+ T", "NKT", "CD8+ T", "CD8+ T", "NKT",     ## 11-15
                       "CD8+ T", "CD8+ T", "NK", "CD4+ T", "CD4+ T",      ## 16-20
                       "CD8+ T", "CD8+ T", "MAIT", "CD4+ T"            ## 21-24
)
Idents(seu_merge_obj_new) <-  seu_merge_obj_new@meta.data$seurat_cluster
names(seurat_obj_newids) <- levels(seu_merge_obj_new)
seu_merge_obj_new <- RenameIdents(seu_merge_obj_new, seurat_obj_newids)
tiff(file = paste0("03_result/03_cluster/all/UMAP_level2_T_NK-.tiff"),
     width = 8, height = 7, units = "in", res = 600, compression = "lzw")
DimPlot(seu_merge_obj_new, reduction = "umap", label = F)
dev.off()

# Save
seu_merge_obj_new@meta.data$cluster <- Idents(seu_merge_obj_new)
cluster_level2 <- seu_merge_obj_new@meta.data %>% 
  rownames_to_column(., var = "cell")
cluster_level2 <- cbind.data.frame(cluster_level2,
                                   Embeddings(seu_merge_obj_new, reduction = "umap"))
save(cluster_level2, file = "03_result/03_cluster/all/cluster_level2_T_NK.RData")
saveRDS(seu_merge_obj_new, 
        file = "03_result/03_cluster/all/seu_merge_obj_new_T_NK.rds")
