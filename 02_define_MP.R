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
CELL_TYPE = c("MP", "DC")

# Set scrna parameters
RESOLUTION <- 0.6
N_FEATURE <- 100
N_COUNT <- 20000
SEED <- 20240528

# Load data
load("03_result/03_cluster/all/seu_merge_obj.RData")
seu_merge_obj_sub <- subset(seu_merge_obj, 
                            subset = cluster %in% CELL_TYPE)


# Integrate samples
seu_split_obj <- SplitObject(seu_merge_obj_sub, split.by = "orig.ident") %>%
  llply(., function(seu_s){
    
    SCTransform(seu_s, variable.features.n = 2000, method = "glmGamPoi")
  })
DefaultAssay(seu_merge_obj_sub) <- "SCT"
featureID <- SelectIntegrationFeatures(object.list = seu_split_obj,
                                       nfeatures = 3000) 
seu_split_obj <- PrepSCTIntegration(object.list = seu_split_obj,
                                     anchor.features = featureID)
seu_merge_obj_new <- merge(seu_split_obj[[1]], 
                           seu_split_obj[2:length(seu_split_obj)],
                           merge.data = TRUE)
DefaultAssay(seu_merge_obj_new)<-"RNA"
seu_merge_obj_new <- JoinLayers(seu_merge_obj_new)
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
tiff(file = "03_result/03_cluster/all/UMAP_raw_level2_MP.tiff",
     width = 8,height = 7,units = "in",res = 600,compression = "lzw")
DimPlot(seu_merge_obj_new, 
        reduction = "umap", 
        label = T) %>% print()
dev.off()


# Output dotplot
feature_list <- list(
  "Mono14_features_list" = c("S100A8", "S100A12", "VCAN", "MNDA", "FCN1", "CD14"),
  "C1Mono16_features_list" = c("C1QA", "C1QB", "C1QC"),
  "Mono16_features_list" = c("FCGR3A", "SMIM25", "TCF7L2", "RHOC", "CDKN1C"),
  "cDC_features_list" = c("FCER1A", "CD1C", "CLEC10A", "CLEC9A", "CST3"),
  "pDC_features_list" = c("LILRA4", "SERPINF1", "IL3RA","CLEC4C"),
  "NK_features_list" = c("KLRF1","KLRD1","NCAM1","CD160", "XCL1", "XCL2"),
  "T_features_list" = c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B"),
  "PLT_features_list" = c("PF4", "PPBP", "GNG11")
)
tiff(file = "03_result/03_cluster/all/dotplot_level2_MP.tiff",
     width = 15,height = 8,units = "in",res = 200,compression = "lzw")
(DotPlot(seu_merge_obj_new, 
         features = unlist(feature_list) %>% unique) + 
    RotatedAxis()) %>% print
dev.off()

# Annotation
seurat_obj_newids <- c("CD14+ Monocyte",                                      ## 0
                       "CD14+ Monocyte", "CD14+ Monocyte", "CD14+ Monocyte", "CD16+ Monocyte", "CD14+ Monocyte",  ## 1-5
                       "CD14+ Monocyte", "CD14+ Monocyte", "cDC", "pDC", "CD14+ Monocyte",  ## 6-10
                       "CD16+ Monocyte", "pDC"               ## 11-12
                       
)
Idents(seu_merge_obj_new) <-  seu_merge_obj_new@meta.data$seurat_clusters
names(seurat_obj_newids) <- levels(seu_merge_obj_new)
seu_merge_obj_new <- RenameIdents(seu_merge_obj_new, seurat_obj_newids)
tiff(file = paste0("03_result/03_cluster/all/UMAP_level2_MP.tiff"),
     width = 8, height = 7, units = "in", res = 600, compression = "lzw")
DimPlot(seu_merge_obj_new, reduction = "umap", label = F)
dev.off()

# Save
seu_merge_obj_new@meta.data$cluster <- Idents(seu_merge_obj_new)
cluster_level2 <- seu_merge_obj_new@meta.data %>% 
  rownames_to_column(., var = "cell")
cluster_level2 <- cbind.data.frame(cluster_level2,
                                   Embeddings(seu_merge_obj_new, reduction = "umap"))
save(cluster_level2, file = "03_result/03_cluster/all/cluster_level2_MP.RData")
saveRDS(seu_merge_obj_new, 
        file = "03_result/03_cluster/all/seu_merge_obj_new_MP.rds")

