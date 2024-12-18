
rm(list = ls())

# Load packages
library(Seurat)
library(harmony)
library(scDblFinder)

library(bigreadr)
library(plyr)
library(dplyr)
library(stringr)
library(tibble)

library(patchwork)

# Set data parameters
# PROJ_PATH="D:\\hivProject\\"
PROJ_PATH="/public/home/biostat03/project/hivProject"
setwd(PROJ_PATH)
SAMPLE_INFO <- fread2("./01_data/sampleInfo.txt")

# Set scrna parameters
MIN_CELLS <- 10
MIN_FEATURE <- 50
RESOLUTION <- 0.6
N_FEATURE <- 100
N_COUNT <- 20000
PERCENT_MT <- 10
SCT_BOOL <- T
DOUBLETS_BOOL <- T
SEED <- 20240528

# Quality control scRNA-seq data 
seu_merge_list <- alply(SAMPLE_INFO[, 1], 1, function(ss){
  
  count_data <- Read10X(paste0("./01_data/rna/", ss, ".matrix/"))
  message("Before QC, the data contains ", dim(count_data)[2], 
          " cells and ", dim(count_data)[1], " genes.\n")
  seu_obj <- CreateSeuratObject(counts = count_data, 
                                min.cells = MIN_CELLS, 
                                min.features = MIN_FEATURE)
  
  ## Remove doublets
  if (DOUBLETS_BOOL){
    
    sce_obj <- as.SingleCellExperiment(seu_obj)
    sce_obj <- scDblFinder(sce_obj, dbr = 0.1, verbose = F)
    singlet_cell <- sce_obj[, sce_obj$scDblFinder.class == "singlet"] 
    seu_obj <- subset(seu_obj, cells = colnames(singlet_cell))
  }
  
  ## Remove cells with too low expression or too high mt genes expression
  seu_obj$orig.ident <- ss
  seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj, pattern = "^MT-")
  seu_obj <- subset(seu_obj, nFeature_RNA >= N_FEATURE &
                      nCount_RNA <= N_COUNT & percent.mt <= PERCENT_MT)
  
  ## Remove MT, miRNA, lncRNA, HB, and ENSG genes
  gene_delete_id <- c(rownames(seu_obj)[grepl("^MT-", rownames(seu_obj))], 
                      rownames(seu_obj)[grepl("^MIR", rownames(seu_obj))],
                      rownames(seu_obj)[grepl("^LINC", rownames(seu_obj))], 
                      # rownames(seu_obj)[grepl("^HB", rownames(seu_obj))], 
                      rownames(seu_obj)[grepl("^ENSG", rownames(seu_obj))])
  seu_obj <- subset(seu_obj, features = setdiff(rownames(seu_obj), gene_delete_id))
  message("After QC, the data contains ", dim(seu_obj@assays$RNA$counts)[2], 
          " cells and ", dim(seu_obj@assays$RNA$counts)[1], " genes.\n")
  
  ## Normalize scRNA
  if (SCT_BOOL == T){
    
    seu_obj <- SCTransform(seu_obj, variable.features.n = 2000, method = "glmGamPoi")
  } else {
    
    seu_obj <- NormalizeData(seu_obj, verbose = FALSE)
    seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)
  }
  seu_obj <- RenameCells(seu_obj, add.cell.id = ss)
  return(seu_obj)
}) 

# Integrate samples
featureID <- SelectIntegrationFeatures(object.list = seu_merge_list, 
                                       nfeatures = 2000)
if (SCT_BOOL == T){
  
  seu_merge_list <- PrepSCTIntegration(object.list = seu_merge_list,
                                       anchor.features = featureID)
  seu_merge_obj <- merge(x= seu_merge_list[[1]],
                         y= seu_merge_list[2:length(seu_merge_list)],
                         merge.data = TRUE)
  DefaultAssay(seu_merge_obj)<-"RNA"
  seu_merge_obj <- JoinLayers(seu_merge_obj)
  use_assay <- "SCT"
} else {
  
  seu_merge_obj <- Reduce("merge", seu_merge_list)
  VariableFeatures(seu_merge_obj) <- featureID
  seu_merge_obj <- ScaleData(seu_merge_obj, 
                             features = featureID,
                             verbose = FALSE)
  DefaultAssay(seu_merge_obj)<-"RNA"
  seu_merge_obj <- JoinLayers(seu_merge_obj)
  use_assay <- "RNA"
}
seu_merge_obj <- RunPCA(seu_merge_obj, 
                        features = featureID, 
                        assay = use_assay, 
                        npcs = 50, 
                        verbose = FALSE)
seu_merge_obj <- RunHarmony(seu_merge_obj, 
                            reduction.use = "pca",
                            reduction.save = "harmony",
                            assay.use = use_assay,
                            group.by.vars = "orig.ident")
eig_val <- (seu_merge_obj@reductions$harmony@stdev)^2
var_explained <- eig_val / sum(eig_val)
pc_clust <- min(which(cumsum(var_explained)>0.9))
cat(paste0("We select ", pc_clust, " PCs to explain 90% variance.\n"))

# Cluster
DefaultAssay(seu_merge_obj) <- use_assay
seu_merge_obj <- FindNeighbors(seu_merge_obj, 
                               reduction = "harmony", 
                               dims = 1: pc_clust) 
seu_merge_obj <- FindClusters(seu_merge_obj,
                              n.iter = 200,
                              algorithm = 1,
                              resolution = RESOLUTION,
                              random.seed = SEED)
seu_merge_obj <- RunUMAP(seu_merge_obj, 
                         n.neighbors = 30,
                         reduction = "harmony",
                         dims = 1: pc_clust,
                         n.epochs = 100,
                         min.dist = 0.3,
                         negative.sample.rate = 5L,
                         seed.use = SEED)

# Output UMAP
tiff(file = "03_result/03_cluster/all/UMAP_raw_level1-.tiff",
     width = 8,height = 7,units = "in",res = 600,compression = "lzw")
DimPlot(seu_merge_obj, 
        reduction = "umap", 
        label = T) %>% print()
dev.off()

# Output dotplot
feature_list <- list(
  "T_features_list" = c("CD2", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B"),
  "NK_features_list" = c("KLRF1","KLRD1","NCAM1","CD160"),
  "MP_features_list" = c("CD14", "CD163", "CD68", "CSF1R", "FCGR3A", "S100A8"),
  "DC_features_list" = c("CD1C", "LILRA4", "CD1E", "FCER1A"),
  "B_features_list" = c("CD79A", "MS4A1", "JSRP1", "MZB1", "CD38"),
  "Gra_features_list" = c("ENPP3", "ITGAM", "FUT4"),
  "PLT_features_list" = c("PF4", "PPBP", "GNG11"),
  "Prolif_features_list" = c("MKI67","TYMS","PCNA"),
  "HB_features_list" = c("HBD", "HBM", "AHSP")
)
tiff(file = "03_result/03_cluster/all/dotplot_level1-.tiff",
     width = 15,height = 8,units = "in",res = 200,compression = "lzw")
(DotPlot(seu_merge_obj, 
         features = unlist(feature_list) %>% unique) + 
    RotatedAxis()) %>% print
dev.off()

# Annotation
seurat_obj_newids <- c("T_NK",                                     ## 0
                       "T_NK", "MP", "T_NK", "T_NK", "T_NK",       ## 1-5
                       "T_NK", "T_NK", "B", "T_NK", "B",           ## 6-10 
                       "T_NK", "MP", "T_NK", "MP", "PLT",          ## 11-15 
                       "T_NK", "T_NK", "DC", "RBC", "Plasma",      ## 16-20
                       "T_NK", "T_NK", "B", "MP", "DC",            ## 21-25
                       "T_NK", "MP", "T_NK", "T_NK", "T_NK",       ## 26-30
                       "MP"                                        ## 31
)
names(seurat_obj_newids) <- levels(seu_merge_obj)
seu_merge_obj <- RenameIdents(seu_merge_obj, seurat_obj_newids)
tiff(file = paste0("03_result/03_cluster/all/UMAP_level1.tiff"),
     width = 8, height = 7, units = "in", res = 600, compression = "lzw")
DimPlot(seu_merge_obj, reduction = "umap", label = F)
dev.off()

# Save
seu_merge_obj@meta.data$cluster <- Idents(seu_merge_obj)
cluster_level1 <- seu_merge_obj@meta.data %>% 
  rownames_to_column(., var = "cell")
cluster_level1 <- cbind.data.frame(cluster_level1,
                                   Embeddings(seu_merge_obj, reduction = "umap"))
save(cluster_level1, file = "03_result/03_cluster/all/cluster_level1.RData")
save(seu_merge_obj, file = "03_result/03_cluster/all/seu_merge_obj.RData")


