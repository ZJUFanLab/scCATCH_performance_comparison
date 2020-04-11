library(SingleCellExperiment)
library(cellassign)
library(scran)

# Mouse kidney 203 cells
# Note: make sure the working diretory is correct!
mouse_kidney_203_cds <- readRDS("mouse_kidney_203_cds.rds")
mouse_kidney_markers <- readRDS("mouse_kidney_markers.rds")

mouse_kidney_203_data <- exprs(mouse_kidney_203_cds)
mouse_kidney_ann <- pData(mouse_kidney_203_cds)
cellname <- rownames(mouse_kidney_ann)
mouse_kidney_ann<- as.data.frame(mouse_kidney_ann[,4],stringsAsFactors = F)
colnames(mouse_kidney_ann) <- 'Cell'
rownames(mouse_kidney_ann) <- cellname
rm(mouse_kidney_203_cds,cellname)

mouse_kidney_markers<- mouse_kidney_markers[rownames(mouse_kidney_markers) %in% rownames(mouse_kidney_203_data),]
cellname<- rep(1,ncol(mouse_kidney_markers))
for (i in 1:ncol(mouse_kidney_markers)) {
  d1<- mouse_kidney_markers[,i]
  if (all(d1 == 0)) {
    cellname[i]<- 0
  }
}
cellname1<- which(cellname == 1)
mouse_kidney_markers<- mouse_kidney_markers[,cellname1]

mouse_kidney_ann1 <- mouse_kidney_ann
mouse_kidney_ann1$Cell <- rownames(mouse_kidney_ann1)

mouse_kidney_ann2 <- data.frame(Gene = rownames(mouse_kidney_203_data),stringsAsFactors = F)
rownames(mouse_kidney_ann2) <- mouse_kidney_ann2$Gene

# Create sce object
mouse_kidney_203_sce <- SingleCellExperiment(assays = list(counts = as.matrix(mouse_kidney_203_data)),colData = mouse_kidney_ann1,rowData = mouse_kidney_ann2)
mouse_kidney_203_sce <- computeSumFactors(mouse_kidney_203_sce)
mouse_kidney_203_sce_Size_factors <- sizeFactors(mouse_kidney_203_sce)

# cellassign
mouse_kidney_203_fit <- cellassign(exprs_obj = mouse_kidney_203_sce[rownames(mouse_kidney_markers),], 
                                   marker_gene_info = as.matrix(mouse_kidney_markers), 
                                   s = mouse_kidney_203_sce_Size_factors, 
                                   learning_rate = 1e-2, 
                                   shrinkage = TRUE,
                                   verbose = TRUE)

mouse_kidney_ann$cellassign<- celltypes(mouse_kidney_203_fit)


# Human islet 1600 cells
human_islet_1600_cds <- readRDS("human_islet_1600_cds.rds")
human_islet_markers <- readRDS("human_islet_markers.rds")

human_islet_1600_data <- exprs(human_islet_1600_cds)
human_islet_ann <- pData(human_islet_1600_cds)
cellname <- rownames(human_islet_ann)
human_islet_ann<- as.data.frame(human_islet_ann[,4],stringsAsFactors = F)
colnames(human_islet_ann) <- 'Cell'
rownames(human_islet_ann) <- cellname
rm(human_islet_1600_cds,cellname)

human_islet_markers<- human_islet_markers[rownames(human_islet_markers) %in% rownames(human_islet_1600_data),]
cellname<- rep(1,ncol(human_islet_markers))
for (i in 1:ncol(human_islet_markers)) {
  d1<- human_islet_markers[,i]
  if (all(d1 == 0)) {
    cellname[i]<- 0
  }
}
cellname1<- which(cellname == 1)
human_islet_markers<- human_islet_markers[,cellname1]

human_islet_ann1 <- human_islet_ann
human_islet_ann1$Cell <- rownames(human_islet_ann1)

human_islet_ann2 <- data.frame(Gene = rownames(human_islet_1600_data),stringsAsFactors = F)
rownames(human_islet_ann2) <- human_islet_ann2$Gene

# Create sce object
human_islet_1600_sce <- SingleCellExperiment(assays = list(counts = as.matrix(human_islet_1600_data)),colData = human_islet_ann1,rowData = human_islet_ann2)
human_islet_1600_sce <- computeSumFactors(human_islet_1600_sce)
human_islet_1600_sce_Size_factors <- sizeFactors(human_islet_1600_sce)

# cellassign
human_islet_1600_fit <- cellassign(exprs_obj = human_islet_1600_sce[rownames(human_islet_markers),], 
                                   marker_gene_info = as.matrix(human_islet_markers), 
                                   s = human_islet_1600_sce_Size_factors, 
                                   learning_rate = 1e-2, 
                                   shrinkage = TRUE,
                                   verbose = TRUE)

human_islet_ann$cellassign<- celltypes(human_islet_1600_fit)


# Human pbmc 3694 cells
human_pbmc_3694_cds <- readRDS("human_pbmc_3694_cds.rds")
human_pbmc_markers <- readRDS("human_pbmc_markers.rds")

human_pbmc_3694_data <- exprs(human_pbmc_3694_cds)
human_pbmc_ann <- pData(human_pbmc_3694_cds)
cellname <- rownames(human_pbmc_ann)
human_pbmc_ann<- as.data.frame(human_pbmc_ann[,4],stringsAsFactors = F)
colnames(human_pbmc_ann) <- 'Cell'
rownames(human_pbmc_ann) <- cellname
rm(human_pbmc_3694_cds,cellname)

human_pbmc_markers<- human_pbmc_markers[rownames(human_pbmc_markers) %in% rownames(human_pbmc_3694_data),]
cellname<- rep(1,ncol(human_pbmc_markers))
for (i in 1:ncol(human_pbmc_markers)) {
  d1<- human_pbmc_markers[,i]
  if (all(d1 == 0)) {
    cellname[i]<- 0
  }
}
cellname1<- which(cellname == 1)
human_pbmc_markers<- human_pbmc_markers[,cellname1]

human_pbmc_ann1 <- human_pbmc_ann
human_pbmc_ann1$Cell <- rownames(human_pbmc_ann1)

human_pbmc_ann2 <- data.frame(Gene = rownames(human_pbmc_3694_data),stringsAsFactors = F)
rownames(human_pbmc_ann2) <- human_pbmc_ann2$Gene

# Create sce object
human_pbmc_3694_sce <- SingleCellExperiment(assays = list(counts = as.matrix(human_pbmc_3694_data)),colData = human_pbmc_ann1,rowData = human_pbmc_ann2)
human_pbmc_3694_sce <- computeSumFactors(human_pbmc_3694_sce)
human_pbmc_3694_sce_Size_factors <- sizeFactors(human_pbmc_3694_sce)

# cellassign
human_pbmc_3694_fit <- cellassign(exprs_obj = human_pbmc_3694_sce[rownames(human_pbmc_markers),], 
                                   marker_gene_info = as.matrix(human_pbmc_markers), 
                                   s = human_pbmc_3694_sce_Size_factors, 
                                   learning_rate = 1e-2, 
                                   shrinkage = TRUE,
                                   verbose = TRUE)

human_pbmc_ann$cellassign<- celltypes(human_pbmc_3694_fit)



