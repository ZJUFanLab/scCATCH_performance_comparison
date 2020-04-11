library(CHETAH)

# Mouse kidney 203 cells
# Note: make sure the working diretory is correct!
mouse_kidney_203_cds <- readRDS("mouse_kidney_203_cds.rds")
mouse_kidney_ref <- readRDS("mouse_kidney_markers.rds")

# All cells as test set
mouse_kidney_203_data <- exprs(mouse_kidney_203_cds)
mouse_kidney_ann <- pData(mouse_kidney_203_cds)
mouse_kidney_umap <- mouse_kidney_ann[,c(1,2)]
cellname <- rownames(mouse_kidney_ann)
mouse_kidney_ann <- as.data.frame(mouse_kidney_ann[,4],stringsAsFactors = F)
colnames(mouse_kidney_ann) <- 'Cell_type'
rownames(mouse_kidney_ann) <- cellname
rm(mouse_kidney_203_cds,cellname)

# Generating reference
cellname <- rep('ref',ncol(mouse_kidney_ref))
cellname <- paste(cellname,1:length(cellname),sep = '_')
mouse_kidney_ref1 <- colnames(mouse_kidney_ref)
names(mouse_kidney_ref1) <- cellname
colnames(mouse_kidney_ref) <- cellname
rm(cellname)

mouse_kidney_ref<- mouse_kidney_ref[rownames(mouse_kidney_ref) %in% rownames(mouse_kidney_203_data),]
cellname<- rep(1,ncol(mouse_kidney_ref))
for (i in 1:ncol(mouse_kidney_ref)) {
  d1<- mouse_kidney_ref[,i]
  if (all(d1 == 0)) {
    cellname[i]<- 0
  }
}
cellname1<- which(cellname == 1)
mouse_kidney_ref<- mouse_kidney_ref[,cellname1]
mouse_kidney_ref1<- mouse_kidney_ref1[cellname1]

# Create sce object
mouse_kidney_ref_sce <- SingleCellExperiment(assays = list(counts = as.matrix(mouse_kidney_ref)),colData = DataFrame(celltypes = mouse_kidney_ref1))
mouse_kidney_203_sce <- SingleCellExperiment(assays = list(counts = as.matrix(mouse_kidney_203_data)),reducedDims = SimpleList(UMAP = as.matrix(mouse_kidney_umap)))

# CHETAH
mouse_kidney_203_sce <- CHETAHclassifier(input = mouse_kidney_203_sce,ref_cells = mouse_kidney_ref_sce)
mouse_kidney_ann$CHETAH_cellmatch_ref <- mouse_kidney_203_sce$celltype_CHETAH

# Error! 
# mouse_kidney_203_sce <- CHETAHclassifier(input = mouse_kidney_203_sce,ref_cells = headneck_ref)


# Human islet 1600 cells
# Note: make sure the working diretory is correct!
human_islet_1600_cds <- readRDS("human_islet_1600_cds.rds")
human_islet_ref <- readRDS("human_islet_markers.rds")

# 1. All cells as test set
human_islet_1600_data <- exprs(human_islet_1600_cds)
human_islet_ann <- pData(human_islet_1600_cds)
human_islet_umap <- human_islet_ann[,c(1,2)]
cellname <- rownames(human_islet_ann)
human_islet_ann <- as.data.frame(human_islet_ann[,4],stringsAsFactors = F)
colnames(human_islet_ann) <- 'Cell_type'
rownames(human_islet_ann) <- cellname
rm(human_islet_1600_cds,cellname)

# Generating reference
cellname <- rep('ref',ncol(human_islet_ref))
cellname <- paste(cellname,1:length(cellname),sep = '_')
human_islet_ref1 <- colnames(human_islet_ref)
names(human_islet_ref1) <- cellname
colnames(human_islet_ref) <- cellname
rm(cellname)

human_islet_ref<- human_islet_ref[rownames(human_islet_ref) %in% rownames(human_islet_1600_data),]
cellname<- rep(1,ncol(human_islet_ref))
for (i in 1:ncol(human_islet_ref)) {
  d1<- human_islet_ref[,i]
  if (all(d1 == 0)) {
    cellname[i]<- 0
  }
}
cellname1<- which(cellname == 1)
human_islet_ref<- human_islet_ref[,cellname1]
human_islet_ref1<- human_islet_ref1[cellname1]

# Create sce object
human_islet_ref_sce <- SingleCellExperiment(assays = list(counts = as.matrix(human_islet_ref)),colData = DataFrame(celltypes = human_islet_ref1))
human_islet_1600_sce <- SingleCellExperiment(assays = list(counts = as.matrix(human_islet_1600_data)),reducedDims = SimpleList(UMAP = as.matrix(human_islet_umap)))

# CHETAH
human_islet_1600_sce <- CHETAHclassifier(input = human_islet_1600_sce,ref_cells = human_islet_ref_sce)
human_islet_ann$CHETAH_cellmatch_ref <- human_islet_1600_sce$celltype_CHETAH

human_islet_1600_sce1 <- CHETAHclassifier(input = human_islet_1600_sce,ref_cells = headneck_ref)
human_islet_ann$CHETAH_ref <- human_islet_1600_sce1$celltype_CHETAH


# Human pbmc 3694 cells
# Note: make sure the working diretory is correct!
human_pbmc_3694_cds <- readRDS("human_pbmc_3694_cds.rds")
human_pbmc_ref <- readRDS("human_pbmc_markers.rds")

# All cells as test set
human_pbmc_3694_data <- exprs(human_pbmc_3694_cds)
human_pbmc_ann <- pData(human_pbmc_3694_cds)
human_pbmc_umap <- human_pbmc_ann[,c(1,2)]
cellname <- rownames(human_pbmc_ann)
human_pbmc_ann <- as.data.frame(human_pbmc_ann[,4],stringsAsFactors = F)
colnames(human_pbmc_ann) <- 'Cell_type'
rownames(human_pbmc_ann) <- cellname
rm(human_pbmc_3694_cds,cellname)

# Generating reference
cellname <- rep('ref',ncol(human_pbmc_ref))
cellname <- paste(cellname,1:length(cellname),sep = '_')
human_pbmc_ref1 <- colnames(human_pbmc_ref)
names(human_pbmc_ref1) <- cellname
colnames(human_pbmc_ref) <- cellname
rm(cellname)

human_pbmc_ref<- human_pbmc_ref[rownames(human_pbmc_ref) %in% rownames(human_pbmc_3694_data),]
cellname<- rep(1,ncol(human_pbmc_ref))
for (i in 1:ncol(human_pbmc_ref)) {
  d1<- human_pbmc_ref[,i]
  if (all(d1 == 0)) {
    cellname[i]<- 0
  }
}
cellname1<- which(cellname == 1)
human_pbmc_ref<- human_pbmc_ref[,cellname1]
human_pbmc_ref1<- human_pbmc_ref1[cellname1]

# Create sce object
human_pbmc_ref_sce <- SingleCellExperiment(assays = list(counts = as.matrix(human_pbmc_ref)),colData = DataFrame(celltypes = human_pbmc_ref1))
human_pbmc_3694_sce <- SingleCellExperiment(assays = list(counts = as.matrix(human_pbmc_3694_data)),reducedDims = SimpleList(UMAP = as.matrix(human_pbmc_umap)))

# CHETAH 
human_pbmc_3694_sce <- CHETAHclassifier(input = human_pbmc_3694_sce,ref_cells = human_pbmc_ref_sce)
human_pbmc_ann$CHETAH_cellmatch_ref <- human_pbmc_3694_sce$celltype_CHETAH

human_pbmc_3694_sce1 <- CHETAHclassifier(input = human_pbmc_3694_sce,ref_cells = headneck_ref)
human_pbmc_ann$CHETAH_ref <- human_pbmc_3694_sce1$celltype_CHETAH



