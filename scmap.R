library(SingleCellExperiment)
library(scmap)

# Mouse kidney 203 cells
# Note: make sure the working diretory is correct!
mouse_kidney_203_cds <- readRDS("mouse_kidney_203_cds.rds")

# All cells as reference set
mouse_kidney_203_data<- exprs(mouse_kidney_203_cds)
mouse_kidney_ann<- pData(mouse_kidney_203_cds)
cellname<- rownames(mouse_kidney_ann)
mouse_kidney_ann<- as.data.frame(mouse_kidney_ann[,4],stringsAsFactors = F)
colnames(mouse_kidney_ann)<- 'cell_type1'
rownames(mouse_kidney_ann)<- cellname
rm(mouse_kidney_203_cds,cellname)

# Create sce object
mouse_kidney_203_sce<- SingleCellExperiment(assays = list(normcounts = as.matrix(mouse_kidney_203_data)),colData = mouse_kidney_ann)
logcounts(mouse_kidney_203_sce) <- log2(normcounts(mouse_kidney_203_sce) + 1)
rowData(mouse_kidney_203_sce)$feature_symbol <- rownames(mouse_kidney_203_sce)
isSpike(mouse_kidney_203_sce, "ERCC") <- grepl("^ERCC-", rownames(mouse_kidney_203_sce))
mouse_kidney_203_sce <- mouse_kidney_203_sce[!duplicated(rownames(mouse_kidney_203_sce)), ]
mouse_kidney_203_sce <- selectFeatures(mouse_kidney_203_sce, suppress_plot = FALSE)

# scmap cell
set.seed(10)
mouse_kidney_203_sce <- indexCell(mouse_kidney_203_sce)
scmapCell_results <- scmapCell(
  mouse_kidney_203_sce, 
  list(
    scmap = metadata(mouse_kidney_203_sce)$scmap_cell_index
  )
)
scmapCell_clusters <- scmapCell2Cluster(
  scmapCell_results, 
  list(
    as.character(colData(mouse_kidney_203_sce)$cell_type1)
  )
)
mouse_kidney_203_scmap<- cbind(mouse_kidney_ann,scmapCell_clusters$scmap_cluster_labs)

# Human islet 1600 cells
human_islet_1600_cds <- readRDS("human_islet_1600_cds.rds")

# All cells as reference set
human_islet_1600_data<- exprs(human_islet_1600_cds)
human_islet_ann<- pData(human_islet_1600_cds)
cellname<- rownames(human_islet_ann)
human_islet_ann<- as.data.frame(human_islet_ann[,4],stringsAsFactors = F)
colnames(human_islet_ann)<- 'cell_type1'
rownames(human_islet_ann)<- cellname
rm(human_islet_1600_cds,cellname)

# Create sce object
human_islet_1600_sce<- SingleCellExperiment(assays = list(normcounts = as.matrix(human_islet_1600_data)),colData = human_islet_ann)
logcounts(human_islet_1600_sce) <- log2(normcounts(human_islet_1600_sce) + 1)
rowData(human_islet_1600_sce)$feature_symbol <- rownames(human_islet_1600_sce)
isSpike(human_islet_1600_sce, "ERCC") <- grepl("^ERCC-", rownames(human_islet_1600_sce))
human_islet_1600_sce <- human_islet_1600_sce[!duplicated(rownames(human_islet_1600_sce)), ]
human_islet_1600_sce <- selectFeatures(human_islet_1600_sce, suppress_plot = FALSE)

# scmap cell
set.seed(12)
human_islet_1600_sce <- indexCell(human_islet_1600_sce)
scmapCell_results <- scmapCell(
  human_islet_1600_sce, 
  list(
    scmap = metadata(human_islet_1600_sce)$scmap_cell_index
  )
)
scmapCell_clusters <- scmapCell2Cluster(
  scmapCell_results, 
  list(
    as.character(colData(human_islet_1600_sce)$cell_type1)
  )
)
human_islet_1600_scmap<- cbind(human_islet_ann,scmapCell_clusters$scmap_cluster_labs)


# Human pbmc 3694 cells
human_pbmc_3694_cds <- readRDS("human_pbmc_3694_cds.rds")

# All cells as reference set
human_pbmc_3694_data<- exprs(human_pbmc_3694_cds)
human_pbmc_ann<- pData(human_pbmc_3694_cds)
cellname<- rownames(human_pbmc_ann)
human_pbmc_ann<- as.data.frame(human_pbmc_ann[,4],stringsAsFactors = F)
colnames(human_pbmc_ann)<- 'cell_type1'
rownames(human_pbmc_ann)<- cellname
rm(human_pbmc_3694_cds,cellname)

# Create sce object
human_pbmc_3694_sce<- SingleCellExperiment(assays = list(normcounts = as.matrix(human_pbmc_3694_data)),colData = human_pbmc_ann)
logcounts(human_pbmc_3694_sce) <- log2(normcounts(human_pbmc_3694_sce) + 1)
rowData(human_pbmc_3694_sce)$feature_symbol <- rownames(human_pbmc_3694_sce)
isSpike(human_pbmc_3694_sce, "ERCC") <- grepl("^ERCC-", rownames(human_pbmc_3694_sce))
human_pbmc_3694_sce <- human_pbmc_3694_sce[!duplicated(rownames(human_pbmc_3694_sce)), ]
human_pbmc_3694_sce <- selectFeatures(human_pbmc_3694_sce, suppress_plot = FALSE)

# scmap cell
set.seed(14)
human_pbmc_3694_sce <- indexCell(human_pbmc_3694_sce)
scmapCell_results <- scmapCell(
  human_pbmc_3694_sce, 
  list(
    scmap = metadata(human_pbmc_3694_sce)$scmap_cell_index
  )
)
scmapCell_clusters <- scmapCell2Cluster(
  scmapCell_results, 
  list(
    as.character(colData(human_pbmc_3694_sce)$cell_type1)
  )
)
human_pbmc_3694_scmap<- cbind(human_pbmc_ann,scmapCell_clusters$scmap_cluster_labs)
