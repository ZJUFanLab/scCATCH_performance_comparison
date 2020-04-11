library(scCATCH)
library(SingleCellExperiment)
library(Seurat)

# Mouse kidney 203 cells
# Note: make sure the working diretory is correct!
mouse_kidney_203_cds <- readRDS("mouse_kidney_203_cds.rds")

# All cells as test set
mouse_kidney_203_data<- exprs(mouse_kidney_203_cds)
mouse_kidney_clu <- pData(mouse_kidney_203_cds)
cell_num <- unique(mouse_kidney_clu$FACS_type)
mouse_kidney_ann <- data.frame(cluster = as.character(1:length(cell_num)), cell_type = cell_num,stringsAsFactors = F)

for (i in 1:length(cell_num)) {
  mouse_kidney_clu[mouse_kidney_clu$FACS_type == cell_num[i],]$FACS_type<- mouse_kidney_ann$cluster[i]
}

cellname<- rownames(mouse_kidney_clu)
mouse_kidney_clu<- as.data.frame(mouse_kidney_clu[,4],stringsAsFactors = F)
colnames(mouse_kidney_clu)<- 'Cell_cluster'
rownames(mouse_kidney_clu)<- cellname
rm(mouse_kidney_203_cds,cellname,cell_num,i)

# Create Seurat object
mouse_kidney_203_Seurat <- CreateSeuratObject(counts = mouse_kidney_203_data,project = 'mouse_kidney_203')
mouse_kidney_203_Seurat <- NormalizeData(object = mouse_kidney_203_Seurat)
Idents(mouse_kidney_203_Seurat) <- factor(mouse_kidney_clu$Cell_cluster)

# scCATCH
mouse_kidney_203_marker_genes <- findmarkergenes(object = mouse_kidney_203_Seurat)
mouse_kidney_ann1<- scCATCH(object = mouse_kidney_203_marker_genes,species = 'Mouse',tissue = 'Kidney')
colnames(mouse_kidney_ann1)[3] <- 'scCATCH_cell_type'
mouse_kidney_ann<- cbind(mouse_kidney_ann,mouse_kidney_ann1)


# Human islet 1600 cells
human_islet_1600_cds <- readRDS("human_islet_1600_cds.rds")

# All cells as test set
human_islet_1600_data<- exprs(human_islet_1600_cds)
human_islet_clu <- pData(human_islet_1600_cds)
cell_num <- unique(human_islet_clu$FACS_type)
human_islet_ann <- data.frame(cluster = as.character(1:length(cell_num)), cell_type = cell_num,stringsAsFactors = F)

for (i in 1:length(cell_num)) {
  human_islet_clu[human_islet_clu$FACS_type == cell_num[i],]$FACS_type<- human_islet_ann$cluster[i]
}

cellname<- rownames(human_islet_clu)
human_islet_clu<- as.data.frame(human_islet_clu[,4],stringsAsFactors = F)
colnames(human_islet_clu)<- 'Cell_cluster'
rownames(human_islet_clu)<- cellname
rm(human_islet_1600_cds,cellname,cell_num,i)

# Create Seurat object
human_islet_1600_Seurat <- CreateSeuratObject(counts = human_islet_1600_data,project = 'human_islet_1600')
human_islet_1600_Seurat <- NormalizeData(object = human_islet_1600_Seurat)
Idents(human_islet_1600_Seurat) <- factor(human_islet_clu$Cell_cluster)

# scCATCH
human_islet_1600_marker_genes <- findmarkergenes(object = human_islet_1600_Seurat)
human_islet_ann1<- scCATCH(object = human_islet_1600_marker_genes,species = 'Human',tissue = c('Pancreas','Pancreatic islet'))
colnames(human_islet_ann1)[3] <- 'scCATCH_cell_type'
human_islet_ann<- cbind(human_islet_ann,human_islet_ann1)


# Human pbmc 3694 cells
human_pbmc_3694_cds <- readRDS("human_pbmc_3694_cds.rds")

# All cells as test set
human_pbmc_3694_data<- exprs(human_pbmc_3694_cds)
human_pbmc_clu <- pData(human_pbmc_3694_cds)
cell_num <- unique(human_pbmc_clu$FACS_type)
human_pbmc_ann <- data.frame(cluster = as.character(1:length(cell_num)), cell_type = cell_num,stringsAsFactors = F)

for (i in 1:length(cell_num)) {
  human_pbmc_clu[human_pbmc_clu$FACS_type == cell_num[i],]$FACS_type<- human_pbmc_ann$cluster[i]
}

cellname<- rownames(human_pbmc_clu)
human_pbmc_clu<- as.data.frame(human_pbmc_clu[,4],stringsAsFactors = F)
colnames(human_pbmc_clu)<- 'Cell_cluster'
rownames(human_pbmc_clu)<- cellname
rm(human_pbmc_3694_cds,cellname,cell_num,i)

# Create Seurat object
human_pbmc_3694_Seurat <- CreateSeuratObject(counts = human_pbmc_3694_data,project = 'human_pbmc_3694')
human_pbmc_3694_Seurat <- NormalizeData(object = human_pbmc_3694_Seurat)
Idents(human_pbmc_3694_Seurat) <- factor(human_pbmc_clu$Cell_cluster)

# scCATCH
human_pbmc_3694_marker_genes <- findmarkergenes(object = human_pbmc_3694_Seurat)
human_pbmc_ann1<- scCATCH(object = human_pbmc_3694_marker_genes,species = 'Human',tissue = c('Blood','Peripheral blood','Bone marrow'))
colnames(human_pbmc_ann1)[3] <- 'scCATCH_cell_type'
human_pbmc_ann<- cbind(human_pbmc_ann,human_pbmc_ann1)


# Human pancreas 2281 cells
human_pancreas_2281_cds <- readRDS("human_pancreas_2281_cds.rds")

# All cells as test set
human_pancreas_2281_data<- exprs(human_pancreas_2281_cds)
human_pancreas_clu <- pData(human_pancreas_2281_cds)
cell_num <- unique(human_pancreas_clu$FACS_type)
human_pancreas_ann <- data.frame(cluster = as.character(1:length(cell_num)), cell_type = cell_num,stringsAsFactors = F)

for (i in 1:length(cell_num)) {
  human_pancreas_clu[human_pancreas_clu$FACS_type == cell_num[i],]$FACS_type<- human_pancreas_ann$cluster[i]
}

cellname<- rownames(human_pancreas_clu)
human_pancreas_clu<- as.data.frame(human_pancreas_clu[,4],stringsAsFactors = F)
colnames(human_pancreas_clu)<- 'Cell_cluster'
rownames(human_pancreas_clu)<- cellname
rm(human_pancreas_2281_cds,cellname,cell_num,i)

# Create Seurat object
human_pancreas_2281_Seurat <- CreateSeuratObject(counts = human_pancreas_2281_data,project = 'human_pancreas_2281')
human_pancreas_2281_Seurat <- NormalizeData(object = human_pancreas_2281_Seurat)
Idents(human_pancreas_2281_Seurat) <- factor(human_pancreas_clu$Cell_cluster)

# scCATCH
human_pancreas_2281_marker_genes <- findmarkergenes(object = human_pancreas_2281_Seurat)
human_pancreas_ann1<- scCATCH(object = human_pancreas_2281_marker_genes,species = 'Human',tissue = c('Pancreas','Pancreatic islet'))
colnames(human_pancreas_ann1)[3] <- 'scCATCH_cell_type'
human_pancreas_ann<- cbind(human_pancreas_ann,human_pancreas_ann1)


# Mouse brain 20679 cells
mouse_brain_20679_cds <- readRDS("mouse_brain_20679_cds.rds")

# All cells as test set
mouse_brain_20679_data<- exprs(mouse_brain_20679_cds)
mouse_brain_clu <- pData(mouse_brain_20679_cds)
cell_num <- unique(mouse_brain_clu$FACS_type)
mouse_brain_ann <- data.frame(cluster = as.character(1:length(cell_num)), cell_type = cell_num,stringsAsFactors = F)

for (i in 1:length(cell_num)) {
  mouse_brain_clu[mouse_brain_clu$FACS_type == cell_num[i],]$FACS_type<- mouse_brain_ann$cluster[i]
}

cellname<- rownames(mouse_brain_clu)
mouse_brain_clu<- as.data.frame(mouse_brain_clu[,4],stringsAsFactors = F)
colnames(mouse_brain_clu)<- 'Cell_cluster'
rownames(mouse_brain_clu)<- cellname
rm(mouse_brain_20679_cds,cellname,cell_num,i)

# Create Seurat object
mouse_brain_20679_Seurat <- CreateSeuratObject(counts = mouse_brain_20679_data,project = 'mouse_brain_20679')
mouse_brain_20679_Seurat <- NormalizeData(object = mouse_brain_20679_Seurat)
Idents(mouse_brain_20679_Seurat) <- factor(mouse_brain_clu$Cell_cluster)

# scCATCH
mouse_brain_20679_marker_genes <- findmarkergenes(object = mouse_brain_20679_Seurat)
mouse_brain_ann1<- scCATCH(object = mouse_brain_20679_marker_genes,species = 'Mouse',tissue = 'Brain')
colnames(mouse_brain_ann1)[3] <- 'scCATCH_cell_type'
mouse_brain_ann<- cbind(mouse_brain_ann,mouse_brain_ann1)


# Human lung 2970 cells
human_lung_2970_cds <- readRDS("human_lung_2970_cds.rds")

# All cells as test set
human_lung_2970_data<- exprs(human_lung_2970_cds)
human_lung_clu <- pData(human_lung_2970_cds)
cell_num <- unique(human_lung_clu$FACS_type)
cell_num<- cell_num[c(2,5,3,7,6,1,4)]
human_lung_ann <- data.frame(cluster = as.character(1:length(cell_num)), cell_type = cell_num,stringsAsFactors = F)

for (i in 1:length(cell_num)) {
  human_lung_clu[human_lung_clu$FACS_type == cell_num[i],]$FACS_type<- human_lung_ann$cluster[i]
}

cellname<- rownames(human_lung_clu)
human_lung_clu<- as.data.frame(human_lung_clu[,4],stringsAsFactors = F)
colnames(human_lung_clu)<- 'Cell_cluster'
rownames(human_lung_clu)<- cellname
rm(human_lung_2970_cds,cellname,cell_num,i)

# Create Seurat object
human_lung_2970_Seurat <- CreateSeuratObject(counts = human_lung_2970_data,project = 'human_lung_2970')
human_lung_2970_Seurat <- NormalizeData(object = human_lung_2970_Seurat)
Idents(human_lung_2970_Seurat) <- factor(human_lung_clu$Cell_cluster)

# scCATCH
human_lung_2970_marker_genes <- findmarkergenes(object = human_lung_2970_Seurat)
human_lung_ann1<- scCATCH(object = human_lung_2970_marker_genes,species = 'Human',tissue = 'Lung')
colnames(human_lung_ann1)[3] <- 'scCATCH_cell_type'
human_lung_ann<- cbind(human_lung_ann,human_lung_ann1)


# Human pbmc 2638 cells
human_pbmc_2638_cds <- readRDS("human_pbmc_2638_cds.rds")

# All cells as test set
human_pbmc_2638_data<- exprs(human_pbmc_2638_cds)
human_pbmc_clu <- pData(human_pbmc_2638_cds)
cell_num <- unique(human_pbmc_clu$FACS_type)
cell_num<- cell_num[c(6,1,3,2,5,7,4,8,9)]
human_pbmc_ann <- data.frame(cluster = as.character(1:length(cell_num)), cell_type = cell_num,stringsAsFactors = F)

for (i in 1:length(cell_num)) {
  human_pbmc_clu[human_pbmc_clu$FACS_type == cell_num[i],]$FACS_type<- human_pbmc_ann$cluster[i]
}

cellname<- rownames(human_pbmc_clu)
human_pbmc_clu<- as.data.frame(human_pbmc_clu[,4],stringsAsFactors = F)
colnames(human_pbmc_clu)<- 'Cell_cluster'
rownames(human_pbmc_clu)<- cellname
rm(human_pbmc_2638_cds,cellname,cell_num,i)

# Create Seurat object
human_pbmc_2638_Seurat <- CreateSeuratObject(counts = human_pbmc_2638_data,project = 'human_pbmc_2638')
human_pbmc_2638_Seurat <- NormalizeData(object = human_pbmc_2638_Seurat)
Idents(human_pbmc_2638_Seurat) <- factor(human_pbmc_clu$Cell_cluster)

# scCATCH
human_pbmc_2638_marker_genes <- findmarkergenes(object = human_pbmc_2638_Seurat)
human_pbmc_ann1<- scCATCH(object = human_pbmc_2638_marker_genes,species = 'Human',tissue = c('Blood','Peripheral blood','Bone marrow'))
colnames(human_pbmc_ann1)[3] <- 'scCATCH_cell_type'
human_pbmc_ann<- cbind(human_pbmc_ann,human_pbmc_ann1)


# Mouse brain 2915 cells
mouse_brain_2915_cds <- readRDS("mouse_brain_2915_cds.rds")

# All cells as test set
mouse_brain_2915_data<- exprs(mouse_brain_2915_cds)
mouse_brain_clu <- pData(mouse_brain_2915_cds)
cell_num <- unique(mouse_brain_clu$FACS_type)
cell_num<- cell_num[c(1,2,3,4,5,7,8,6,9,10)]
mouse_brain_ann <- data.frame(cluster = as.character(1:length(cell_num)), cell_type = cell_num,stringsAsFactors = F)

for (i in 1:length(cell_num)) {
  mouse_brain_clu[mouse_brain_clu$FACS_type == cell_num[i],]$FACS_type<- mouse_brain_ann$cluster[i]
}

cellname<- rownames(mouse_brain_clu)
mouse_brain_clu<- as.data.frame(mouse_brain_clu[,4],stringsAsFactors = F)
colnames(mouse_brain_clu)<- 'Cell_cluster'
rownames(mouse_brain_clu)<- cellname
rm(mouse_brain_2915_cds,cellname,cell_num,i)

# Create Seurat object
mouse_brain_2915_Seurat <- CreateSeuratObject(counts = mouse_brain_2915_data,project = 'mouse_brain_2915')
mouse_brain_2915_Seurat <- NormalizeData(object = mouse_brain_2915_Seurat)
Idents(mouse_brain_2915_Seurat) <- factor(mouse_brain_clu$Cell_cluster)

# scCATCH
mouse_brain_2915_marker_genes <- findmarkergenes(object = mouse_brain_2915_Seurat)
mouse_brain_ann1<- scCATCH(object = mouse_brain_2915_marker_genes,species = 'Mouse',tissue = 'Brain')
colnames(mouse_brain_ann1)[3] <- 'scCATCH_cell_type'
mouse_brain_ann<- cbind(mouse_brain_ann,mouse_brain_ann1)


# Mouse brain 3918 cells
mouse_brain_3918_cds <- readRDS("mouse_brain_3918_cds.rds")

# All cells as test set
mouse_brain_3918_data<- exprs(mouse_brain_3918_cds)
mouse_brain_clu <- pData(mouse_brain_3918_cds)
cell_num <- unique(mouse_brain_clu$FACS_type)
mouse_brain_ann <- data.frame(cluster = as.character(1:length(cell_num)), cell_type = cell_num,stringsAsFactors = F)

for (i in 1:length(cell_num)) {
  mouse_brain_clu[mouse_brain_clu$FACS_type == cell_num[i],]$FACS_type<- mouse_brain_ann$cluster[i]
}

cellname<- rownames(mouse_brain_clu)
mouse_brain_clu<- as.data.frame(mouse_brain_clu[,4],stringsAsFactors = F)
colnames(mouse_brain_clu)<- 'Cell_cluster'
rownames(mouse_brain_clu)<- cellname
rm(mouse_brain_3918_cds,cellname,cell_num,i)

# Create Seurat object
mouse_brain_3918_Seurat <- CreateSeuratObject(counts = mouse_brain_3918_data,project = 'mouse_brain_3918')
mouse_brain_3918_Seurat <- NormalizeData(object = mouse_brain_3918_Seurat)
Idents(mouse_brain_3918_Seurat) <- factor(mouse_brain_clu$Cell_cluster)

# scCATCH
mouse_brain_3918_marker_genes <- findmarkergenes(object = mouse_brain_3918_Seurat)
mouse_brain_ann1<- scCATCH(object = mouse_brain_3918_marker_genes,species = 'Mouse',tissue = 'Brain')
colnames(mouse_brain_ann1)[3] <- 'scCATCH_cell_type'
mouse_brain_ann<- cbind(mouse_brain_ann,mouse_brain_ann1)
