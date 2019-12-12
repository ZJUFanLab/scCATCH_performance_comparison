library(SingleR)
library(Biobase)
# Note: SingleR package was intalled from github by (devtools::install_github('dviraran/SingleR'))

# Mouse kidney 203 cells
# Note: make sure the working diretory is correct!
mouse_kidney_203_cds <- readRDS("mouse_kidney_203_cds.rds")
mouse_kidney_markers <- readRDS("mouse_kidney_markers.rds")

# All cells as test set
mouse_kidney_203_data <- exprs(mouse_kidney_203_cds)
mouse_kidney_ann <- pData(mouse_kidney_203_cds)
cellname <- rownames(mouse_kidney_ann)
mouse_kidney_ann <- as.data.frame(mouse_kidney_ann[,4],stringsAsFactors = F)
colnames(mouse_kidney_ann) <- 'cell_type'
rownames(mouse_kidney_ann) <- cellname
rm(mouse_kidney_203_cds,cellname)

# Creating a new reference data set
celltypes <- colnames(mouse_kidney_markers)
cellmaintypes <- colnames(mouse_kidney_markers)
mouse_kidney_ref <- list(name = 'CellMatch',data = as.matrix(mouse_kidney_markers),
                        types = celltypes,main_types = cellmaintypes)
mouse_kidney_ref$de.genes <- CreateVariableGeneSet(as.matrix(mouse_kidney_markers),celltypes,200)
mouse_kidney_ref$de.genes.main <- CreateVariableGeneSet(as.matrix(mouse_kidney_markers),cellmaintypes,300)
rm(celltypes,cellmaintypes)

# Create singleR object
mouse_kidney_203_singleR <- CreateSinglerObject(counts = as.matrix(mouse_kidney_203_data),
                                                project.name = 'mouse_kidney_203',
                                                species = 'Mouse',
                                                ref.list = list(immgen, mouse_kidney_ref, mouse.rnaseq))

mouse_kidney_203_singleR <- mouse_kidney_203_singleR$singler

mouse_kidney_203_singleR1 <- mouse_kidney_203_singleR[[1]]
mouse_kidney_203_singleR2 <- mouse_kidney_203_singleR[[2]]
mouse_kidney_203_singleR3 <- mouse_kidney_203_singleR[[3]]

mouse_kidney_203_singleR1 <- mouse_kidney_203_singleR1$SingleR.single.main
mouse_kidney_203_singleR2 <- mouse_kidney_203_singleR2$SingleR.single.main
mouse_kidney_203_singleR3 <- mouse_kidney_203_singleR3$SingleR.single.main
mouse_kidney_ann$immgen <- mouse_kidney_203_singleR1$labels1
mouse_kidney_ann$cellmatch <- mouse_kidney_203_singleR2$labels1
mouse_kidney_ann$mouse.rnaseq <- mouse_kidney_203_singleR3$labels1
rm(mouse_kidney_203_singleR1,mouse_kidney_203_singleR2,mouse_kidney_203_singleR3)


# Human islet 1600 cells
human_islet_1600_cds <- readRDS("human_islet_1600_cds.rds")
human_islet_markers <- readRDS("human_islet_markers.rds")

# All cells as test set
human_islet_1600_data <- exprs(human_islet_1600_cds)
human_islet_ann <- pData(human_islet_1600_cds)
cellname <- rownames(human_islet_ann)
human_islet_ann <- as.data.frame(human_islet_ann[,4],stringsAsFactors = F)
colnames(human_islet_ann) <- 'cell_type'
rownames(human_islet_ann) <- cellname
rm(human_islet_1600_cds,cellname)

# Creating a new reference data set
celltypes <- colnames(human_islet_markers)
cellmaintypes <- colnames(human_islet_markers)
human_islet_ref<- list(name = 'CellMatch',data = as.matrix(human_islet_markers),
                        types = celltypes,main_types = cellmaintypes)
human_islet_ref$de.genes <- CreateVariableGeneSet(as.matrix(human_islet_markers),celltypes,200)
human_islet_ref$de.genes.main <- CreateVariableGeneSet(as.matrix(human_islet_markers),cellmaintypes,300)
rm(celltypes,cellmaintypes)

# Create singleR object
human_islet_1600_singleR <- CreateSinglerObject(counts = as.matrix(human_islet_1600_data),
                                                project.name = 'human_islet_1600',
                                                species = 'Human',
                                                ref.list = list(hpca, human_islet_ref, blueprint_encode))

human_islet_1600_singleR <- human_islet_1600_singleR$singler

human_islet_1600_singleR1 <- human_islet_1600_singleR[[1]]
human_islet_1600_singleR2 <- human_islet_1600_singleR[[2]]
human_islet_1600_singleR3 <- human_islet_1600_singleR[[3]]

human_islet_1600_singleR1 <- human_islet_1600_singleR1$SingleR.single.main
human_islet_1600_singleR2 <- human_islet_1600_singleR2$SingleR.single.main
human_islet_1600_singleR3 <- human_islet_1600_singleR3$SingleR.single.main
human_islet_ann$hpca <- human_islet_1600_singleR1$labels1
human_islet_ann$cellmatch <- human_islet_1600_singleR2$labels1
human_islet_ann$blueprint_encode <- human_islet_1600_singleR3$labels1
rm(human_islet_1600_singleR1,human_islet_1600_singleR2,human_islet_1600_singleR3)


# Human pbmc 3694 cells
human_pbmc_3694_cds <- readRDS("human_pbmc_3694_cds.rds")
human_pbmc_markers <- readRDS("human_pbmc_markers.rds")

# All cells as test set
human_pbmc_3694_data <- exprs(human_pbmc_3694_cds)
human_pbmc_ann <- pData(human_pbmc_3694_cds)
cellname <- rownames(human_pbmc_ann)
human_pbmc_ann <- as.data.frame(human_pbmc_ann[,4],stringsAsFactors = F)
colnames(human_pbmc_ann) <- 'cell_type'
rownames(human_pbmc_ann) <- cellname
rm(human_pbmc_3694_cds,cellname)

# Creating a new reference data set
celltypes <- colnames(human_pbmc_markers)
cellmaintypes <- colnames(human_pbmc_markers)
d1<- grep(x = cellmaintypes,pattern = 'T cell')
cellmaintypes[d1]<- 'T cell'
d1<- grep(x = cellmaintypes,pattern = 'T helper')
cellmaintypes[d1]<- 'T cell'
d1<- grep(x = cellmaintypes,pattern = 'B cell')
cellmaintypes[d1]<- 'B cell'
d1<- grep(x = cellmaintypes,pattern = 'monocyte')
cellmaintypes[d1]<- 'Monocyte'
d1<- grep(x = cellmaintypes,pattern = 'stromal cell')
cellmaintypes[d1]<- 'Stromal cell'
d1<- grep(x = cellmaintypes,pattern = 'stroma cell')
cellmaintypes[d1]<- 'Stromal cell'
d1<- grep(x = cellmaintypes,pattern = 'dendritic cell')
cellmaintypes[d1]<- 'Dendritic cell'
d1<- grep(x = cellmaintypes,pattern = 'macrophage')
cellmaintypes[d1]<- 'Macrophage'
d1<- grep(x = cellmaintypes,pattern = 'endothelial cell')
cellmaintypes[d1]<- 'Endothelial cell'

human_pbmc_ref<- list(name = 'CellMatch',data = as.matrix(human_pbmc_markers),
                       types = celltypes,main_types = cellmaintypes)
human_pbmc_ref$de.genes <- CreateVariableGeneSet(as.matrix(human_pbmc_markers),celltypes,200)
human_pbmc_ref$de.genes.main <- CreateVariableGeneSet(as.matrix(human_pbmc_markers),cellmaintypes,300)
rm(celltypes,cellmaintypes)

# Create singleR object
human_pbmc_3694_singleR <- CreateSinglerObject(counts = as.matrix(human_pbmc_3694_data),
                                                project.name = 'human_pbmc_3694',
                                                species = 'Human',
                                                ref.list = list(hpca, human_pbmc_ref, blueprint_encode))

human_pbmc_3694_singleR <- human_pbmc_3694_singleR$singler

human_pbmc_3694_singleR1 <- human_pbmc_3694_singleR[[1]]
human_pbmc_3694_singleR2 <- human_pbmc_3694_singleR[[2]]
human_pbmc_3694_singleR3 <- human_pbmc_3694_singleR[[3]]

human_pbmc_3694_singleR1 <- human_pbmc_3694_singleR1$SingleR.single.main
human_pbmc_3694_singleR2 <- human_pbmc_3694_singleR2$SingleR.single.main
human_pbmc_3694_singleR3 <- human_pbmc_3694_singleR3$SingleR.single.main
human_pbmc_ann$hpca <- human_pbmc_3694_singleR1$labels1
human_pbmc_ann$cellmatch <- human_pbmc_3694_singleR2$labels1
human_pbmc_ann$blueprint_encode <- human_pbmc_3694_singleR3$labels1
rm(human_pbmc_3694_singleR1,human_pbmc_3694_singleR2,human_pbmc_3694_singleR3)

