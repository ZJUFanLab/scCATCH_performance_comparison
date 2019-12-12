library(garnett)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

# Mouse kidney 203 cells
# Note: make sure the working diretory is correct!
mouse_kidney_203_cds <- readRDS("mouse_kidney_203_cds.rds")
mmKidney_markers <- 'mmKidney_markers.txt'

# All cells as training set
# Train and get classifier
mouse_kidney_203_classifier <- train_cell_classifier(cds = mouse_kidney_203_cds,
                                                     marker_file = mmKidney_markers,
                                                     db = org.Mm.eg.db,
                                                     cds_gene_id_type = "SYMBOL",
                                                     num_unknown = 50,
                                                     marker_file_gene_id_type = "SYMBOL")
# Classify mouse kidney 203 cells
mouse_kidney_203_cds <- classify_cells(mouse_kidney_203_cds, mouse_kidney_203_classifier,
                                       db = org.Mm.eg.db,
                                       cluster_extend = TRUE,
                                       cds_gene_id_type = "SYMBOL")



# Human islet 1600 cells
human_islet_1600_cds <- readRDS("human_islet_1600_cds.rds")
hsIslet_markers <- 'hsIslet_markers.txt'

# All cells as training set
# Train and get classifier
human_islet_1600_classifier <- train_cell_classifier(cds = human_islet_1600_cds,
                                                     marker_file = hsIslet_markers,
                                                     db = org.Hs.eg.db,
                                                     cds_gene_id_type = "SYMBOL",
                                                     num_unknown = 50,
                                                     marker_file_gene_id_type = "SYMBOL")
# Classify human islet 1600 cells
human_islet_1600_cds <- classify_cells(human_islet_1600_cds, human_islet_1600_classifier,
                                       db = org.Hs.eg.db,
                                       cluster_extend = TRUE,
                                       cds_gene_id_type = "SYMBOL")



# Human PBMC 3694 cells
human_pbmc_3694_cds <- readRDS('human_pbmc_3694_cds.rds')
hsPBMC_markers <- 'hsPBMC_markers.txt'

# All cells as training set
# Train and get classifier
human_pbmc_3694_classifier <- train_cell_classifier(cds = human_pbmc_3694_cds,
                                                     marker_file = hsPBMC_markers,
                                                     db = org.Hs.eg.db,
                                                     cds_gene_id_type = "SYMBOL",
                                                     num_unknown = 50,
                                                     marker_file_gene_id_type = "SYMBOL")
# Classify human pbmc 3694 cells
human_pbmc_3694_cds <- classify_cells(human_pbmc_3694_cds, human_pbmc_3694_classifier,
                                       db = org.Hs.eg.db,
                                       cluster_extend = TRUE,
                                       cds_gene_id_type = "SYMBOL")


