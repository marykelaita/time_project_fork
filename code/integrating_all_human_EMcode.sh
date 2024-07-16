#!/bin/bash

#SBATCH --job-name=initial_combined
#SBATCH --time=36:00:00
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem-per-cpu=32G
#SBATCH --account=pi-gilad

module unload R
module unload gcc

module load R/4.2.0
module load gcc/10.2.0

Rscript -e '

library(Seurat)

load("/project/gilad/emaan/time_project/data/GH_batch2_QC.RData")
load("/project/gilad/emaan/time_project/data/GH_batch1_QC.RData")
load("/project/gilad/emaan/time_project/data/seurat_QC_2024.RData")

#Seurat objects being used
NA18856 <- subset(seurat_QC_2024, subset = donor_id == "NA18856")
NA18855 <- subset(seurat_QC_2024, subset = donor_id == "NA18855")
NA19160 <- subset(seurat_QC_2024, subset = donor_id == "NA19160")
H1 <- subset(GH_batch1_QC, subset = donor_ID == "H1")
H2 <- subset(GH_batch1_QC, subset = donor_ID == "H2")
H3 <- subset(GH_batch2_QC, subset = donor_ID == "H3")
H5 <- subset(GH_batch2_QC, subset = donor_ID == "H5")

#Creating a list of the Seurat objects
individual_list <- list(NA18856, NA18855, NA19160, H1,H2,H3,H5)

# normalize and identify variable features for each dataset independently
individual_list <- lapply(X = individual_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

#Select features that are repeatedly variable across data sets for integration
features <- SelectIntegrationFeatures(object.list = individual_list)
cell.anchors <- FindIntegrationAnchors(object.list = individual_list, anchor.features = features)

#Creating the integrated data set
cells.combined <- IntegrateData(anchorset = cell.anchors)
DefaultAssay(cells.combined) <- "integrated"

save(cells.combined, file = "/project/gilad/emaan/time_project/data/integrating_EM_GH/cells.combined_integrated_DA_reg.RData")

'