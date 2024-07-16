#!/bin/bash

#SBATCH --job-name=UMAP_KB_combined
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=10G
#SBATCH --account=pi-gilad

module unload R
module unload gcc

module load R/4.2.0
module load gcc/10.2.0

Rscript -e '

library(Seurat)


filter <- ".filterC"
genes <- "Var"
dimred <- "cca"
kweigh <- 50
regout <- "reg"

load("/project/gilad/emaan/time_project/data/seurat_QC_2024.RData")
load("/project/gilad/emaan/time_project/data/GH_batch1_QC.RData")
load("/project/gilad/emaan/time_project/data/GH_batch2_QC.RData")
load("/project/gilad/emaan/time_project/code/pilot.humans.RData")

NA18856 <- subset(seurat_QC_2024, subset = donor_id == "NA18856")
NA18855 <- subset(seurat_QC_2024, subset = donor_id == "NA18855")
NA19160 <- subset(seurat_QC_2024, subset = donor_id == "NA19160")
H1 <- subset(GH_batch1_QC, subset = donor_ID == "H1")
H2 <- subset(GH_batch1_QC, subset = donor_ID == "H2")
H3 <- subset(GH_batch2_QC, subset = donor_ID == "H3")
H5 <- subset(GH_batch2_QC, subset = donor_ID == "H5")

obj <- c(NA18856, NA18855, NA19160, H1, H2, H3, H5, pilot.humans)


if (genes=="Var") {
  obj.features <- SelectIntegrationFeatures(object.list = obj)
  length(obj.features)
}

obj <- lapply(X=obj, FUN=function(x) {
  x <- ScaleData(x, features=obj.features, verbose=FALSE)
  x <- RunPCA(x, features=obj.features, verbose=FALSE)
})

if (dimred=="cca") {
  obj.anchors <- FindIntegrationAnchors(object.list=obj,
                                        normalization.method="LogNormalize",
                                        anchor.features=obj.features,
                                        reduction="cca")
}

integrate <- IntegrateData(anchorset=obj.anchors,
                           normalization.method="LogNormalize",
                           k.weight=kweigh)


if (regout=="reg") {
  integrate <- ScaleData(integrate, vars.to.regress=c("nCount_RNA","percent.mt"))
}

integrate <- RunPCA(object=integrate,
                    npcs=100,
                    verbose=FALSE)


pva <- integrate@reductions$pca@stdev^2/integrate@reductions$pca@misc$total.variance
ndim <- length(which(pva>=0.001))
print(ndim)


integrate <- RunUMAP(integrate,
                     dims=1:ndim)

save(integrate, file = "/project/gilad/emaan/time_project/data/integrate_indv_KB_var_cca_reg.RData")

'