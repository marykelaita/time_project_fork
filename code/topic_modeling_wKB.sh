#!/bin/bash

#SBATCH --job-name=initial_combined_wu
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

library(fastTopics)
library(Seurat)
library(Matrix)

load("/project/gilad/emaan/time_project/data/seurat_QC_2024.RData")
load("/project/gilad/emaan/time_project/data/GH_batch1_QC.RData")
load("/project/gilad/emaan/time_project/data/GH_batch2_QC.RData")
load("/project/gilad/emaan/time_project/code/pilot.humans.RData")

pilot.ecotderm <- subset(pilot.humans, subset = labels %in% c("Early Ectoderm 1", "Early Ectoderm 2", "Early Endoderm 1", "Early Endoderm 2"))
colnames(pilot.ecotderm@meta.data)[8] <- "condition"

load("/project/gilad/emaan/time_project/data/wu_combined_rev.RData")

merged_data <- merge(seurat_QC_2024, y = c(GH_batch1_QC, GH_batch2_QC,wu_combined_rev,pilot.ecotderm))

#Extract raw count matrix from seurat object and get it in correct format for fastTopics
#need to fit the model to the count data (unintegrated)
raw_counts <- merged_data@assays$RNA@counts
#remove genes without any counts in droplets
raw_counts <- raw_counts[rowSums(raw_counts > 0) > 0,] 
#get into correct orientation (barcodes x features)
raw_counts <- t(raw_counts)
dim(raw_counts)

fit <- fit_poisson_nmf(raw_counts,k = 4,numiter = 100)
saveRDS(fit, "/project/gilad/emaan/time_project/data/topicModel_k4_n100_wu_kb.rds")

fit5 <- fit_poisson_nmf(raw_counts,k = 5,numiter = 100)
saveRDS(fit5, "/project/gilad/emaan/time_project/data/topicModel_k5_n100_wu_kb.rds")

'