library(Seurat)

load("/project2/gilad/emaan/time_project/cells.combined.RData")

pilot.humans <- readRDS("/project2/gilad/kenneth/share/pilot.humans.rds")
  
anchors <- FindTransferAnchors(
  reference = pilot.humans,
  query = cells.combined,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:50
)

pilot.humans=RunUMAP(pilot.humans,return.model=TRUE,reduction.name = "new.UMAP", dims=1:50)

seurat_revised_map <- MapQuery(
  anchorset = anchors,
  query = cells.combined,
  reference = pilot.humans,
  refdata = pilot.humans$labels,
  reference.reduction = "pca", 
  reduction.model = "new.UMAP"
)

save(seurat_revised_map, file = "seurat_revised_map_2.RData")

