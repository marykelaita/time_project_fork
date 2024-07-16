library(Seurat)


filter <- ".filterC"
genes <- "Var"
dimred <- "cca"
kweigh <- 50
regout <- "reg"

load("/project/gilad/emaan/time_project/code/pilot.humans.RData")
wu_combined_rev <- readRDS("/project/gilad/ghousman/chondro-human-chimp/hc-chondro-time/chondro-time-evo/data/external/wu_combined_rev.rds")

obj <- c(pilot.humans,wu_combined_rev)


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

save(integrate, file = "/project/gilad/emaan/time_project/data/integrate_Wu_KB_var_cca_reg.RData")
