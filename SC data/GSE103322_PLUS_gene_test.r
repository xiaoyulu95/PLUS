
#############################
#loading data
###########################
counts0<-as.matrix(read.table("GSE103322_malignant_cells.txt"))
MAT_seurat<- CreateSeuratObject(counts = counts0, project = "MAT", min.cells = 3, min.features = 200)
names(attributes(MAT_seurat))
dim(attributes(MAT_seurat)$assays$RNA)

MAT_seurat[["percent.mt"]] <- PercentageFeatureSet(MAT_seurat, pattern = "^MT-")

MAT_seurat <- NormalizeData(MAT_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
MAT_seurat <- NormalizeData(MAT_seurat)
MAT_seurat <- FindVariableFeatures(MAT_seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MAT_seurat)
MAT_seurat <- ScaleData(MAT_seurat, features = all.genes)


############################
library(Seurat)
load("GSE103322_PLUS_gene_test.RData")
##############################
MAT_seurat <- RunPCA(MAT_seurat, features = meta.var.genes)

MAT_seurat <- FindNeighbors(MAT_seurat, dims = 1:15)
MAT_seurat <- FindClusters(MAT_seurat, resolution = 0.5)


MAT_seurat <- RunUMAP(MAT_seurat, dims = setdiff(c(1:10),5))
#DimPlot(MAT_seurat, reduction = "pca")
#DimHeatmap(MAT_seurat, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(MAT_seurat, dims = 1:15, cells = 500, balanced = TRUE)

DimPlot(MAT_seurat, reduction = "umap")
Idents(MAT_seurat)<-ident3
DimPlot(MAT_seurat, reduction = "umap",pt.size =1)

Idents(MAT_seurat)<-attributes(MAT_seurat)$meta.data[,1]
DimPlot(MAT_seurat, reduction = "umap",pt.size =1)


###############
MAT_seurat <- RunPCA(MAT_seurat, features = all.var.genes)
MAT_seurat <- FindNeighbors(MAT_seurat, dims = 1:15)
MAT_seurat <- FindClusters(MAT_seurat, resolution = 0.5)

MAT_seurat <- RunUMAP(MAT_seurat, dims =setdiff(1:25,c(6)))
#DimPlot(MAT_seurat, reduction = "pca")
#DimHeatmap(MAT_seurat, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(MAT_seurat, dims = 1:15, cells = 500, balanced = TRUE)

DimPlot(MAT_seurat, reduction = "umap")
Idents(MAT_seurat)<-ident3
DimPlot(MAT_seurat, reduction = "umap",pt.size =1)

Idents(MAT_seurat)<-attributes(MAT_seurat)$meta.data[,1]
DimPlot(MAT_seurat, reduction = "umap",pt.size =1)
