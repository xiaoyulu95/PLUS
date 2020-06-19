
#############################
#loading data
###########################
MAT<-as.matrix(read.table("GSE75688_malignant_cells.txt"))
MAT_seurat<-CreateSeuratObject(counts = MAT, project = "MAT", min.cells = 3, min.features = 200)
names(attributes(MAT_seurat))
dim(attributes(MAT_seurat)$assays$RNA)
MAT_seurat <- NormalizeData(MAT_seurat)
all.genes <- rownames(MAT_seurat)
MAT_seurat <- ScaleData(MAT_seurat, features = all.genes)

###################################
MAT_seurat <- RunPCA(MAT_seurat, features = all.var.genes)

MAT_seurat <- FindNeighbors(MAT_seurat, dims = 1:12)
MAT_seurat <- FindClusters(MAT_seurat, resolution = 0.9) #higher resulution gives bigger cluster number
  
MAT_seurat <- RunUMAP(MAT_seurat, dims = 1:12,perplexity = 20)
DimPlot(MAT_seurat, reduction = "umap", label=T)

Idents(object = MAT_seurat) <-   as.factor(MM1)
DimPlot(MAT_seurat,pt.size=1.5,label = FALSE)

Idents(object = MAT_seurat) <- PATIENT
DimPlot(MAT_seurat,pt.size=1.5,label = FALSE)

#######
MAT_seurat <- RunPCA(MAT_seurat, features = meta.var.genes)

MAT_seurat <- FindNeighbors(MAT_seurat, dims = 1:12) 
MAT_seurat <- FindClusters(MAT_seurat, resolution = 0.9) #higher resulution gives bigger cluster number
   
MAT_seurat <- RunUMAP(MAT_seurat, dims = 1:12,perplexity = 20)
DimPlot(MAT_seurat, reduction = "umap", label=T)

Idents(object = MAT_seurat) <-   as.factor(MM1)
DimPlot(MAT_seurat,pt.size=1.5,label = FALSE)

Idents(object = MAT_seurat) <- PATIENT
DimPlot(MAT_seurat,pt.size=1.5,label = FALSE)




