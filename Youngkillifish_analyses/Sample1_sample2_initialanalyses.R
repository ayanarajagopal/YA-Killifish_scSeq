####Sample 1 initial run####
datadir.young_new1 <- '~/rstudio_jan2020/20jan_afterintegration/Young_old_separate/Young_Isoseq_latest/July2020/18thjuly/filtered_matrix_july2021'
list.files(datadir.young_new1)
young_isoseq_matrix_new1 <- Read10X(datadir.young_new1, unique.features = T)
View(rownames(young_isoseq_matrix_new))
rownames(young_isoseq_matrix_new) <- toupper(rownames(young_isoseq_matrix_new))
Ribo.young_new1 <- grep(pattern = "^RP[SL]*", x = rownames(young_isoseq_matrix_new1), value = FALSE) # Select row indices and not ERCC names 
Ribo.young_new1
percent.Ribo.Young.new <- Matrix::colSums(young_isoseq_matrix_new[Ribo.young_new, ])/Matrix::colSums(young_isoseq_matrix_new1)
percent.Ribo.Young.new1
young_isoseq_matrix_new1 <- young_isoseq_matrix_new1[-Ribo.young_new1, ]
young_isoseq_matrix_new1

youngisoseq$log10GenesPerUMI <- log10(youngisoseq$nFeature_RNA) / log10(youngisoseq$nCount_RNA)

#set up young object
youngisoseq1 <- CreateSeuratObject(counts=young_isoseq_matrix_new1, project= "Young_KILLI_isoseq1", meta.data = data.frame(percent.ribo = percent.Ribo.Young.new1), min.cells=2)
dim(youngisoseq1)
table(Idents(youngisoseq1))

Mito.young_new1 <- grep(pattern = "^MT-", x = rownames(young_isoseq_matrix_new1), value = F)
Mito.young_new1
C<-GetAssayData(object = youngisoseq1, slot = "counts")
percent.mito <- Matrix::colSums(C[Mito.young_new1, ])/Matrix::colSums(C)
percent.mito
youngisoseq <- AddMetaData(youngisoseq, percent.mito, col.name = "percent.mt")

youngisoseq <- AddMetaData(youngisoseq, percent.Ribo.Young.new1, col.name = "percent.ribo")
youngisoseq <- AddMetaData(youngisoseq, youngisoseq$nCount_RNA, col.name = "nUMI")
youngisoseq <- AddMetaData(youngisoseq, youngisoseq$nFeature_RNA, col.name = "nGene")
VlnPlot(youngisoseq, features = c("percent.mt","percent.ribo", "nGene","nUMI"), group.by="orig.ident", ncol = 4)

grep ("^MT-", rownames(youngisoseq1[["RNA"]]),value = T)

par(mfrow = c(1, 2))
FeatureScatter(object = youngisoseq,group.by="orig.ident", feature1 = "nUMI", feature2 = "percent.mt")
FeatureScatter(object = youngisoseq,group.by="orig.ident", feature1 = "percent.ribo", feature2 = "nGene")
FeatureScatter(object = youngisoseq,group.by="orig.ident", feature1 = "nUMI", feature2 = "nGene")

#Filter bad quality cells
youngisoseq <- subset(youngisoseq, subset=nFeature_RNA > 200 & percent.mt < 5)
youngisoseq <- SCTransform(youngisoseq,verbose = T)
dim(youngisoseq)
youngisoseq <- RunPCA(object = youngisoseq, verbose = FALSE)
ElbowPlot(youngisoseq, ndims=40)

#Determine the K-nearest neighbor graph
youngisoseq <- FindNeighbors(object = youngisoseq, 
                             dims = 1:20)

# Determine the clusters for various resolutions                                
youngisoseq <- FindClusters(object = youngisoseq,
                            resolution = 0.5)

#Calc. of TSNE
youngisoseq <- RunTSNE(object = youngisoseq, perplexity=30, reduction= "pca",dims = 1:20)
youngisoseq <- RunUMAP(object = youngisoseq, reduction= "pca", dims = 1:20)

DimPlot(youngisoseq, reduction = "tsne", pt.size = 1, label=T)

####Sample 2 initial run####
datadir.young_sample2 <- 'C:/Users/u0129074/Documents/Killi_injury_ScSeq/6w_killifish_sample2/outs/filtered_feature_bc_matrix'
list.files(datadir.young_sample2)
young_matrix_S2 <- Read10X(datadir.young_sample2, unique.features = T)

rownames(young_matrix_S2) <- toupper(rownames(young_matrix_S2))
Ribo.young_S2 <- grep(pattern = "^RP[SL]*", x = rownames(young_matrix_S2), value = FALSE) # Select row indices and not ERCC names 
Ribo.young_S2
percent.Ribo.Young.S2 <- Matrix::colSums(young_matrix_S2[Ribo.young_S2, ])/Matrix::colSums(young_matrix_S2)
percent.Ribo.Young.S2
young_matrix_S2 <- young_matrix_S2[-Ribo.young_S2, ]
young_matrix_S2


#set up Seurat object
youngisoseq_S2 <- CreateSeuratObject(counts=young_matrix_S2, project= "Young_KILLI_isoseq_S2", meta.data = data.frame(percent.ribo = percent.Ribo.Young.S2), min.cells=3, min.features=200)
dim(youngisoseq_S2)

Mito.young_S2 <- grep(pattern = "^MT-", x = rownames(young_matrix_S2), value = F)
Mito.young_S2
C_S2 <-GetAssayData(object = youngisoseq_S2, slot = "counts")
percent.mito_S2 <- Matrix::colSums(C_S2[Mito.young_S2, ])/Matrix::colSums(C_S2)
percent.mito_S2

View(rownames(youngisoseq_S2))
rb.genes_S2 <- rownames(youngisoseq_S2)[grep("^RP[SL]",rownames(youngisoseq_S2))]

youngisoseq_S2 <- AddMetaData(youngisoseq_S2, percent.mito_S2, col.name = "percent.mt")
youngisoseq_S2 <- AddMetaData(youngisoseq_S2, percent.Ribo.Young.S2, col.name = "percent.ribo")
youngisoseq_S2 <- AddMetaData(youngisoseq_S2, youngisoseq_S2$nCount_RNA, col.name = "nUMI")
youngisoseq_S2 <- AddMetaData(youngisoseq_S2, youngisoseq_S2$nFeature_RNA, col.name = "nGene")

VlnPlot(youngisoseq_S2, features = c("percent.mt","percent.ribo", "nGene","nUMI"), group.by="orig.ident", ncol = 4)

View(youngisoseq_S2@meta.data)
#young_S2_obj <- saveRDS(youngisoseq_S2, file = "C:/Users/u0129074/Documents/ScRNAseq_paper/youngisoseq_S2.rds")

youngisoseq_S2 <- subset(youngisoseq_S2, subset=nFeature_RNA > 200 & percent.mt < 5)
youngisoseq_S2 <- SCTransform(youngisoseq_S2, verbose = T)
dim(youngisoseq_S2)
youngisoseq_S2 <- RunPCA(object = youngisoseq_S2, verbose = FALSE)

ElbowPlot(youngisoseq_S2, ndims=40)

# Determine the K-nearest neighbor graph
youngisoseq_S2 <- FindNeighbors(object = youngisoseq_S2, 
                                dims = 1:30)

# Determine the clusters for various resolutions                                
youngisoseq_S2 <- FindClusters(object = youngisoseq_S2,
                               resolution = 1)

#Calc. of TSNE
youngisoseq_S2 <- RunTSNE(object = youngisoseq_S2, perplexity=35, reduction= "pca",dims = 1:30)
youngisoseq_S2 <- RunUMAP(object = youngisoseq_S2, reduction= "pca", dims = 1:30)

DimPlot(youngisoseq_S2, reduction = "tsne", pt.size = 1, label=T)
