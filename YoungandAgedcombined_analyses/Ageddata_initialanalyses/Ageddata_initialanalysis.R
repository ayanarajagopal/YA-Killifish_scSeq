library(swne)

agedisoseq <- CreateSeuratObject(counts=aged_isoseq_matrix_new, project= "aged_KILLI_isoseq", meta.data = data.frame(percent.ribo = percent.Ribo.aged.new), min.cells=2)
?CreateSeuratObject

if(!require(remotes)){ install.packages("remotes") }
remotes::install_github("linxihui/NNLM")
remotes::install_github("yanwu2014/swne")

Mito.aged_new <- grep(pattern = "^MT-", x = rownames(aged_isoseq_matrix_new), value =F )
Mito.aged_new
C<-GetAssayData(object = agedisoseq, slot = "counts")
percent.mito.aged.new <- Matrix::colSums(C[Mito.aged_new, ])/Matrix::colSums(C)
percent.mito.aged.new
agedisoseq <- AddMetaData(agedisoseq, percent.mito.aged.new, col.name = "percent.mt")
agedisoseq <- AddMetaData(agedisoseq, percent.Ribo.aged.new, col.name = "percent.ribo")
agedisoseq <- AddMetaData(agedisoseq, agedisoseq$nCount_RNA, col.name = "nUMI")
agedisoseq <- AddMetaData(agedisoseq, agedisoseq$nFeature_RNA, col.name = "nGene")
#youngisoseq <- AddMetaData(youngisoseq, percent.Mito.Young.new, col.name = "percent.mito")


DE_SASP_vs_all <- FindMarkers(agedisoseq, ident.1 = "SASP cells", ident.2="other", logfc.threshold = 0, test.use = "MAST", only.pos = FALSE, min.pct = 0.0)
dim(agedisoseq)
#agedisoseq[["percent.mt"]] <- PercentageFeatureSet(agedisoseq, pattern="^MT-")
VlnPlot(agedisoseq, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)

agedisoseq <- subset(agedisoseq, subset=nFeature_RNA > 200 & percent.mt < 5)
agedisoseq <- SCTransform(agedisoseq, verbose = T)

agedisoseq <- RunPCA(object = agedisoseq, verbose = FALSE)
ElbowPlot(agedisoseq, ndims=30)

# Determine the K-nearest neighbor graph
agedisoseq <- FindNeighbors(object = agedisoseq, 
                            dims = 1:20)

# Determine the clusters for various resolutions                                
agedisoseq <- FindClusters(object = agedisoseq,
                           resolution = 0.8)

clusters <- agedisoseq$seurat_clusters
clusters
#Calc. of TSNE
agedisoseq <- RunTSNE(object = agedisoseq, perplexity=90, reduction= "pca",dims = 1:20)
agedisoseq <- RunUMAP(object = agedisoseq, reduction= "pca", dims = 1:20)
## Run SWNE
norm.counts.a <- ExtractNormCounts(agedisoseq, obj.type = "seurat", rescale.method = "log")
dim(norm.counts.a)
var.genes.a <- VariableFeatures(agedisoseq)


FindNumFactors(agedisoseq,k.range=seq(2, 30, 2), do.plot=T)
k.err.a <- FindNumFactors(norm.counts.a[var.genes.a,], k.range = seq(2, 30, 2), n.cores = 8, do.plot = T)
k.err.a
?FindNumFactors
swne_embed <- RunSWNE(agedisoseq, k = 24, alpha.exp=1.7, genes.embed = JVH_SASPs)
## Plot SWNE
PlotSWNE(swne_embed, alpha.plot = 0.4, sample.groups = clusters, do.label=T, label.size = 3.5, pt.size = 1.5, show.legend = F, seed = 42)

genes.to.label.clusters <- c("CLU","EPD","APOA1","APOEB","AIF1","PFN1","LCP1","RLBP1","HSPB1","BCAN","CNDP1","WNT8B","HEPACAM","S100B","GLUL (2 OF 2)","SLC1A2","PCNA","STMN1A","HMGB2A","CD9B","OLIG1","CX43","FABP7A","TAGLN","PDGFRB","KLF2","SOX17","STAB2","FLT4","ELAVL4","SLC17A6B","SV2A","MAP2","SYPA","SYP","GABRA1","CCK(1OF2)","CAMK2N2","CAMK2B1","ELAVL3","MAP1B","NEUROD2","DCX","TUBB5")


?PlotSWNE
#killi.integrated <- FindNeighbors(killi.integrated, reduction = "pca", dims = 1:15)

DimPlot(agedisoseq, reduction = "swne", group.by="orig.ident",pt.size = 1.2)
#DimPlot(agedisoseq, reduction = "tsne",label=T,cols=col_new1, pt.size = 1.1,label.size = 5, repel=T) + theme_test() 
DimPlot(agedisoseq, reduction = "tsne",label=T, pt.size = 1.2,label.size = 5) + theme_test()
?SCTransform
DefaultAssay(agedisoseq) = "RNA"
DefaultAssay(agedisoseq)
table(Idents(agedisoseq))

get_gene_names(aged_isoseq_matrix_new)

agedisoseq <- NormalizeData(agedisoseq, verbose = FALSE)

agedisoseq_markers <- FindAllMarkers(agedisoseq, min.pct = 0.25, only.pos = F, logfc.threshold = 0.25)

allgenes <-FindAllMarkers(agedisoseq, min.pct = 0.0, only.pos = F, logfc.threshold = 0.0)

agedisoseq_markers
aged_isoseq_cluster_top <- agedisoseq_markers %>% group_by(cluster) %>% top_n(n=500, wt = avg_log2FC)
write.table(aged_isoseq_cluster_top, file="C:/Users/u0129074/Documents/ScRNAseq_paper/YoungvsAged/march2022/Allcells_1.0_20pcs_fin.txt", sep="\t", col.names=TRUE)

agedisoseq@meta.data

agedisoseq.Ncbi.res0.8.new <- c("mN1","mN2","Trans-N","MG1","InterC","mN3","mN4","Ih-mN5","NGP","Vas1","Vas2","Astro-RG1","EPD-RG4","OPC/OD","Astro-RG2","NE-RG3","MG2","MG3","MG4","NBN","MG5")
agedisoseq.Ncbi.res0.8.fresh <- c("Neuronal1","Neuronal2","Trans-N","MG1","InterC","Neuronal3","Neuronal4","Neuronal5","NGP","Vas1","Vas2","Glial1","Glial4","OPC/OD","Glial4","Glial3","MG2","MG3","MG4","NBN","MG5")

names(agedisoseq.Ncbi.res0.8.fresh) <- levels(agedisoseq)
agedisoseq <- RenameIdents(agedisoseq, agedisoseq.Ncbi.res0.8.fresh)

agedisoseq$CellType <- Idents(agedisoseq)
agedisoseq$CellType
agedisoseq$seurat_clusters <- Idents(agedisoseq)
agedisoseq$seurat_clusters

table(Idents(agedisoseq))

write.table(agedisoseq@assays$RNA@counts, file="C:/Users/u0129074/Documents/ScRNAseq_paper/YoungvsAged/march2022/Allgenenames.txt", sep="\t", col.names=TRUE)

#Subsetting the microglial clusters
MGs_subset <- subset(agedisoseq, idents=c("Trans-N","MG1","MG2","MG3","MG4","MG5"))

View(rownames(MGs_subset))
MGs_subset <- SCTransform(MGs_subset, verbose = TRUE )
MGs_subset <- RunPCA(MGs_subset, verbose=FALSE)
ElbowPlot(MGs_subset, ndims = 30)
MGs_subset<- FindNeighbors(object = MGs_subset, 
                             dims = 1:20)
youngisoseq$seurat_clusters
# Determine the clusters for various resolutions                                
MGs_subset <- FindClusters(object = MGs_subset,
                             resolution = 0.5)
#Idents(object = killi.integrated) <- "RNA_snn_res.0.8"

#Calc. of TSNE
MGs_subset <- RunTSNE(object = MGs_subset, reduction= "pca", dims=1:20)
MGs_subset <- RunUMAP(object = MGs_subset, reduction= "pca", dims = 1:20)

saveRDS(MGs_subset, "~/ScRNAseq_paper/Redo_november2020/june2021/MGsubsetting_new.Rds") 

DimPlot(MGs_subset, reduction = "tsne", group.by="orig.ident",pt.size = 1.2)
#DimPlot(agedisoseq, reduction = "tsne",label=T,cols=col_new1, pt.size = 1.1,label.size = 5, repel=T) + theme_test() 
DimPlot(MGs_subset, reduction = "tsne",label=T, pt.size = 1.2,label.size = 5) + theme_test()
DefaultAssay(MGs_subset) = "RNA"

table(Idents(MGs_subset))

MGs_subset <- NormalizeData(MGs_subset, verbose = FALSE)

MGs_subset_markers <- FindAllMarkers(MGs_subset, min.pct = 0.25, only.pos = F, logfc.threshold = 0.25)
MGs_subset_markers
MGs_subset_cluster_top <- MGs_subset_markers %>% group_by(cluster) %>% top_n(n=500, wt = avg_log2FC)
write.table(MGs_subset_cluster_top, file="C:/Users/u0129074/Documents/ScRNAseq_paper/YoungvsAged/march2022/MGsubclusters_transN_0.5_20pcs.txt", sep="\t", col.names=TRUE)

agedisoseq@meta.data

MGs.Ncbi.res0.5.new <- c("mN1","mN2","Trans-N","MG1","InterC","mN3","mN4","Ih-mN5","NGP","Vas1","Vas2","Astro-RG1","EPD-RG4","OPC/OD","Astro-RG2","NE-RG3","MG2","MG3","MG4","NBN","MG5")

MGsonly_idents = c("Trans-N","mN","MG1","MG2","MG3","MG4","MG5")

names(MGsonly_idents) <- levels(MGs_subset)
MGs_subset <- RenameIdents(MGs_subset, MGsonly_idents)

MGs_subset$CellType <- Idents(MGs_subset)
MGs_subset$CellType
MGs_subset$seurat_clusters <- Idents(MGs_subset)
MGs_subset$seurat_clusters
table(Idents(MGs_subset))


names(agedisoseq.Ncbi.res0.8.new) <- levels(agedisoseq)
agedisoseq <- RenameIdents(agedisoseq, agedisoseq.Ncbi.res0.8.new)

agedisoseq$CellType <- Idents(agedisoseq)
agedisoseq$CellType
agedisoseq$seurat_clusters <- Idents(agedisoseq)
agedisoseq$seurat_clusters

table(Idents(agedisoseq))

#Subsetting the microglial clusters
MGs_subset <- subset(agedisoseq, idents=c("MG1","MG2","MG3","MG4","MG5"))
MGs_subset
View(rownames(MGs_subset))
MGs_subset <- SCTransform(MGs_subset, verbose = TRUE )
MGs_subset <- RunPCA(MGs_subset, verbose=FALSE)
ElbowPlot(MGs_subset, ndims = 30)
MGs_subset<- FindNeighbors(object = MGs_subset, 
                           dims = 1:20)
youngisoseq$seurat_clusters
# Determine the clusters for various resolutions                                
MGs_subset <- FindClusters(object = MGs_subset,
                           resolution = 0.5)
#Idents(object = killi.integrated) <- "RNA_snn_res.0.8"

#Calc. of TSNE
MGs_subset <- RunTSNE(object = MGs_subset, reduction= "pca", dims=1:20)
MGs_subset <- RunUMAP(object = MGs_subset, reduction= "pca", dims = 1:20)

saveRDS(MGs_subset, "~/ScRNAseq_paper/Redo_november2020/june2021/MGsubsetting.Rds") 

DimPlot(MGs_subset, reduction = "tsne", group.by="orig.ident",pt.size = 1.2)
#DimPlot(agedisoseq, reduction = "tsne",label=T,cols=col_new1, pt.size = 1.1,label.size = 5, repel=T)
DimPlot(MGs_subset, reduction = "tsne",label=T, pt.size = 1.5,cols=c("#555555","gray","red","#ff3de1","#39FF14"))+ theme_test()
table(Idents(MGs_subset))
DefaultAssay(MGs_subset) = "RNA"

FeaturePlot(agedisoseq, reduction = "tsne", features="cluster_id_2", label=TRUE, label.size=3.5) +theme(axis.title.x = element_blank(), axis.title.y = element_blank())& NoAxes()

table(Idents(MGs_subset))
A <- FeaturePlot(agedisoseq, reduction = "swne", features="SEN_Mayo", label=TRUE, label.size=3.5)+
  scale_colour_gradient2(midpoint=0.2, low="blue", mid="grey", high="red", space="Lab")+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())& NoAxes()
B <- FeaturePlot(agedisoseq, reduction = "umap", features="Inhouse_SASPS", label=TRUE, label.size=3.5)+
  scale_colour_gradient2(midpoint=0.4, low="blue", mid="grey", high="red", space="Lab")+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())& NoAxes()
B1 <- VlnPlot(agedisoseq, "MIF", log=TRUE)+theme(axis.title.x = element_blank()) + ylab("Log UMI")& NoLegend()
B2 <- VlnPlot(agedisoseq, "IL1B", log=TRUE)+theme(axis.title.x = element_blank()) + ylab("Log UMI")& NoLegend()

top5_markers_MG <- Extract_Top_Markers(marker_dataframe = MGs_subset_markers, num_genes = 10, named_vector = FALSE,
                                    make_unique = TRUE)
Clustered_DotPlot(seurat_object = MGs_subset, features = top5_markers_MG, k=5)
?Clustered_DotPlot
print(A | B | B1/B2)
dim(MGs_subset)
MGs_subset <- NormalizeData(MGs_subset, verbose = FALSE)

MGs_subset_markers <- FindAllMarkers(MGs_subset, min.pct = 0.25, only.pos = F, logfc.threshold = 0.25)
MGs_subset_markers
MGs_subset_cluster_top <- MGs_subset_markers %>% group_by(cluster) %>% top_n(n=500, wt = avg_log2FC)
write.table(MGs_subset_cluster_top, file="C:/Users/u0129074/Documents/ScRNAseq_paper/YoungvsAged/march2022/MGsubclusters_0.5_20pcs_fin.txt", sep="\t", col.names=TRUE)

MGs_top5 <- MGs_subset_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
head(MGs_top5)

Jolien_saspfactors = c("TP53","MDM2","CDKN1BA","CDKN1B","CXCL8","MAP2K6","RRAS2","BRF1B","AREGB","IGFBP3","IGFBP4","IGFBP5","PRODH (2 OF 2)","SMURF2","PGF","CTNNB1","ANGPTL2","LMNB1","FAS")
DotPlot(agedisoseq,cols=c("green","purple"),split.by ="sample",features= c("GPR17","SOX10")) + RotatedAxis()

VlnPlot(MGs_subset,features=Jolien_saspfactors,cols=c("green","purple","dark blue","yellow","grey","blue","steelblue"))
DotPlot(MGs_subset,cols=c("green","purple"), features=c("AIF1","C1QA","C1QB","CD68","CD3E","HBA1","ITGAL","CORO1A","ARPC1B","APOEB","CD74","PFN1","MCP1","CCL3","CCL4","MRC1")) + RotatedAxis()
VlnPlot(MGs_subset,cols=c("green","purple","dark blue","yellow","grey","blue","steelblue"), features=c("AIF1","C1QA","C1QB","CD68","CD3E","ITGAL","CORO1A","ARPC1B","APOEB","CD74","PFN1","MRC1","LCP1")) + RotatedAxis()
VlnPlot(MGs_subset,cols=c("green","purple","dark blue","yellow","grey","blue","steelblue"), features=c("PTPRC","ITGAM","AIF1","C1QA","C1QB","CTSS","CD14","CSFR3","ARGLU1","FAM46A","ISG15 (2 OF 2)","IFIT3","MRC1","TNFA","CD83","EGR2")) + RotatedAxis()
VlnPlot(MGs_subset,cols=c("green","purple","dark blue","yellow","grey","blue","steelblue"), features=c("TNFSF18","CCL8","TFRC (2 OF 2)","FCGBP","GPR84","PCNA","CDC20","BIRC5A","FCN1","VCANA","CD3E","GZMB","IL7R","FGFBP2 (1 OF 2)","CCR7","CD79A","FABP7A","MBPB","SNAP25 (1 OF 2)","HBAA1","TFA","CD68","CTSBA","CD74","PFN1","SNCB","CD63","GRIA2","NRXN1A")) + RotatedAxis()

#Regulators of TFs
VlnPlot(MGs_subset, cols=c("green","purple","dark blue","yellow","grey","blue","steelblue"), features=c("KLF2","ZNF821","OAS2","OAS3","GRHL1","NR4A1","KLF10","EGR2","PCNA","ORC6","FEN1","RFX2")) + RotatedAxis()

#Cell surface/membrane genes
VlnPlot(MGs_subset, cols=c("green","purple","dark blue","yellow","grey","blue","steelblue"), features=c("LY6E","GPR85","IFITM5","IRAK2","FLT1","PLXNA2","CD68","EMC9","SLC16A1")) + RotatedAxis()

#DAM microglia
VlnPlot(MGs_subset, cols=c("green","purple","dark blue","yellow","grey","blue","steelblue"), features=c("CST7","APOEB","LPL","SPP2","ITGAX","CSF1","CSF1RA","","CLEC7A","IGF1","AXL","CD63","CTSD","TYROBP","CTSB","TREM2","CTSL","CD9","CTSZ","P2RY12","TMEM119","P2RY13","SELPLG","TXNIP","CCR5")) + RotatedAxis()
DotPlot(MGs_subset, cols=c("green","steelblue"), features=c("LOC107396233","APOEB","LPL","SPP2","ITGAX","CORO1A","CSF1RA","CLEC7A","IGF1","AXL","CD63","CTSD","TYROBP","CTSB","TREM2","CTSL","CD9","CTSZ","P2RY12","TMEM119","P2RY13","SELPLG","TXNIP","CCR5")) + RotatedAxis()

#Early response MG
VlnPlot(MGs_subset, cols=c("green","purple","dark blue","yellow","grey","blue","steelblue"), features=c("TOP2A","RRM2","HELLS","GMNN","SPC25","NUSAP","KNSTRN","KIF22","AURKB","C3","C4B","CFB","AXL","CLEC7A","LGALS3","CD74","IFITM3","IRF7","RSAD2")) + RotatedAxis()
DotPlot(MGs_subset, cols=c("green","purple"), features=c("TOP2A","RRM2B","HELLS","GMNN","SPC25","NUSAP1","KNSTRN","KIF22","AURKB","C3","C4","CFB","AXL","CLEC7A","LGALS3","CD74","IFITM3","IRF7","RSAD2")) + RotatedAxis()

cc.genes.updated.2019
#cell cycle analysis
cellcycle_MG <- CellCycleScoring(MGs_subset, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes) 
cellcycle_MG[[]]
as_tibble(cellcycle_MG[[]]) %>%
  ggplot(aes(Phase)) + geom_bar()
cellcycle_MG@meta.data %>%
  group_by(CellType,Phase) %>%
  count() %>%
  group_by(CellType) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=CellType,y=percent, split,fill=Phase)) +
  geom_col() + RotatedAxis()+
  ggtitle("Percentage of cell cycle phases per cluster")

#Cellular senescence analysis

SenMayoGS <- c()

#aged PCs sensecence data analysis
afish_lineage_new
clusters_aPCs <- afish_lineage_new$seurat_clusters
swne_embed.aPCs <- RunSWNE(afish_lineage_new, k = 16, alpha.exp=1.4, genes.embed = MGs_top5)
## Plot SWNE
colr_aPC<- ExtractSWNEColors(swne_embed.aPCs, clusters_aPCs, seed=42)
colr_aPC <- colr_PC
PlotSWNE(swne_embed.aPCs, alpha.plot = 0.4, colors.use=, sample.groups = clusters_aPCs, do.label=T, label.size = 3.5, pt.size = 1.5, show.legend = F, seed = 42)


#onlyMGs
clusters_MGs <- MGs_subset$seurat_clusters
var.genes.MGs <- VariableFeatures(MGs_subset)
View(MGs_top5)
mgsmarkers_list <- MGs_top5 <-
DefaultAssay(MGs_subset) <- "SCT"
MGmarkers_1_TM <- c("AIF1","C1QA","C1QB","CD68","CD3E","HBA1","ITGAL","CORO1A","ARPC1B","APOEB","CD74","PFN1","MCP1","CCL3","CCL4","MRC1")
MGmarkers_2_DAM <- c("CST7","APOEB","LPL","SPP2","ITGAX","CSF1","CSF1RA","CLEC7A","IGF1","AXL","CD63","CTSD","TYROBP","CTSB","TREM2","CTSL","CD9","CTSZ","P2RY12","TMEM119","P2RY13","SELPLG","TXNIP","CCR5")
MGmarkers_3_ER <- c("TOP2A","RRM2","HELLS","GMNN","SPC25","NUSAP","KNSTRN","KIF22","AURKB","C3","C4B","CFB","AXL","CLEC7A","LGALS3","CD74","IFITM3","IRF7","RSAD2")
MGmarkers_4_TFr <- c("KLF2","ZNF821","OAS2","OAS3","GRHL1","NR4A1","KLF10","EGR2","PCNA","ORC6","FEN1","RFX2")
swne.MGs <- RunSWNE(MGs_subset, k = 16, alpha.exp=1.4, genes.embed = MGmarkers_2_DAM)
PlotSWNE(swne.MGs, alpha.plot = 0.4, sample.groups = clusters_MGs, do.label=T, label.size = 3.5, pt.size = 1.5, show.legend = F, seed = 42)

gene.use <- "PCNA"
norm.counts.agedMGs <- ExtractNormCounts(MGs_subset, obj.type = "seurat", rescale.method = "log")

gene.expr <- norm.counts.agedMGs[gene.use,]
FeaturePlotSWNE(swne.MGs, gene.expr, gene.use, alpha.plot = 0.4, label.size = 3.5, pt.size = 1.25)

table(Idents(MGs_subset))
CellCycleScoring(agedisoseq, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE) -> data
data[[]]
as_tibble(data[[]]) %>%
  ggplot(aes(Phase)) + geom_bar()
data@meta.data %>%
  group_by(CellType,Phase) %>%
  count() %>%
  group_by(CellType) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=CellType,y=percent, split,fill=Phase)) +
  geom_col() + RotatedAxis()+
  ggtitle("Percentage of cell cycle phases per cluster")

