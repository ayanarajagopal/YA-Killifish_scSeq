BiocManager::install("escape")
BiocManager::install("dittoSeq")
remotes::install_github("crazyhottommy/scclusteval")
library(scclusteval)

datadir.aged <- '~/rstudio_jan2020/20jan_afterintegration/Young_old_separate/18w_killifish_18july/outs/filtered_feature_bc_matrix/'
list.files(datadir.aged) 
aged_isoseq_matrix_new <- Read10X(data.dir=datadir.aged, unique.features = T)
View(rownames(aged_isoseq_matrix_new))
rownames(aged_isoseq_matrix_new) <- toupper(rownames(aged_isoseq_matrix_new))
View(rownames(aged_isoseq_matrix_new))
Ribo.aged_new <- grep(pattern = "^RP[SL]*", x = rownames(aged_isoseq_matrix_new), value = FALSE) # Select row indices and not ERCC names 
Ribo.aged_new
percent.Ribo.aged.new <- Matrix::colSums(aged_isoseq_matrix_new[Ribo.aged_new, ])/Matrix::colSums(aged_isoseq_matrix_new)
percent.Ribo.aged.new
aged_isoseq_matrix_new <- aged_isoseq_matrix_new[-Ribo.aged_new, ]
aged_isoseq_matrix_new

#set up aged object
agedisoseq <- CreateSeuratObject(counts=aged_isoseq_matrix_new, project= "aged_KILLI_isoseq", meta.data = data.frame(percent.ribo = percent.Ribo.aged.new), min.cells=2)
?CreateSeuratObject

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
                            resolution = 0.5)

#Calc. of TSNE
agedisoseq <- RunTSNE(object = agedisoseq, reduction= "pca",dims = 1:20)
agedisoseq <- RunUMAP(object = agedisoseq, reduction= "pca", dims = 1:20)

#killi.integrated <- FindNeighbors(killi.integrated, reduction = "pca", dims = 1:15)

DimPlot(youngisoseq, reduction = "tsne", group.by="orig.ident",pt.size = 1)
DimPlot(agedisoseq, reduction = "tsne",label=T,cols=col_new1, pt.size = 1.1,label.size = 5, repel=T) + theme_test() 
DimPlot(agedisoseq, reduction = "tsne",label=T, pt.size = 1.1,label.size = 5) + theme_test()
?SCTransform
DefaultAssay(agedisoseq) = "RNA"
youngisoseq
dim(agedisoseq)
dim(youngisoseq)

raw.data_YandA = list("YoungK" =youngisoseq,"AgedK"=agedisoseq)
agedisoseq$sample <- "AgedK"
youngisoseq$sample <- "YoungK"

combined.features_YandA <- SelectIntegrationFeatures(object.list = raw.data_YandA, nfeatures = 3000)
raw.data_YandA <- PrepSCTIntegration(object.list = raw.data_YandA, anchor.features = combined.features_YandA)

#referencing the ref dataset: zfish
#reference_dataset <- which(names(raw.data) == "Young_zebrafish")
#names(raw.data)

int.anchors_YandA <- FindIntegrationAnchors(object.list = raw.data_YandA, 
                                      anchor.features = combined.features_YandA, reduction="cca", normalization.method = "SCT")
View(int.anchors_YandA@anchor.features)

killi_youngandaged_integrated <- IntegrateData(anchorset = int.anchors_YandA, normalization.method = "SCT", dims = 1:30)

DefaultAssay(killi_youngandaged_integrated) <- "integrated"

View((killi_youngandaged_integrated))
killi_youngandaged_integrated <- ScaleData(killi_youngandaged_integrated, verbose = FALSE)
killi_youngandaged_integrated <- RunPCA(object = killi_youngandaged_integrated, verbose = FALSE)
ElbowPlot(killi_youngandaged_integrated, ndims=40)

# Determine the K-nearest neighbor graph
killi_youngandaged_integrated <- FindNeighbors(object = killi_youngandaged_integrated, 
                                       dims = 1:20)

# Determine the clusters for various resolutions                                
killi_youngandaged_integrated <- FindClusters(object = killi_youngandaged_integrated,
                                      resolution = 0.6)
#Idents(object = killi.integrated) <- "RNA_snn_res.0.8"

#Calc. of TSNE
killi_youngandaged_integrated <- RunUMAP(object = killi_youngandaged_integrated, reduction= "pca", dims = 1:20)

DimPlot(killi_youngandaged_integrated, reduction = "tsne", group.by="orig.ident",pt.size = 0.3)
DimPlot(killi_youngandaged_integrated, reduction = "tsne", pt.size = 1.0, split.by="sample",label = T,label.size=3) + theme_test()
plot_grid(p1, p2)
dev.off()
?DimPlot
#normalize and setting RNA assay
DefaultAssay(killi_youngandaged_integrated) <- "RNA"
killi_youngandaged_integrated <- NormalizeData(killi_youngandaged_integrated, verbose = T)

#expression_table
youngandaged_isoseq_combinedclusters <- FindAllMarkers(killi_youngandaged_integrated, min.pct = 0.25, only.pos = F, logfc.threshold = 0.25, test.use="MAST")
head(youngandaged_isoseq_combinedclusters)
young_aged_isoseq_markers <- youngandaged_isoseq_combinedclusters %>% group_by(cluster) %>% top_n(n=1000, wt = avg_log2FC)
young_aged_top30_markers <- youngandaged_isoseq_combinedclusters %>% group_by(cluster) %>% top_n(n=30, wt = avg_log2FC)
head(young_aged_top30_markers)
write.table(young_aged_isoseq_markers, file="C:/Users/u0129074/Documents/ScRNAseq_paper/Redo_november2020/Aug2021/y_a_markers_res0.6_new.txt", sep="\t", col.names=TRUE)
write.table(young_aged_top30_markers, file="C:/Users/u0129074/Documents/ScRNAseq_paper/Redo_november2020/Aug2021/y_a_markers_top30_res1.txt", sep="\t", col.names=TRUE)

#Markers for proliferation clusters only
prol_clustermarkers <- FindConservedMarkers(killi_youngandaged_integrated, idents = c("NGP","RG1","RG2","RG3","RG4","Inter-NGP","Inter-RG"), grouping.var = sample)
prol_clustermarkers

RG1_changes <- read.csv(file="C:/Users/u0129074/Documents/ScRNAseq_paper/Redo_november2020/march2021/Analysis/rg1_top30.csv")
RG3_changes <- read.csv(file="C:/Users/u0129074/Documents/ScRNAseq_paper/Redo_november2020/march2021/GOanalysis/rg3_markers_yanda.csv")

#Regional markers: Cosacak et al
DefaultAssay(killi_youngandaged_integrated) <- "integrated"
PCs_newanno_1$
young_aged_isoseq_combinedclusters_top
scale.data = GetAssayData(object=killi_youngandaged_integrated[["SCT"]], slot="counts")
scale.data <- as.matrix(x = scale.data + 1)
killi_youngandaged_integrated <- SetAssayData(object=killi_youngandaged_integrated, slot="SCT", new.data= scale.data, assay="RNA")

top10 <- youngandaged_isoseq_combinedclusters %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(killi_youngandaged_integrated, features = top10$gene,   draw.lines = TRUE,
          lines.width = NULL,
          group.bar.height = 0.02) + scale_fill_gradientn(colors = c("green","black","magenta"))

youngfish.Ncbi.res0.4.ids <- c("ImN1","RG1","ImmC", "mN1", "mN2", "IPC","InterN","EC","NA","NGP","mN3","RG2","RG3","ImN2","RG4","OPC","PeriC")
youngfish.Ncbi.res0.6.ids <- c("ImN1","mN1","InterN","mN2","ImmC","RG1","mN3","ImN2","InterC","EC","PeriC","NGP","RG3","RG4","RG2","OPC","PVC")
YvsA.Ncbi.res0.6.n.ids <- c("Trans-N","Ex-mN1","Ex-mN2","NBN","Ex-mN3","MG","Inter-NGP","Inh-mN4","Inter-RG","Vas1","Vas2","EPD-RG4","NA","NE-RG3","Astro-RG1","OPC","NGP","Astro-RG2","Vas3")
YvsA.Ncbi.res0.6.n.ids_JVH <- c("ImN","MN","MN","ImN","MN","MG","Inter-NGP","MN","Inter-RG","Vas1","Vas2","EPD","NA","RG","RG","OD","NGP","RG","Vas3")

names(YvsA.Ncbi.res0.6.n.ids_JVH) <- levels(killi_youngandaged_integrated)
killi_youngandaged_integrated <- RenameIdents(killi_youngandaged_integrated, YvsA.Ncbi.res0.6.n.ids_JVH)

killi_youngandaged_integrated$CellType <- Idents(killi_youngandaged_integrated)
killi_youngandaged_integrated$CellType
killi_youngandaged_integrated$seurat_clusters <- Idents(killi_youngandaged_integrated)
killi_youngandaged_integrated$seurat_clusters


agedisoseq$CellType <- Idents(agedisoseq)
agedisoseq$CellType
agedisoseq$seurat_clusters <- Idents(agedisoseq)
agedisoseq$seurat_clusters


(killi_young_integrated)

##Cell type per cluster barcode
write.table(killi_youngandaged_integrated@active.ident, file='C:/Users/u0129074/Documents/ScRNAseq_paper/Redo_november2020/Cell_Cluster.tsv', quote=FALSE, sep='\t', col.names = FALSE)

plot(colSums(GetAssayData(killi_youngandaged_integrated,assay="integrated")>5))

proliferation_clusters <- subset(killi_youngandaged_integrated, idents = c("NGP","Astro-RG1","Astro-RG2","NE-RG3","EPD-RG4","Inter-NGP","Inter-RG"))
VlnPlot(proliferation_clusters, features = c("PCNA","FABP7A","GLUL (2 OF 2)","CX43","S100B","SLC1A2B","SOX2","SOX3","GFAP","HEPACAM","WNT8B","HES4.1","APOA1","NDRG4","S100A10A"), stack=T, split.by = "sample") + theme_test() + RotatedAxis()
DotPlot(proliferation_clusters, features = c("SOX2","SOX3","MEX3A","PCNA","GLUL (2 OF 2)","APOA1","EPD","CX43","S100B","FABP7A","SLC1A2B","GFAP","HEPACAM","WNT8B","HES4.1","NDRG4","S100A10A"), cols=c("steel blue","darkorange1"), split.by = "sample") + theme_few() + RotatedAxis() + coord_flip()
DotPlot(proliferation_clusters, features = c("SCARA3","S100A10A","FABP7A","NEUROG1","HMGB2A","PCNA","CX43","SLC1A3B","SLC1A2B","S100B","GLUL (2 OF 2)","HEPACAM","WNT8B","VIM","HES1","APOA1","EPD","ID4","HES5 (3 OF 9)","SOX4"), cols=c("steel blue","darkorange1"), scale = F, split.by = "sample") + theme_test() + RotatedAxis() + coord_flip()
DotPlot(proliferation_clusters, features = c("SOX2","GFAP","NOTCH3","HEY1","JUN","HES1","HES4","HES5 (2 OF 2)","VIM","NES","KDM6B","ASCL1B","LGALS2"), cols=c("steel blue","darkorange1"), split.by = "sample") + theme_few() + RotatedAxis() + coord_flip()


?DotPlot
proliferation_clusters
#Forisoseq paper
DotPlot(killi_youngandaged_integrated, features = c("GPR17","RGS9B","OPTN","TMEM5","TAB1","GPR85","FUS","HSPA4B","ARFIP2B","SCAF1","MAPRE3A","MAP7D2B","DCUN1D4","PRRT1","CAMK2N2","GAP43","TUBB5","REM2","LOC107384102","LOC107384764","MRPS15","ZC2HC1A","LOC107391117","UPF3A","MEF2D","TOM1L2 (1 OF 2)","KIRREL","TLDC1"), split.by = "sample") + theme_test() + RotatedAxis()
DotPlot(proliferation_clusters, features = c("RGS9B","OPTN","TMEM5","TAB1","GPR85","FUS","HSPA4B","ARFIP2B","SCAF1","MAPRE3A","MAP7D2B","DCUN1D4","PRRT1","CAMK2N2","GAP43","TUBB5","REM2","LOC107384102","LOC107384764","MRPS15","ZC2HC1A","LOC107391117","UPF3A","MEF2D","TOM1L2 (1 OF 2)","KIRREL","TLDC1"), split.by = "sample") + theme_test() + RotatedAxis()

DotPlot(neuronal_clusters, features = c("EIF1AXB","EIF3C","EIF3H","EIF4B","EIF4E3","EIF4EB","EIF4G3"), dot.scale = 6, split.by = "sample") + theme_test() + RotatedAxis()
DotPlot(killi_youngandaged_integrated, features = c("EIF1AXB","EIF3C","EIF3H","EIF4B","EIF4E3","EIF4EB","EIF4G3"),dot.scale = 4, split.by = "sample", ) + theme_test() + RotatedAxis()

DotPlot(oligs, features=c("OLIG1","OLIG2","SOX10","CD9B","SOX2","SOX8","GPR17"), cols=c("#FA340D","#005b96"), split.by="sample") + theme_test()
oligs <- subset(killi_youngandaged_integrated, idents = c("OD","RG"))

PCs_only <- subset(killi_youngandaged_integrated, idents= c("NGP","NE-RG3","Astro-RG1","Astro-RG2","EPD-RG4"))
DotPlot(proliferation_clusters, features = c("NEUROG1","NOTCH1A","NOTCH1B","DLA","DLD","ASCL1B","ID4","PCNA","KIAA0101","STMN1A","MKI67","HES5 (3 OF 9)","HES5(2OF2)","HES3","HES1","EOMESA","HMGA1","HMGN1","HMGN3","HMGB2A","NOTCH3","SCARA3.3","S100A10A"),cols = c("green", "purple"), split.by ="sample",cluster.idents = F) + RotatedAxis()
?ggthemes()
NC_clusters <- subset(killi_youngandaged_integrated, idents = c("NGP","ImN1","mN1","mN2","InterN","InterC"))
DotPlot(neuronal_clusters,features= c("TUBB5","DCX","STMN1B","ELAVL3","ELAVL4","CAMK2N2","CAMK2B","GABRA1","DIO3","GAD2","GAD1B","MAP2","SV2A","SOX2","PCNA","KIAA0101","HMGB2A","STMN1A","FABP7A","GLUL (2 OF 2)","CX43","VIM","S100B","SLC1A2","GFAP","EMX2","ID4","CLU","MKI67","EOMESA"), cols=c("green","purple"),split.by ="sample")+ RotatedAxis()

Neuronal_clusters <- subset(killi_youngandaged_integrated, idents = c("ImN1","ImN2","mN1","mN2","InterN","IPC"))
DotPlot(Neuronal_clusters,features= c("SLC17A8","SLC17A6B",,"PI4K2A","SNAP25 (1 OF 2)","OLA.13679","SNAP23","APPB","GPM6AB","FZD6","FZD3A","FZD3B","FZD2","SRGAP2A","TMEM125B","NYAP1","FMR1","SMO"), cols=c("green","purple"),split.by ="sample") + RotatedAxis()
VlnPlot(Neuronal_clusters,features= c("SLC17A6B","SNAP25 (1 OF 2)","SNAP25 (2 OF 2)","SNAP23","APPB","GPM6AB","SYT2 (1 OF 2)"), cols=c("green","purple"),split.by ="sample") + RotatedAxis()
VlnPlot(Neuronal_clusters,features= c("GAP43","SV2A","PROX1B","DCX","SLC7A11"), stack=T,cols=c("dark blue","turquoise"),adjust=T,split.by ="sample",flip=T) + RotatedAxis()

VlnPlot(killi_youngandaged_integrated, features= c("SOX2","PCNA","HES3"), cols=c("green","purple"),split.by ="sample") + RotatedAxis()
VlnPlot(killi_youngandaged_integrated, features= c("UCP2"), cols=c("green","purple"),split.by ="sample") + RotatedAxis()

FeaturePlot(killi_youngandaged_integrated,reduction="tsne",cols=c("grey","blue"),features= "UCP2",split.by ="sample",label = T) + RotatedAxis()
gc()
Heatmap
FeaturePlot(killi_youngandaged_integrated,reduction="tsne",cols=c("yellow","blue"),features= c("CD68","AIF1","C1QA","LCP1"),split.by ="sample",label = T) + theme_test() + RotatedAxis()


install.packages("enrichR")
BiocManager::install("goseq")
library(goseq)
library(enrichR)
listEnrichrDbs()
??enrichR
top30_alias <- read.csv(file="C:/Users/u0129074/Documents/ScRNAseq_paper/Redo_november2020/march2021/onlyyoung_DEgenes.csv")
head(rg3_changes_new)
websiteLive <- TRUE

suppressPackageStartupMessages(library(escape))
suppressPackageStartupMessages(library(dittoSeq))
suppressPackageStartupMessages(library(SingleCellExperiment))
#GSEA
GS.celltypes <- getGeneSets(library = "C8")
GS.hallmarks <- getGeneSets(library = "H")

YandA.seurat <- enrichIt(obj = killi_youngandaged_integrated, gene.sets = GS.celltypes, groups = 500, cores = 8)
YandA.seurat_h <- enrichIt(obj = killi_youngandaged_integrated, gene.sets = GS.hallmarks, groups = 1000, cores = 4)

YandA.seurat <- enrichIt(obj = killi_youngandaged_integrated, gene.sets = GS.celltypes, groups = 500, cores = 8)


YandA.seurat_h
killi_youngandaged_integrated <- Seurat::AddMetaData(killi_youngandaged_integrated, YandA.seurat_h)
newcolors <- colorRampPalette(c("#0348A6", "#7AC5FF", "#C6FDEC", "#FFB433", "#FF4B20"))
newcolors
dittoHeatmap(killi_youngandaged_integrated, genes = NULL, metas = names(YandA.seurat), 
             annot.by = "sample", 
             fontsize = 7, 
             cluster_cols = TRUE,
             heatmap.colors = newcolors(50))





dittoHeatmap(killi_youngandaged_integrated,
                          metas = c("HALLMARK_APOPTOSIS", "HALLMARK_DNA_REPAIR", "HALLMARK_P53_PATHWAY"), 
                          annot.by = "sample", 
                          fontsize = 7,
                          heatmap.colors = newcolors(10))

multi_dittoPlot(killi_youngandaged_integrated, vars = c("HALLMARK_APOPTOSIS", "HALLMARK_DNA_REPAIR", "HALLMARK_P53_PATHWAY"), 
                group.by = "sample", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Enrichment Scores", 
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018")
if (websiteLive) {
  enriched_Y_aged <- enrichr(c(rg3_changes_new$gene), dbs)
}

if (websiteLive) enriched_Y_aged[["GO_Molecular_Function_2018"]]

if (websiteLive) plotEnrich(enriched_Y_aged[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")



complexity.per.cell <- apply(killi_youngandaged_integrated@assays,
                             2, function(x) sum(x>0))
VlnPlot(complexity.per.cell)

FeatureScatter(Neuronal_clusters,feature1 = "GAP43","NEUROD2")
batch.effect <- technical_batches[row.names(killi_youngandaged_integrated@meta.data)]


# Get batches based on cell names
technical_batches <- sapply(colnames(killi_youngandaged_integrated@meta.data),
                            FUN=function(x){strsplit(x,"_",fixed=TRUE)[[1]][1]})

# Turn to numbers and add cell names to them
technical_batches <- as.numeric(as.factor(technical_batches))
names(technical_batches) <- colnames(bipolar.seurat.raw@raw.data)

?DotPlot
DefaultAssay(killi_youngandaged_integrated) <- "RNA"
RG1.markers <- FindConservedMarkers(killi_youngandaged_integrated, ident.1 = "RG1", grouping.var = "sample", verbose = FALSE)
RG2.markers <- FindConservedMarkers(killi_youngandaged_integrated, ident.1 = "RG2", grouping.var = "sample", verbose = FALSE)
RG3.markers <- FindConservedMarkers(killi_youngandaged_integrated, ident.1 = "NE-RG3", grouping.var = "sample", verbose = FALSE)
RG4.markers <- FindConservedMarkers(killi_youngandaged_integrated, ident.1 = "RG4", grouping.var = "sample", verbose = FALSE)
NGP.markers <- FindConservedMarkers(killi_youngandaged_integrated, ident.1 = "NGP", grouping.var = "sample", verbose = FALSE)


ImmC.markers <- FindConservedMarkers(killi_youngandaged_integrated, ident.1 = "ImmC", grouping.var = "sample", verbose = FALSE)
ImN1.markers <-FindConservedMarkers(killi_youngandaged_integrated, ident.1 = "ImN1", grouping.var = "sample", verbose = FALSE)
ImN2.markers <-FindConservedMarkers(killi_youngandaged_integrated, ident.1 = "ImN2", grouping.var = "sample", verbose = FALSE)
mN1.markers <-FindConservedMarkers(killi_youngandaged_integrated, ident.1 = "mN1", grouping.var = "sample", verbose = FALSE)
mN2.markers <-FindConservedMarkers(killi_youngandaged_integrated, ident.1 = "mN2", grouping.var = "sample", verbose = FALSE)

write.table(NGP.markers, file="C:/Users/u0129074/Documents/ScRNAseq_paper/Redo_november2020/Aug2021/forGOanalysis/NGP_consmarkers.txt" ,sep="\t", col.names=TRUE)
write.table(RG3.markers, file="C:/Users/u0129074/Documents/ScRNAseq_paper/Redo_november2020/sep2021/forGOanalysis/NE-RG3_consmarkers.txt" ,sep="\t", col.names=TRUE)

rg3_changes_new <- read.csv(file="C:/Users/u0129074/Documents/ScRNAseq_paper/Redo_november2020/march2021/GOanalysis/RG3_onlygenesnames.csv")
View(RG3.markers)
head(rg3_changes_new)
write.table(mN1.markers, file="C:/Users/u0129074/Documents/rstudio_jan2020/20jan_afterintegration/Integration_withY&A/Mincells_2_new/matureN1diff_markers_new.txt", sep="\t", col.names=TRUE)

#differential expression
killi_youngandaged_integrated$celltype.aged <- paste(Idents(killi_youngandaged_integrated), killi_youngandaged_integrated$AgedK, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.stim"
b.interferon.response <- FindMarkers(immune.combined, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
head(b.interferon.response, n = 15)

RG1_diffexpr <- FindMarkers(killi_youngandaged_integrated, ident.1 = "Astro-RG1", verbose = TRUE)
RG1_diffexpr
RG3_diffexpr <- FindMarkers(killi_youngandaged_integrated, ident.1 = "NE-RG3_AgedK", ident.2 = "NE-RG3_YoungK", verbose = FALSE)
tail(killi_youngandaged_integrated)
head(b.interferon.response, n = 15)
killi_youngandaged_integrated$sample
FeaturePlot(killi_youngandaged_integrated, features = c("HEPACAM","WNT8B","PCNA"), split.by = "sample", max.cutoff = 3,
            cols = c("grey", "red"))

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
RG1.cells <- subset(killi_youngandaged_integrated, idents = "RG1")
Idents(RG1.cells) <- "sample"
avg.RG1.cells <- as.data.frame(log1p(AverageExpression(RG1.cells, verbose = FALSE)$RNA))
avg.RG1.cells$gene <- rownames(avg.RG1.cells)
write.table(avg.RG1.cells, file="C:/Users/u0129074/Documents/ScRNAseq_paper/Redo_november2020/Aug2021/forGOanalysis/RG1_Diff.txt",sep="\t")
avg.RG1.cells
RG2.cells <- subset(killi_youngandaged_integrated, idents = "RG2")
Idents(RG2.cells) <- "sample"
avg.RG2.cells <- as.data.frame(log1p(AverageExpression(RG2.cells, verbose = FALSE)$RNA))
avg.RG2.cells$gene <- rownames(avg.RG2.cells)
write.table(avg.RG2.cells, file="C:/Users/u0129074/Documents/ScRNAseq_paper/Redo_november2020/Aug2021/forGOanalysis/RG2_Diff.txt",sep="\t")

RG3.cells <- subset(killi_youngandaged_integrated, idents = "RG3")
Idents(RG3.cells) <- "sample"
avg.RG3.cells <- as.data.frame(log1p(AverageExpression(RG3.cells, verbose = FALSE)$RNA))
avg.RG3.cells$gene <- rownames(avg.RG3.cells)
write.table(avg.RG3.cells, file="C:/Users/u0129074/Documents/ScRNAseq_paper/Redo_november2020/Aug2021/forGOanalysis/RG3_Diff.txt",sep="\t")
avg.RG3.cells

RG4.cells <- subset(killi_youngandaged_integrated, idents = "RG4")
Idents(RG4.cells) <- "sample"
avg.RG4.cells <- as.data.frame(log1p(AverageExpression(RG4.cells, verbose = FALSE)$RNA))
avg.RG4.cells$gene <- rownames(avg.RG4.cells)
write.table(avg.RG4.cells, file="C:/Users/u0129074/Documents/ScRNAseq_paper/Redo_november2020/Aug2021/forGOanalysis/RG4_Diff.txt",sep="\t")
log1p(10)
NGP.cells <- subset(killi_youngandaged_integrated, idents = "NGP")
Idents(NGP.cells) <- "sample"
avg.NGP.cells <- as.data.frame(log1p(AverageExpression(NGP.cells, verbose = FALSE)$RNA))
avg.NGP.cells$gene <- rownames(avg.NGP.cells)
write.table(avg.NGP.cells, file="C:/Users/u0129074/Documents/ScRNAseq_paper/Redo_november2020/Aug2021/forGOanalysis/NGP_Diff.txt",sep="\t")


NBN.cells <- subset(killi_youngandaged_integrated, idents = "ImN2")
Idents(NBN.cells) <- "sample"
avg.NBN.cells <- as.data.frame(log1p(AverageExpression(NBN.cells, verbose = FALSE)$RNA))
avg.NBN.cells$gene <- rownames(avg.NBN.cells)

#Extracting number of cells per cluster
library(data.table)
library(magrittr)

#All cell clusters
md_all <- killi_youngandaged_integrated@meta.data %>% as.data.table


#For proliferatiom clusters
md_pc <- proliferation_clusters@meta.data %>% as.data.table
## count the number of cells per unique combinations of "Sample" and "seurat_clusters"
md_pc
md_pc[, .N, by = c("sample", "CellType")]
pie(table(Idents(pr)))+ scale_color_aaas()

#For neuronal clusters
neuronal_clusters@meta.data
md_n <- neuronal_clusters@meta.data %>% as.data.table
## count the number of cells per unique combinations of "Sample" and "seurat_clusters"
md_n
md_n[, .N, by = c("sample", "CellType")]
pie(table(Idents(proliferation_clusters)))+ scale_color_aaas()

cxfvtbgyhju#col_new= c("#FA340D","lightskyblue","darkcyan","#38577A","#FAB693","lightsteelblue","powderblue","#D80202")
slices=table(Idents(proliferation_clusters))
slices
pie(table(Idents(proliferation_clusters)), border="grey")
pct <- round(slices/sum(slices)*100)
lbls= c("NGP","RG1","RG2","RG3","RG4","Inter-NGP","Inter-RG")

lbls <- paste(lbls, pct)
lbls <- paste(lbls,"%",sep="")
pie(slices)

proliferation_clusters
genes.to.label = c("CX43", "FABP7A","SLC1A2", "FAM167AB", "GLUL (2 OF 2)", "SOX2", "HEPACAM","WNT8B", "VIM", "HES1","HES5 (3 OF 9)","HES5(2OF2)","SOX3","EPD")
genes_new = c("NEUROG1","NOTCH1A","NOTCH1B","DLA","DLD","ASCL1A","ASCL1B","ID4","PCNA","KIAA0101","HMGB2A","STMN1A","MKI67","HES5 (3 OF 9)","HES5(2OF2)","HES3","HES1","EOMESA","HMGA1","HMGN1","HMGN3","HMGB2A","NOTCH3")
genes_to_labl = c("GAP43","DCX","STMN1B","NEUROD2")
p1 <- ggplot(avg.RG1.cells, aes(YoungK, AgedK)) + geom_point(colour="chartreuse2") + ggtitle("RG1 differences")
p1 <- LabelPoints(plot = p1, points = genes.to.label.RG1, repel = TRUE)
p2 <- ggplot(avg.RG2.cells, aes(YoungK, AgedK)) + geom_point(colour="pink") + ggtitle("RG2 differences")
p2 <- LabelPoints(plot = p2, points = genes.to.label.RG2, repel = TRUE)
p3 <- ggplot(avg.RG3.cells, aes(YoungK, AgedK)) + geom_point(colour="violet") + ggtitle("RG3 differences")
p3 <- LabelPoints(plot = p3, points = genes.to.label.RG3, repel = TRUE)
p4 <- ggplot(avg.RG4.cells, aes(YoungK, AgedK)) + geom_point(colour="cornflowerblue") + ggtitle("RG4 differences")
p4 <- LabelPoints(plot = p4, points = genes.to.label.RG4, repel = TRUE)

p1+p2+p3+p4
p5 <- ggplot(avg.NGP.cells, aes(YoungK, AgedK)) + geom_point(colour="royalblue1") + ggtitle("NGP differences")
p5 <- LabelPoints(plot = p5, points = genes_new, repel = TRUE)
p5

p1 <- ggplot(avg.RG1.cells, aes(YoungK, AgedK)) + geom_point(colour="chartreuse2") + ggtitle("RG1 differences")
p1 <- LabelPoints(plot = p1, points = sen_core_genes, repel = TRUE)
p2 <- ggplot(avg.RG2.cells, aes(YoungK, AgedK)) + geom_point(colour="pink") + ggtitle("RG2 differences")
p2 <- LabelPoints(plot = p2, points =sen_core_genes , repel = TRUE)
p3 <- ggplot(avg.RG3.cells, aes(YoungK, AgedK)) + geom_point(colour="violet") + ggtitle("RG3 differences")
p3 <- LabelPoints(plot = p3, points = sen_core_genes, repel = TRUE)
p4 <- ggplot(avg.RG4.cells, aes(YoungK, AgedK)) + geom_point(colour="cornflowerblue") + ggtitle("RG4 differences")
p4 <- LabelPoints(plot = p4, points = sen_core_genes, repel = TRUE)
p4
p5 <- ggplot(avg.NGP.cells, aes(YoungK, AgedK)) + geom_point(colour="royalblue1") + ggtitle("NGP differences")
p5 <- LabelPoints(plot = p5, points = sen_core_genes, repel = TRUE)
p1
p3

DotPlot(proliferation_clusters, features =sen_core_genes, cols = c("green", "purple"), split.by ="sample",cluster.idents = F) + RotatedAxis()
DotPlot(proliferation_clusters, features =sen_effector_genes, cols = c("green", "purple"), split.by ="sample",cluster.idents = F) + RotatedAxis()

DefaultAssay(killi_youngandaged_integrated) <- "RNA"
Jolien_saspfactors = c("TP53","MDM2","CDKN1BA","TNFA","CCL20","CXCR3","CXCL8","IFNB2","FGF7","MMP9","LMNB1","FAS","OLA.22466")
Jolien_saspfactors = c("TP53","MDM2","CDKN1BA","CDKN1B","CXCL8","MAP2K6","RRAS2","BRF1B","AREGB","IGFBP3 (2 OF 2)","IGFBP4","IGFBP5 (2 OF 2)","PRODH (2 OF 3)","SMURF2","PGF","CTNNB1","ANGPTL2","LMNB1","FAS")
Jolien_saspfactors_new = c("MDM2","CDKN1BA","CDKN1B","CXCL8","MAP2K6","RRAS2","BRF1B","AREGB","IGFBP3 (2 OF 2)","IGFBP4","IGFBP5 (2 OF 2)","PRODH (2 OF 3)","SMURF2","CTNNB1","ANGPTL2","LMNA (2 OF 2)","LMNB1","FAS")

data("killi_youngandaged_integrated")
head(AverageExpression(object = Jolien_celltypes))

GetAssayData(object = killi_youngandaged_integrated, slot = "counts")
killi_youngandaged_integrated <- SetAssayData(object = killi_youngandaged_integrated, slot = "scale.data")


Jolien_celltypes <- subset(killi_youngandaged_integrated,  idents = c("NGP","RG","EPD"))
Jolien_celltypes$nGene
proliferation_clusters
DotPlot(Jolien_celltypes, features =Jolien_saspfactors_new) + RotatedAxis()
VlnPlot(Jolien_celltypes$sample, features =Jolien_saspfactors_new, cols = c("blue", "grey", "blue"), stack=T)

?DotPlot
saveRDS(proliferation_clusters, "~/ScRNAseq_paper/Redo_november2020/june2021/Y&afish_PCs_final2.Rds")

killi_youngandaged_integrated@active.ident

genes.to.label.NBN= c("GAP43","ELAVL3","OLA.22466","STMN1B","NEUROD2","MARCKS")
p6 <- ggplot(avg.NBN.cells, aes(YoungK, AgedK)) + geom_point(colour="cornflowerblue") + ggtitle("NBN differences")
p6 <- LabelPoints(plot = p6, points = genes.to.label.NBN, repel = TRUE)

#senescence-related genes expression

sen_core_genes = c("CDKN2A", "BMI1", "CDKN2B","PDRG1","TP53","LMNB1", "HMGA1", "CHEK1", "CHEK2", "PRODH (2 OF 3)", "TNFRSF10B", "CDKN1A","DAO")

sen_effector_genes = c("PPP1CAA", "AHCY", "BRF1B", "MAP2K3", "MAP2K6", "SMURF2", "TGFB1I1", "SRSF1 (1 OF 2)", "ANGPTL2")

saspfactors = c("CTNNB1","CXCL12","CXCL14","IGFBP2","IGFBP3 (2 OF 2)","IGFBP4","IGFBP5B","IGFBP6B","IL15L","IL1B","MIF","MMP11A","MMP13","MMP14","PGF","PLAT","TIMP2A","NFU-G-1-023773","KITLGA","SERPINE2","TNFRSF1A","HGF","NRG1","EREG","AREGB")
saspfactors_CZ =c("LAMB1","PPP3CCB", "SLC25A6", "RRAS2", "HRAS", "MAPK1", "MAPK3", "MAP2K2", "SIRT2", "VDAC2")



all_sen_genes =c("CDKN2A", "BMI1", "CDKN2B","PDRG1","LMNB1", "HMGA1", "CHEK1", "CHEK2", "PRODH (2 OF 3)", "TNFRSF10B", "CDKN1A","DAO","PPP1CAA", "AHCY", "BRF1B", "MAP2K3", "MAP2K6", "SMURF2", "TGFB1I1", "SRSF1 (1 OF 2)", "ANGPTL2")
genes.to.label.RG1 = c("CX43", "FABP7A","SLC1A2", "FAM167AB","GRM7", "GLUL (2 OF 2)","SLC1A2B","SCARA3.3","S100A10A")
p1 <- ggplot(avg.RG1.cells, aes(YoungK, AgedK)) + geom_point(colour="chartreuse2") + ggtitle("RG1 differences")
p1 <- LabelPoints(plot = p1, points = genes.to.label.RG1, repel = TRUE)
genes.to.label.RG2 = c("CX43", "FABP7A","SLC1A2","S100A10A","SCARA3.3")
genes.to.label.RG3 = c("HES1", "FABP7A","VIM","SOX2","PCNA","HEPACAM","WNT8B")
genes.to.label.RG4 = c("EPD","AHSG","APOA1","APOA2","CLU","FABP7A","SCARA3.3","S100A10A")
p1
p2
p3
p4

p5 <- ggplot(avg.NBN.cells, aes(YoungK, AgedK)) + geom_point(colour="royalblue1") + ggtitle("NBN differences")
p5 <- LabelPoints(plot = p5, points = sen_core_genes, repel = TRUE)
p5
??VolcanoPlot
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
browseVignettes("EnhancedVolcano")
avg.NBN.cells$YoungK
View(avg.NBN.cells)
EnhancedVolcano(avg.NBN.cells,
                lab = rownames(avg.NBN.cells),
                x = 'YoungK',
                y = Y,pointSize = 3.0,
                labSize = 6.0)

p <- ggplot(data=, aes(x=avg_log2FC, y=-log10(p_val))) + geom_point() + theme_minimal()
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p2

de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]

ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()



FeaturePlot(killi_youngandaged_integrated, reduction="tsne", features=c("GAP43","MARCKS"),cols=c("green", "magenta"),blend=T, split.by = "sample", label=T)
?FeaturePlot

RG1.cells
FeatureScatter(object = killi_youngandaged_integrated, cells="RG1",feature1 = 'FABP7A', feature2 = 'CX43')
VlnPlot(killi_youngandaged_integrated, stack=T, split.by ="sample", combine=T, pt.size = 1.2, features =c("WNT8B","PCNA"))
VlnPlot(neuronal_clusters, stack=T, split.by ="sample", combine=T, pt.size = 0, features =c("GAP43","ELAVL3"))

neuronal_clusters <- subset(killi_youngandaged_integrated, idents = c("mN1","mN2","mN1","ImN2","mN3","mN4","ImN1"))
VlnPlot(neuronal_clusters, features=c("SYP","MAP2","SV2A","ADCYAP1B","ADCYAP1R1A","SLC17A6A","SLC17A6B","SLC17A7","SLC6A1 (1 OF 2)","GABBR1","GABRB2","GRIA2","GRIN1A","LOC107384443","GAD2","PVALB"),cols=barplotcol_new,log=F, sort=F,fill.by="ident",flip=T, stack=T,combine=T)


RGonly_clusters <- subset(killi_youngandaged_integrated, idents = c("RG1","RG2","RG3","RG4"))
DotPlot(RGonly_clusters, features = c("SLC1A2","FABP7A","GLUL (2 OF 2)","CX43","S100B","SLC1A2","SOX2","SOX3","GFAP","HEPACAM","WNT8B","HES4.1","APOA1","NDRG4","CLU"),cols = c("green", "purple"), split.by ="sample",cluster.idents = F) + theme_clean() + RotatedAxis()

#Senescence in cells
AUC
cells_rankings_a <- AUCell_buildRankings(aged_isoseq_matrix_new, nCores=1, plotStats=TRUE)
dim(aged_isoseq_matrix_new)
cells_rankings_a
aging_geneset <- GeneSet(all_sen_genes, setName="aging_geneset")
aging_cells_AUC <- AUCell_calcAUC(all_sen_genes, cells_rankings_a)
save(aging_cells_AUC, file="aged_cells_AUC.RData")

set.seed(123)
agedcells_assignment <- AUCell_exploreThresholds(aging_cells_AUC, plotHist=TRUE, assign=TRUE) 
dev.off()

#AUC in young killifish brain
cells_rankings_y <- AUCell_buildRankings(young_isoseq_matrix_new, nCores=1, plotStats=TRUE)
dim(young_isoseq_matrix_new)
cells_rankings_y
yg_geneset <- GeneSet(sen_core_genes, setName="young_geneset")
yg_cells_AUC <- AUCell_calcAUC(all_sen_genes, cells_rankings_y)
save(yg_cells_AUC, file="young_cells_AUC.RData")

set.seed(123)
youngcells_assignment <- AUCell_exploreThresholds(yg_cells_AUC, plotHist=TRUE, assign=TRUE) 
dev.off()

library(plotly)
cells= c("NE-RG3",	"Inter-NGP", "EPD-RG4",	"Inter-RG",	"Astro-RG1","NGP","Astro-RG2")
No_Young = c(18.51851852,
             16.52421652,
             13.67521368,
             26.92307692,
             9.829059829,
             9.544159544,
             4.985754986)

No_Aged = c(8.282208589,
            38.46625767,
            11.65644172,
            18.65030675,
            8.343558282,
            6.319018405,
            8.282208589)
fig <- plot_ly(data, x =~No_Young , y = ~cells, type = 'bar', name = 'Young',orientation = 'h')
fig <- fig %>% add_trace(x = ~No_Aged, name = 'Aged')
fig <- fig %>% layout(xaxis = list(title = 'Cell proportion (%)', tickfont = list(size = 15)), yaxis=list(tickfont = list(size = 15)), barmode = 'group')


fig
?plot_ly

p1 <- DotPlot(proliferation_clusters,features = c("TMSB4X","ASCL1","SOX4","HMGB2B","MARCKS","NEUROG1","NOTCH1B","DLA","DLD","FGFR1A","FGFR1 (2 OF 2)","FGFR2","FGFR3","FGFR4","ID1","ID2","ID3","ID4","LOC107388996","BMP1","HES1","HES2","HES4","HES5 (3 OF 9)","HES5 (2 OF 2)","HES6","HMGA1","HMGN1","HMGN3","HMGB2A","NOTCH3"),cols=c("steel blue","darkorange1"), split.by="sample") + RotatedAxis() + coord_flip()
p2 <- DotPlot(proliferation_clusters, features = c("NEUROG1","NOTCH1B","DLA","DLD","FGFR1A","FGFR1 (2 OF 2)","FGFR2","FGFR3","FGFR4"), cols=c("steel blue","darkorange1"), split.by="sample") + RotatedAxis() + coord_flip()
p4 <- DotPlot(proliferation_clusters, features=c("ID1","ID2","ID3","ID4","ASCL1","LOC107388996","BMP1"),cols=c("steel blue","darkorange1"), split.by="sample") + RotatedAxis() + coord_flip()
p3 <- DotPlot(proliferation_clusters, features = c("HES1","HES2","HES4","HES5 (3 OF 9)","HES5 (2 OF 2)","HES6","HMGA1","HMGN1","HMGN3","HMGB2A","NOTCH3"), cols=c("steel blue","darkorange1"), split.by="sample") + RotatedAxis() + coord_flip()
p1+p2+p3+p4

p3


DotPlot(proliferation_clusters, features = c("HES4"), cols=c("steel blue","darkorange1"), split.by="sample") + RotatedAxis() + coord_flip()

DotPlot(killi_youngandaged_integrated, features = c("SLC7A11"), cols=c("steel blue","darkorange1"), split.by="sample") + RotatedAxis() + coord_flip()

DotPlot(proliferation_clusters, features = c("SLC7A11","SLC6A11"), cols=c("steel blue","darkorange1"), split.by="sample") + RotatedAxis() + coord_flip()
DotPlot(neuronal_clusters, features = c("SLC7A11","SLC6A11"), cols=c("steel blue","darkorange1"), split.by="sample") + RotatedAxis() + coord_flip()



proliferation_clusters$CellType
