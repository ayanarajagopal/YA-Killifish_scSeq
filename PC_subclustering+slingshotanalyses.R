require(scales)
library(scCustomize)
library("Nebulosa")
library(scCustomize)
library(ComplexHeatmap)
library(dittoSeq)
library("DropletUtils")
library(dplyr)
library(WriteXLS)
BiocManager::install(c("ComplexHeatmap", "dittoSeq", "DropletUtils", "Nebulosa"))
#Optional CRAN packages
install.packages(c("ggpubr", "hdf5r", "rliger"))

##NGP + RG all types sub-clustering
PCs_integrated <- subset(killi_young_integrated, idents = c("Astro-RG1","Astro-RG2","NE-RG3","EPD-RG4","NGP","Intercell.PC"))
PCs_integrated
DefaultAssay(PCs_integrated) <- "integrated"
View(rownames(PCs_integrated))
#PCs_integrated <- SCTransform(PCs_integrated, verbose = TRUE )
PCs_integrated <- RunPCA(PCs_integrated, verbose=FALSE)
ElbowPlot(PCs_integrated, ndims = 40)
PCs_integrated<- FindNeighbors(object = PCs_integrated, 
                               dims = 1:20)

dim(PCs_integrated)
PCs_integrated$seurat_clusters
# Determine the clusters for various resolutions                                
PCs_integrated <- FindClusters(object = PCs_integrated,
                               resolution = 0.5)

#Calc. of TSNE/UMAP
PCs_integrated <- RunTSNE(object = PCs_integrated, reduction= "pca", dims=1:20)
PCs_integrated <- RunUMAP(object = PCs_integrated, reduction= "pca", dims = 1:20)
#spread=1, min.dist = 0.4
DefaultAssay <- "RNA"

col_new_ridgeplot_allPCs = c("lightsteelblue","#FAB693", "gray","#F4BB44", "#FFFF8F","#FDDA0D","#E4D00A", "#FA340D", "powderblue","darkcyan")

p1<- DimPlot(PCs_integrated, reduction = "tsne", pt.size = 1.8,label=T, split.by="sample")
p2 <- DimPlot(PCs_integrated, reduction = "tsne", group.by="CellType",pt.size = 1.8,label=T, order=T, cols=col_new_ridgeplot_allPCs) + theme_test()
p1+p2

table(Idents(PCs_integrated))

DefaultAssay(PCs_integrated) <- "RNA"
PCs_integrated <- NormalizeData(PCs_integrated, verbose = FALSE)
DefaultAssay(PCs_integrated)
?FindAllMarkers
PCs_integrated_markers <- FindAllMarkers(PCs_integrated, min.pct = 0.25, only.pos = F, logfc.threshold = 0.25)

head(PCs_integrated_clustersmarkers)
youngPCs_PCs_cluster_top200 <- PCs_integrated_cons_markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=200, wt = avg_log2FC)
youngPCs_allS_cluster_top
youngPCs_top100genes <- PCs_integrated_clustersmarkers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=100, wt = avg_log2FC)
youngPCs_top100genes
youngPCs_top10genes <- PCs_integrated_clustersmarkers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=10, wt = avg_log2FC)
youngPCs_top10genes

#####Find all markers #####
write_xlsx(youngPCs_PCs_cluster_top200, col_names = T, path="C:/Users/u0129074/Documents/rstudio_jan2020/20jan_afterintegration/Young_old_separate/Young_Isoseq_latest/July2020/young_aggr_oct2022/Analysis/youngPCs_clubbed_top200genes_res1.7_fin.xlsx")
write_xlsx(youngPCs_top100genes_renamed, col_names = T, path="C:/Users/u0129074/Documents/rstudio_jan2020/20jan_afterintegration/Young_old_separate/Young_Isoseq_latest/July2020/young_aggr_oct2022/Analysis/youngPCs_top100genes_res0.5_new.xlsx")

#Cell type annotation
names(youngPCs_all_0.5_new) <- levels(PCs_integrated)
PCs_integrated <- RenameIdents(PCs_integrated, youngPCs_all_0.5_new)

#Cell type identification

####Annnotating the cell types#####
youngPCs_all_0.5_new <- c("NGP.1","Intercell.1","NGP.2","Intercell.2","Intercell.3","Astro-RG1","NE-RG3","Intercell.4","Astro-RG2","EPD-RG4")

names(youngPCs_all_0.5_new) <- levels(PCs_integrated)
PCs_integrated <- RenameIdents(PCs_integrated, youngPCs_all_0.5_new)
PCs_integrated$CellType <- Idents(PCs_integrated)
PCs_integrated$CellType
PCs_integrated$seurat_clusters <- Idents(PCs_integrated)
PCs_integrated$seurat_clusters

table(Idents(PCs_integrated))

# Reorder identity classes
levels(x = PCs_integrated)
#> [1] "B" "A" "C"
levels(x = PCs_integrated) <- c("NE-RG3","NGP.1","NGP.2","Intercell.1","Intercell.2","Intercell.3","Intercell.4","Astro-RG1","Astro-RG2","EPD-RG4")
levels(x = PCs_integrated)
DotPlot(PCs_integrated, features= c("ZIC2","LOC107375911","SOX2","SOX3","SOX4","HES5 (3 OF 9)","HES5 (2 OF 2)","PCNA","HMGB2A","NEUROG1","STMN1A","ASCL1B","GLUL (2 OF 2)","SLC1A2B","CX43","NDRG4","FABP7A","EPD","APOA1"), cols=c("green","purple")) + RotatedAxis()+
  theme(text = element_text(face = "plain"),
        axis.text.x=element_text(angle=45, hjust=1, size=12),
        axis.title = element_text(size=12),
        axis.title.y.right = element_text(size = 12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        axis.line = element_line(size=1))


saveRDS(PCs_integrated, "C:/Users/u0129074/Documents/ScRNAseq_paper/youngkillifish_PCs_mar2023.rds") 

##ccphases analysis
PCs_integrated <- CellCycleScoring(PCs_integrated, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
?CellCycleScoring
PCs_integrated[[]]
as_tibble(PCs_integrated[[]]) %>%
  ggplot(aes(Phase)) + geom_bar()
PCs_integrated@meta.data %>%
  group_by(CellType,Phase) %>%
  count() %>%
  group_by(CellType) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=CellType,y=percent, split,fill=Phase), cols=c("blue","#ffa10a","#00b300")) +
  geom_col() + RotatedAxis()+
  ggtitle("Percentage of cell cycle phases per cluster")

####Lineage inference analysis with slingshot####
library(slingshot)
library(reticulate)
library(scVelo)
library(ggbeeswarm)
library(ggthemes)
library(viridis)
library(SingleCellExperiment)
library(tradeSeq)
library(dplyr)
library(Matrix)
library(scCustomize)
library(qs)
install.packages("scCustomize")
BiocParallel::register(BiocParallel::SerialParam())

PCs_lineages_Yall <- subset(x=PCs_integrated, idents= c("NGP.1","Intercell.1","NGP.2","Intercell.2","Intercell.3","Astro-RG1","NE-RG3","Intercell.4","Astro-RG2","EPD-RG4"))
ElbowPlot(PCs_lineages_Yall, ndims=30)
PCs_lineages_Yall <- RunTSNE(PCs_lineages_Yall,dims = 1:30)
PCs_lineages_Yall

counts_data_Y_allS= GetAssayData(object=PCs_lineages_Yall[["RNA"]], slot="counts")
counts_data_Y_allS

counts(PCs_lineages_Yall)
assays(PCs_lineages_Yall)
head(scale_data_Y_allS)

FeaturePlot(PCs_lineages_Yall, features=c("SOX2","PCNA","MKI67","ELAVL3","S100B","ASCL1B"),label=T )
col_new_ridgeplot_allPCs = c("lightsteelblue","#FAB693", "gray","#F4BB44", "#FFFF8F","#FDDA0D","#E4D00A", "#FA340D", "powderblue","darkcyan")
DimPlot(PCs_lineages_Yall, reduction="tsne", label=T, cols=col_new_ridgeplot_allPCs, pt.size=1.2)

####Save the objects as separate matrices for slingshot input###

dimred_umap <- PCs_lineages_Yall@reductions$umap@cell.embeddings
dimred_tsne <- PCs_lineages_Yall@reductions$tsne@cell.embeddings

clustering_PC <- PCs_lineages_Yall$seurat_clusters
#counts_PCs <- PCs_lineages_Yall@assays$RNA@counts[PCs_lineages_Yall@assays$RNA@var.features, ]

PCs_lineages_Yall@assays$integrated@var.features

PCs_allS_lineages <- getLineages(data = dimred_tsne,
                                 clusterLabels = clustering_PC,
                                 #end.clus = c("11","7","10","9","5"), #define how many branches/lineages to consider
                                 start.clus = "NE-RG3") #define where to start the trajectories

head(pal[clustering_PC], n=10)

curves_PC <- getCurves(PCs_allS_lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves_PC

#Show all lineages
slingLineages(curves_PC)

par(mfrow=c(1,1))
plot(dimred_tsne[,1:2], col = pal[clustering_PC],  cex=.5,pch = 16)
for(i in levels(clustering_PC)){ 
  text( mean(dimred_tsne[clustering_PC==i,1]),
        mean(dimred_tsne[clustering_PC==i,2]), labels = i,font = 2) }
plot(dimred_tsne, col = pal[clustering_PC],  cex=0.5,pch = 16)
lines(SlingshotDataSet(curves_PC), lwd = 3, col = 'black')
dev.off()

####Differential expression analysis across lineages: nbGAM analysis#####
# Removed some genes to speed up the computations
filt_counts_PCs <- counts_data_Y_allS[rowSums(counts_data_Y_allS > 5) > ncol(counts_data_Y_allS)/100, ]
dim(filt_counts_PCs)

#Run Fitgam
sce_PC <- fitGAM(counts = as.matrix(filt_counts_PCs), pseudotime=pdtime_PC,  sds = curves_PC, verbose=T)
plotGeneCount(curves_PC, filt_counts_PCs, clusters = clustering_PC, models = sce_PC)

head(filt_counts_PCs)
plotGeneCount(curves_PC, counts=scale_data_Y_allS, gene="HES1")

library(dplyr)
plot_differential_expression <- function(feature_id) {
  feature_id <- pseudotime_association_PCs %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
  cowplot::plot_grid(plotGeneCount(curves_PC, filt_counts_PCs, gene = feature_id[1], clusters = clustering_PC, models = sce_PC) + ggplot2::theme(legend.position = "none"), 
                     plotSmoothers(curves_PC, as.matrix(counts_data_Y_allS), gene = feature_id[1]))
}
###Genes that change with pseudotime
pseudotime_association_PCs <- associationTest(sce_PC)
pseudotime_association_PCs$fdr <- p.adjust(pseudotime_association_PCs$pvalue, method = "fdr")
pseudotime_association_PCs <- pseudotime_association_PCs[order(pseudotime_association_PCs$pvalue), ]
pseudotime_association_PCs$feature_id <- rownames(pseudotime_association_PCs)
?plotGeneCount
feature_id <- pseudotime_association_PCs %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)

##diff expr analysis
startRes <- startVsEndTest(sce_PC, pseudotimeValues=c(0,30))
oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce_PC)[oStart[39]]
sigGeneStart
plotSmoothers(sce_PC, counts_data_Y_allS, gene =sigGeneStart )
plotGeneCount(curves_PC, counts_data_Y_allS, gene =sigGeneStart)

endRes <- diffEndTest(sce_PC)
head(endRes)
o <- order(endRes$waldStat, decreasing = TRUE)
head(o)
sigGene <- names(sce_PC)[o[5]]
plotSmoothers(sce_PC, counts_data_Y_allS, sigGene)
plotGeneCount(curves_PC, counts_data_Y_allS, gene = sigGene)

###pattern test###
patternRes <- patternTest(sce_PC)
oPat <- order(patternRes$waldStat, decreasing = TRUE)
head(rownames(patternRes)[oPat], n=30)
plotSmoothers(sce_PC, counts_data_Y_allS, gene = "HES5 (3 OF 9)")
plotGeneCount(curves_PC, counts_data_Y_allS, gene = "CX43")

plotSmoothers(sce_PC, counts_data_Y_allS, gene = "ID4")
plotGeneCount(curves_PC, counts_data_Y_allS, gene = "ID4")

earlyDERes <- earlyDETest(sce_PC, knots = c(1, 2),l2fc = log2(1.5))
oEarly <- order(earlyDERes$waldStat, decreasing = TRUE)
head(rownames(earlyDERes)[oEarly], n=50)
plotSmoothers(sce_PC, counts_data_Y_allS, gene = rownames(earlyDERes)[oEarly][3])
plotGeneCount(curves_PC, counts_data_Y_allS, gene = rownames(earlyDERes)[oEarly][3])
