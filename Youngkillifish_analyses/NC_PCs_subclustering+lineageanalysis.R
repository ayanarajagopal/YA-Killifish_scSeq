####Neuron clustering and annotation####
NNs_newanno_allS <- subset(killi_young_integrated, idents = c("Intercell.NC", "mN1", "ImN1", "mN2", "mN3", "ImN2", "Inter-cell.PC", "mN4", "mN5", "mN6", "mN7", "NGP", "mN8", "NE-RG3"))
NNs_newanno_allS
#NNs_newanno_allS <- SCTransform(NNs_newanno_allS,verbose=T)
na <-is.na(rownames(NNs_newanno_allS))
DefaultAssay(NNs_newanno_allS) <- "integrated"
NNs_newanno_allS <- RunPCA(NNs_newanno_allS, verbose=FALSE)
ElbowPlot(NNs_newanno_allS, ndims=40)

NNs_newanno_allS <- FindNeighbors(object = NNs_newanno_allS, 
                                  dims = 1:30)

# Determine the clusters for various resolutions                                
NNs_newanno_allS <- FindClusters(object = NNs_newanno_allS,
                                 resolution = 0.9)
#Calc. of TSNE and perplexity optimization
NNs_newanno_allS <- RunTSNE(object = NNs_newanno_allS, reduction= "pca", dims=1:30)
NNs_newanno_allS <- RunUMAP(object = NNs_newanno_allS, reduction= "pca", dims = 1:30)

DimPlot(NNs_newanno_allS, reduction = "umap",label=T, pt.size = 1.5) + theme_test()
DimPlot(NNs_newanno_allS, reduction = "tsne",label=T, pt.size = 1.5) + theme_test()

DefaultAssay(NNs_newanno_allS) <- "RNA"
NNs_newanno_allS <- NormalizeData(NNs_newanno_allS, verbose = FALSE)

NNssubclusters_markers1_allS <- FindAllMarkers(NNs_newanno_allS, min.pct = 0.25, only.pos = F, logfc.threshold = 0.25)
head(NNssubclusters_markers1_allS)
NNs_newanno1_cluster_top200_allS <- NNssubclusters_markers1_allS %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=200, wt =avg_log2FC )
View(NNs_newanno1_cluster_top1_allS)
write.table(NNs_newanno1_cluster_top1_allS, file="C:/Users/u0129074/Documents/ScRNAseq_paper/Redo_november2020/june2021/NCs_newanno_subclusters_res0.8_Perplexity18_fin.txt", sep="\t", col.names=TRUE)
write_xlsx(NNs_newanno1_cluster_top200_allS, col_names = T, path="C:/Users/u0129074/Documents/rstudio_jan2020/20jan_afterintegration/Young_old_separate/Young_Isoseq_latest/July2020/young_aggr_oct2022/Analysis/youngNCs_top200genes_markers_res0.9_mar23.xlsx")

saveRDS(NNs_newanno_allS,"C:/Users/u0129074/Documents/ScRNAseq_paper/youngkillifish_NCPCs_mar2023.rds") 

pie(table(Idents(NNs_newanno_allS_onlyN)))

DotPlot(NNs_newanno_allS,idents=c("ImN1", "ImN2", "Ex-mN1", "ImN3", "Ex-mN2", "ImN4", "Ih-mN1", "Ex-mN3", "ImN5", "Ex-mN4", "Ex-mN5", "Ih-mN2", "Ih-mN3", "ImN6", "ImN7", "ImN8", "ImN9"), features= c("MEX3A","MIBP2","NEUROD2","ELAVL3","EOMESA","BHLHE22","SLC17A7","SLC17A6B","MAP2","SV2A","OLA.22466","STMN1B","DLX1","DLX2","DLX5","GAD2","GAD1B"), cols=c("green","purple"),cluster.idents=T)+ RotatedAxis()

NNs_newanno_allS_onlyN <- subset(NNs_newanno_allS,idents=c("ImN1", "ImN2", "Ex-mN1", "ImN3", "Ex-mN2", "ImN4", "Ih-mN1", "Ex-mN3", "ImN5", "Ex-mN4", "Ex-mN5", "Ih-mN2", "Ih-mN3", "ImN6", "ImN7", "ImN8", "ImN9") )
DimPlot(NNs_newanno_allS_onlyN, reduction="tsne",label=T, order=T,pt.size = 2) + theme_classic()
#Suppl table_NCs only
DefaultAssay(NNs_newanno_allS_onlyN)
NNs_newanno_allS_onlyN_markers <- FindAllMarkers(NNs_newanno_allS_onlyN, min.pct = 0.25, only.pos = F, logfc.threshold = 0.25)
onlyNeu_clusters_top200_allS <- NNs_newanno_allS_onlyN_markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=200, wt =avg_log2FC )
write.table(onlyNeu_clusters_top200_allS, file="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/Subclustering/Neus_only/Neusonly_17clusters_markers", sep="\t", col.names=TRUE)


#Suppl fig 3
top5_markers_NC <- Extract_Top_Markers(marker_dataframe = NNssubclusters_markers1_allS, num_genes = 5, named_vector = FALSE,
                                       make_unique = TRUE)
?Extract_Top_Markers
Clustered_DotPlot(seurat_object = NNs_newanno_allS, features = top5_markers_NC,colors_use_idents=my_color_palette_mNs)
NCs_PCs_allS <- subset(NNs_newanno_allS, idents=c("ImN1", "imN2", "mN1", "imN3", "Olf-mN2", "Intercell.1", "Ih-mN3", "Ih-imN4", "mN4", "mN5", "NGP.1", "Ih-mN6", "Ih-mN7", "PVALB+ mN6", "Intercell.2", "Ih-ImN2", "NB.1", "Olf-imN3", "NGP.2", "mN7", "NE-RG3"))

my_color_palette_mNs <- c("#D89000", "#C09B00", "#A3A500", "#39B600", "#00BF7D", "#00C1A3","#00BAE0", "#35A2FF")
VlnPlot(NNs_newanno_allS, features=c("SYP","MAP2","SV2A","PVALB","CCK (1 OF 2)","NPY","ADCYAP1B","ADCYAP1R1A","SLC17A6A","SLC17A6B","SLC17A7","SLC6A1 (1 OF 2)","SLC32A1 (1 OF 2)","GABBR1","GABRB2","GRIA2","GRIN1A","GAD1B","GAD2"),idents=c("Ex-mN1", "Ex-mN2", "Ex-mN3","Ih-mN1","Ex-mN4", "Ih-mN2","Ih-mN3","Ex-mN5","Ex-mN6"), log=F, sort=F, cols=my_color_palette_mNs, fill.by="ident",flip=T, stack=T,combine=T)
FeaturePlot(NNs_newanno_allS, features=c("EOMESA","DLX1"), order=T, reduction="tsne", label=T, pt.size = 1.)

##All NC-PC annotations##
NCs_res_0.9_mar <- c("ImN1", "ImN2", "Ex-mN1", "ImN3", "Ex-mN2", "ImN4", "Ih-mN1", "Ex-mN3", "ImN5", "NGP.1", "Ex-mN4", "Ex-mN5", "Ih-mN2", "Ih-mN3", "ImN6", "Intercell.PC", "ImN7", "ImN8", "Intercell.NC", "NGP.2", "ImN9", "NE-RG3")
NCs_res_0.9_dec23 <- c("Neu1", "Neu2", "Neu3", "Neu4", "Neu5", "Neu6", "Neu7", "Neu8","Neu9", "Neu10", "Neu11", "Neu12","Neu13", "Neu14", "Neu15", "Neu16","Neu17")


names(NCs_res_0.9_dec23) <- levels(NNs_newanno_allS_onlyN)
NNs_newanno_allS_onlyN <- RenameIdents(NNs_newanno_allS_onlyN, NCs_res_0.9_dec23)
pie(table(Idents(NNs_newanno_allS_onlyN)))
DotPlot(NNs_newanno_allS_onlyN, features= c("MEX3A","MIBP2","NEUROD2","ELAVL3","EOMESA","BHLHE22","OLA.22466","STMN1B","GAP43","SYT1","SV2A","SLC17A7","SLC17A6B","DLX1","DLX2","DLX5A","GAD2","GAD1B"), cols=c("green","purple"),cluster.idents=T)+ RotatedAxis() + theme_classic()
#####Lineage inference analyses#####
NCs_allS_lins <- subset(x=NNs_newanno_allS, idents=c("ImN1", "ImN2", "Ex-mN1", "ImN3", "Ex-mN2", "ImN4", "Ih-mN1", "Ex-mN3", "ImN5", "NGP.1", "Ex-mN4", "Ex-mN5", "Ih-mN2", "Ih-mN3", "ImN6", "Intercell.PC", "ImN7", "ImN8", "Intercell.NC", "NGP.2", "ImN9", "NE-RG3"))

ElbowPlot(NCs_allS_lins, ndims=50)
NCs_allS_lins <- RunTSNE(NCs_allS_lins,dims = 1:30)
dimred_NC_tsne_new <- NCs_allS_lins@reductions$tsne@cell.embeddings
clustering_NC_new <- NCs_allS_lins$seurat_clusters
clustering_NC_new
counts_NC_new <- as.matrix(NCs_allS_lins@assays$RNA@counts[NCs_allS_lins@assays$RNA@var.features, ])
NCs_allS_lineages_new <- getLineages(data = dimred_NC_tsne_new,
                                      clusterLabels = clustering_NC_new, start.clus="NE-RG3") #define how many branches/lineages to consider
curves_NCs_new <- getCurves(NCs_allS_lineages_new, thresh = 0.01, approx_points = 400, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
DimPlot(NCs_allS_lins, reduction = "tsne", pt.size=1.4, label=T) + NoLegend()
+ lines(SlingshotDataSet(curves_NCs_new), lwd = 3, col = 'black')

#All induced lineages
slingLineages(curves_NCs_new)

#plotting lineages
par(mfrow=c(1,2))
plot(dimred_NC_tsne_new[,1:2], col = pal[clustering_NC_new],  cex=0.5,pch = 16)
for(i in levels(clustering_NC_new)){ 
  text(mean(dimred_NC_tsne_new[pal[clustering_NC_new]==i,1]),
       mean(dimred_NC_tsne_new[pal[clustering_NC_new]==i,2]), labels = i,font = 1) }
plot(dimred_NC_tsne_new, col = pal[clustering_NC_new],cex=0.5, pch = 18)
lines_NC<- lines(SlingshotDataSet(curves_NCs_new), lwd = 3, col = 'black')
dev.off()
lines(SlingshotDataSet(curves_NCs_new), lwd = 3, col = 'black')
?lines

identities_new <- levels(Idents(NCs_allS_lins))

my_color_palette_new <- hue_pal()(length(identities_new))
my_color_palette_new
DimPlot(object = NCs_allS_lins, label = T, reduction="tsne", pt.size=1.2) + NoLegend()+
  scale_color_manual(values = my_color_palette_mod)+
  lines(SlingshotDataSet(curves_NCs_new), lwd = 3, col = 'black')
my_color_palette_mod <- c("#F8766D" ,"#EC8239", "#DB8E00", "#C79800", "#AEA200" ,"#8FAA00" ,"#64B200", "#00B81B", "#00BD5C" ,"#FAB693", "#00C1A7", "#00BFC4" ,"#00BADE", "#00B2F3", "#00A6FF", "#FFFF8F", "#B385FF", "#D874FD", "#FDDA0D", "gray" ,"#FF63B6", "lightsteelblue")

####Differential expression analysis across lineages; nbGAM analysis#####
# Removing some genes to speed up the computations
counts_data_NCs_new= GetAssayData(object=NCs_allS_lins[["RNA"]], slot="counts")
dim(counts_data_NCs_new)

filt_counts_NCs_new <- counts_data_NCs_new[rowSums(counts_data_NCs_new > 5) > ncol(counts_data_NCs_new)/100, ]
dim(filt_counts_NCs_new)
topgenes_new <- names(sort(gam.pval_NC, decreasing = FALSE))[1:50] 

##FITGAM
pdtime_NC_new <- slingPseudotime(curves_NCs_new)
slingLineages(curves_NCs_new)
sce_NC_new <- fitGAM(counts = as.matrix(filt_counts_NCs_new), pseudotime=pdtime_NC_new, sds = curves_NCs_new, verbose=T)
head(filt_counts_NCs_new)

plotGeneCount(curves_NCs_new, filt_counts_NCs_new, clusters = clustering_NC_new, models = sce_NC_new)

###pattern test###
patternRes_N <- patternTest(sce_NC_new)
oPat_N <- order(patternRes_N$waldStat, decreasing = TRUE)
head(rownames(patternRes_N)[oPat_N], n=100)
plotSmoothers(sce_NC_new, counts_data_NCs_new, gene = "GAD1B")
p1 <- plotGeneCount(curves_NCs_new, counts_data_NCs_new, gene = "TBR1")
p2 <- plotGeneCount(curves_NCs_new, counts_data_NCs_new, gene = "GAD1B")
p1+p2

earlyDERes_N <- earlyDETest(sce_NC_new, knots = c(4,5),l2fc = log2(1.5))
oEarly_N <- order(earlyDERes_N$waldStat, decreasing = TRUE)
head(rownames(earlyDERes_N)[oEarly_N], n=50)
plotSmoothers(sce_NC_new, counts_data_Y_NCs, gene = rownames(earlyDERes_N)[oEarly_N][37])
plotGeneCount(curves_NCs_new, counts_data_Y_NCs, gene = rownames(earlyDERes_N)[oEarly_N][37])
