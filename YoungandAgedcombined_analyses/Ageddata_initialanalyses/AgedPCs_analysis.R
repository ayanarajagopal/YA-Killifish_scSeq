####Aging PCs data only processing#####
table(Idents(killi_youngandaged_integrated_final))
DimPlot(agedisoseq,reduction = "tsne", pt.size = 1.8,label=T)
Plot_Density_Joint_Only(agedisoseq, c("PCNA","GLUL (2 OF 2)"), reduction="tsne") + theme_classic()

agedPCs <- subset(agedisoseq,idents=c("InterC","NE-RG3","NGP","Astro-RG1","EPD-RG4","Astro-RG2"))
agedPCs <- RunPCA(object = agedPCs, verbose = FALSE)
ElbowPlot(agedPCs, ndims=40)

# Determine the K-nearest neighbor graph
agedPCs <- FindNeighbors(object = agedPCs,dims = 1:20)

# Determine the clusters for various resolutions                                
agedPCs <- FindClusters(object = agedPCs,resolution = 0.8)

#Idents(object = killi.integrated) <- "RNA_snn_res.0.8"
cluster_ids_0.6 <- Idents(killi_youngandaged_integrated_final)
cluster_ids_0.6

#Calc. of TSNE and UMAP
agedPCs <- RunUMAP(object = agedPCs, reduction= "pca", dims = 1:20)
agedPCs <- RunTSNE(object = agedPCs, reduction= "pca", dims = 1:20)

DimPlot(agedPCs, reduction = "tsne", pt.size = 1,label.size=8,label=T)
DefaultAssay(agedPCs) <- "RNA"
agedPCs <- NormalizeData(agedPCs, verbose = FALSE)

AgedPCs_integrated_markers <- FindAllMarkers(agedPCs, min.pct = 0.25, test.use="MAST", only.pos = F, logfc.threshold = 0.25)
head(AgedPCs_integrated_markers)
agedPCs_PCs_cluster_top200 <- AgedPCs_integrated_markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=200, wt = avg_log2FC)
youngPCs_allS_cluster_top
youngPCs_top100genes <- PCs_integrated_clustersmarkers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=100, wt = avg_log2FC)
youngPCs_top100genes
youngPCs_top10genes <- PCs_integrated_clustersmarkers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=10, wt = avg_log2FC)
youngPCs_top10genes
write_xlsx(agedPCs_PCs_cluster_top200, col_names = T, path="C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/Subclustering/agedPCs_clubbed_top200genes_res0.5_fin.xlsx")
DotPlot(agedPCs, features= genes.to.label.clusters, scale=T, cols=c("green","purple")) + RotatedAxis() + theme(axis.text.y = element_text(size=12, color="black")) + theme(axis.text.x = element_text(size=10, color="black")) + theme(axis.text.y = element_text(size=20, color="black")) +  theme(axis.line.x = element_line(size=2, color="black")) + theme(axis.line.y = element_line(size=2, color="black"))

###### Aged PC specific analyses
agedPCs_annotations <- c("Astro-RG1","Intercell.1","Intercell.2","NGP.2","EPD-RG4","NGP.1","Astro-RG2","NE-RG3","Intercell.3","Intercell.4")
names(agedPCs_annotations) <- levels(agedPCs)
agedPCs <- RenameIdents(agedPCs, agedPCs_annotations)

agedPCs$CellType <- Idents(agedPCs)
agedPCs$CellType
agedPCs$seurat_clusters <- Idents(onlyrevelantPCs)
agedPCs$seurat_clusters
table(Idents(agedPCs))

agedPC_data <- CellCycleScoring(agedPCs, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
agedPC_data[[]]
as_tibble(agedPC_data[[]]) %>%
  ggplot(aes(Phase)) + geom_bar()
agedPC_data@meta.data %>%
  group_by(CellType,Phase) %>%
  count() %>%
  group_by(CellType) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=CellType,y=percent, split,fill=Phase)) +
  geom_col() + RotatedAxis()+theme(text=element_text(size=18))

table(Idents(agedPCs))
?ReorderIdent
agedPCs$CellType <- factor(agedPCs$CellType,levels=c("NGP.1","Intercell.1","NGP.2","Intercell.2","Intercell.3","Astro-RG1","NE-RG3","Intercell.4","EPD-RG4","Astro-RG2"))
Idents(agedPCs)