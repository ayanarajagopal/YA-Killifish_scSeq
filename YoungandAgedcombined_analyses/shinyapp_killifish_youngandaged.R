library(Seurat)
library(ShinyCell)
library(shinyapps)

##Shiny App for single cell sequencing analysis
reqPkg = c("shiny", "shinyhelper", "data.table", "Matrix", "DT", "hdf5r", 
           "reticulate", "ggplot2", "gridExtra", "magrittr", "ggdendro")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}

devtools::install_github("SGDDNB/ShinyCell")
devtools::install_github("shinyapps")

# Load whole dataset 1 (~2 GB)
yfish_allS_allcells = readRDS("C:/Users/u0129074/Documents/ScRNAseq_paper/Final_versions/CellReports/GEO_submission_Yfish_Sep2022/Codes_scseq_submission/youngkillifish_allcells.rds") 
yfish_allS_allcells
scConf1_allS = createConfig(yfish_allS_allcells)
showLegend(scConf1_allS)
showOrder(scConf1_allS)

#makeShinyApp(yfish_data, scConf2, gene.mapping = TRUE,
#shiny.title = "Kfish all cells") 
remove(scConf1_allS)
scConf1_allS = delMeta(scConf1_allS, c("orig.ident", "SCT_snn_res.0.8"))
scConf1_allS = modMetaName(scConf1_allS, meta.to.mod = c("nUMI", "nGene", "percent.mt"), 
                           new.name = c("No. UMIs", "No. detected genes",
                                        "% MT genes"))
?modMetaName
scConf1_allS = modColours(scConf1_allS, meta.to.mod = "library", 
                          new.colours= c("black", "darkorange", "blue", "pink2"))
makeShinyFiles(yfish_allS_allcells, scConf1_allS, gex.assay = "RNA", gex.slot = "data",
               gene.mapping = TRUE,shiny.prefix = "sc1",
               shiny.dir = "C:/Users/u0129074/Documents/shinyAppMulti_allS1/",
               default.gene1 = "PCNA", default.gene2 = "CX43",
               default.multigene = c("PCNA","SOX2","CX43","SLC1A2"),
               default.dimred = c("tSNE1", "tSNE2"))

yfish_PCs_data <- readRDS("C:/Users/u0129074/Documents/ScRNAseq_paper/Final_versions/CellReports/GEO_submission_Yfish_Sep2022/Codes_scseq_submission/youngkillifish_PCs_mar2023.rds")

scConf2 = createConfig(yfish_PCs_data)
scConf2 = delMeta(scConf2, c("orig.ident", "RNA_snn_res.0.5"))
showLegend(scConf2)
showOrder(scConf2)
scConf2
scConf2 = modMetaName(scConf2,  meta.to.mod = c("nUMI", "nGene", "percent.mt"), 
                      new.name = c("No. UMIs", "No. detected genes", "% MT genes"))
scConf2 = modColours(scConf2, meta.to.mod = "library", 
                     new.colours= c("black", "blue", "purple"))
showOrder(scConf2)
makeShinyFiles(yfish_PCs_data, scConf2, gex.assay = "RNA", gex.slot = "data",
               gene.mapping = TRUE,shiny.prefix = "sc2",
               shiny.dir = "C:/Users/u0129074/Documents/shinyAppMulti_allS1/",
               default.gene1 = "PCNA", default.gene2 = "CX43",
               default.multigene = c("PCNA","SOX2","CX43","SLC1A2B"),
               default.dimred = c("tSNE1", "tSNE2"))

yfish_NCs_PCs_data <- readRDS("C:/Users/u0129074/Documents/ScRNAseq_paper/Final_versions/CellReports/GEO_submission_Yfish_Sep2022/Codes_scseq_submission/youngkillifish_NCPCs_mar2023.Rds") 
scConf3 = createConfig(yfish_NCs_PCs_data)
showLegend(scConf3)
showOrder(scConf3)
scConf3 = modMetaName(scConf3,  meta.to.mod = c("nUMI", "nGene", "percent.mt"), 
                      new.name = c("No. UMIs", "No. detected genes", "% MT genes"))

??makeShinyFiles
makeShinyFiles(yfish_NCs_PCs_data, scConf3, gex.assay = "RNA", gex.slot = "data",
               gene.mapping = TRUE, shiny.prefix = "sc3",
               shiny.dir = "C:/Users/u0129074/Documents/shinyAppMulti_allS/",
               default.gene1 = "PCNA", default.gene2 = "ELAVL3",
               default.multigene = c("PCNA","SLC17A6A","CX43","MAP2","MEX3A","EOMESA"),
               default.dimred = c("tSNE1", "tSNE2"))


yfish_afish_allcells_data <- readRDS("C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/aged/killi_young_aged_merged.rds")
rm(yfish_afish_allcells_data)
scConf4 = createConfig(yfish_afish_allcells_data)
showLegend(scConf4)
showOrder(scConf4)

#makeShinyApp(yfish_data, scConf2, gene.mapping = TRUE,
#shiny.title = "Kfish all cells") 

scConf4 = delMeta(scConf4, c("orig.ident"))
scConf4 = modMetaName(scConf4, meta.to.mod = c("nUMI", "nGene", "percent.mt"), 
                           new.name = c("No. UMIs", "No. detected genes",
                                        "% MT genes"))
?modColours
showOrder(scConf4)
scConf4 = modColours(scConf4, meta.to.mod = "library", 
                          new.colours= c("black", "darkorange", "blue", "pink2"))
View(scConf4)
checkConfig(scConf4,yfish_afish_allcells_data)
?modColours
makeShinyFiles(yfish_afish_allcells_data, scConf4, gex.assay = "integrated", gex.slot = "data",
               gene.mapping = TRUE,shiny.prefix = "sc4",
               shiny.dir = "C:/Users/u0129074/Documents/shinyAppMulti_allS1/",
               default.gene1 = "PCNA", default.gene2 = "CX43",
               default.multigene = c("PCNA","SOX2","CX43","SLC1A2B"),
               default.dimred = c("tSNE1", "tSNE2"))

######Young and aged progenitors######
yfish_afish_allPCs <- readRDS("C:/Users/u0129074/Documents/ScRNAseq_paper/Extranaalyses/aged/killi_young_aged_PCs_merged.rds")

scConf5 = createConfig(yfish_afish_allPCs)
showLegend(scConf5)
showOrder(scConf5)

#makeShinyApp(yfish_data, scConf2, gene.mapping = TRUE,
#shiny.title = "Kfish all cells") 

#scConf5 = delMeta(scConf5, c("orig.ident"))
scConf5 = modMetaName(scConf5, meta.to.mod = c("nUMI", "nGene", "percent.mt"), 
                      new.name = c("No. UMIs", "No. detected genes",
                                   "% MT genes"))
?modColours
showOrder(scConf4)
scConf5 = modColours(scConf5, meta.to.mod = "library", 
                     new.colours= c("black", "darkorange", "blue", "pink2"))
View(scConf4)
checkConfig(scConf5,yfish_afish_allPCs)
?modColours
makeShinyFiles(yfish_afish_allPCs, scConf5, gex.assay = "integrated", gex.slot = "data",
               gene.mapping = TRUE,shiny.prefix = "sc5",
               shiny.dir = "C:/Users/u0129074/Documents/shinyAppMulti_allS1/",
               default.gene1 = "PCNA", default.gene2 = "CX43",
               default.multigene = c("PCNA","SOX2","CX43","SLC1A2B"),
               default.dimred = c("tSNE1", "tSNE2"))

citation = list(
  author  = "Ayana R., Caroline Z. et al.",
  title   = "Single cell sequencing unravels the cellular diversity that shapes neuro- and gliogenesis in the fast aging killifish (N. furzeri) brain")
makeShinyCodesMulti(shiny.title = "Killifish telencephalon", shiny.footnotes = citation,
  shiny.prefix = c("sc1", "sc2","sc4","sc5"),
  shiny.headers = c("Youngkillifish_allcells", "Young_PCs","Youngvsaged_allcells","Youngvsaged_allPCs"), 
  shiny.dir = "C:/Users/u0129074/Documents/shinyAppMulti_allS1/") 
?makeShinyCodesMulti
rsconnect::setAccountInfo(name='ayana-rajagopal',
                          token='6AA402F0AD4EC3DCCE3EED515B96E6FD',
                          secret='FTZtvlYhUqR72GC41gHyl5pq83H0nAB53kACIMMH')
rsconnect::deployApp('C:/Users/u0129074/Documents/shinyAppMulti_allS1/')
y