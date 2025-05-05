####################################
# Written by Dr. Noa Novershtern 
# https://github.com/hannalab/Mouse_TFSEM
# Run on R.4.3.1
#####################################
WORKING_DIR = "<enter working dir here>" 
setwd(WORKING_DIR)
library(Seurat)
library(Matrix)

##### Load reference data, including expression data and cell metadata #####
# Taken from Proks, M., Salehin, N. & Brickman, J.M. Deep learning-based models for preimplantation 
# mouse and human embryos based on single-cell RNA sequencing. Nat Methods 22, 207–216 (2025).
#####
brickman_data <- read.csv("brickman_expression_matrix.txt",header=T,sep="\t",row.names = 1)
brickman_meta <- read.csv("brickman_cell_metadata.csv",header=T)
brickman <- CreateSeuratObject(counts = brickman_data, project = "ATLAS", min.cells = 3, min.features = 200)
D_96h.data<-readRDS("D_priming.rds")
C_48h.data<-readRDS("C_48h.rds")
B_24h.data<-readRDS("B_24h.rds")
A_2i_LIF.data<-readRDS("A_2i_LIF.rds")

######### Normalization, UMAP and Clustering ##############
dims=1:10
brickman <- NormalizeData(brickman)
brickman <- FindVariableFeatures(brickman, selection.method = "vst", nfeatures = 2000)
brickman <- ScaleData(brickman, verbose = FALSE)
brickman <- RunPCA(brickman, npcs = 30, verbose = FALSE)
brickman <- RunUMAP(brickman, reduction = "pca", dims = dims,min.dist = 0.1)
brickman <- FindNeighbors(brickman, reduction = "pca", dims = dims)
brickman <- FindClusters(brickman, resolution = 0.1)#seq(0,1.2,0.1))

brickman$ct = brickman_meta$ct
brickman$stage = brickman_meta$timepoint  
brickman_colors = c("#F8766D","#E58700","#C99800","#A3A500","#B983FF","#6BB100","#619CFF",
              "#E76BF3","red3","#00BA38","#FF67A4","#00B0F6","red","#00BF7D","#00BCD8")

DimPlot(brickman, reduction = "umap", group.by="ct",label=T,pt.size = 0.5,repel=2,cols=brickman_colors) & NoAxes()

######## Combined UMAP analysis, with anchoring ############
#brickman$dataset <- "Atlas"

new_sub.list <- c(brickman,A_2i_LIF.data,B_24h.data,C_48h.data,D_96h.data)
new_sub.list <- lapply(X = new_sub.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

dims=1:10
anchors <- FindIntegrationAnchors(object.list = new_sub.list, dims = dims,anchor.features = 2000)
combined <- IntegrateData(anchorset = anchors, dims = dims)
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = dims,min.dist = 0.01)

##### Define colors #####
pluri_color = "#FFCC99"
te_color = "#339999"
pre_color = "#CC6666"
twoc_color = "#9933CC"
my_colors = c(pluri_color,pre_color,te_color,twoc_color)
my_colors_ct = c("#F8766D","#E58700","#C99800","#A3A500","#B983FF","#6BB100","#619CFF",
              "#E76BF3","red3","#00BA38","#FF67A4","#00B0F6","red","#00BF7D",my_colors,"grey")

combined$ct[combined$orig.ident=="SEM2"]= combined$type[combined$orig.ident=="SEM2"]
combined$ct[combined$ct=="2CL"]= "TwoCL"

combined$sample[combined$sample=="2i_LIF"]= "A_2i_LIF"
combined$sample[combined$sample=="24h"]= "B_24h"
combined$sample[combined$sample=="48h"]= "C_48h"
combined$sample[combined$sample=="96h"]= "D_96hh"
combined$sample[combined$orig.ident=="ATLAS"]= "ref"

DimPlot(combined, reduction = "umap", group.by="ct",label=T,pt.size = 1,repel=2,alpha = 1,cols=my_colors_ct,split.by="sample") & NoAxes()
