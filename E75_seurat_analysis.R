####################################
# Written by Dr. Noa Novershtern 
# https://github.com/hannalab/Mouse_TFSEM
# Run on R.4.3.1
#####################################
#WORKING_DIR = "<enter your working dir here>" 
#setwd(WORKING_DIR)
library(Seurat)
library(Matrix)
library(ggplot2)
library(stringr)

###### Define my colors ########
# EPI colors
cluster8_col = "#E58700"
cluster5_col = "#FFCC99"
cluster10_col = "orange"

# TE colors
cluster1_col = "#00BA38"
cluster7_col = "#6BB100"
cluster6_col = "#339999"
cluster3_col = "#A3A500"
cluster14_col = "darkgreen"
  
# PrE colors
cluster0_col = "#F8766D"
cluster11_col = "red3"
cluster12_col = "red"
cluster13_col = "#CC6666"

# blood colors
cluster9_col = "#619CFF"
cluster15_col = "#00B0F6"

# mesoderm colors
cluster2_col = "#B983FF"
cluster4_col = "#E76BF3"

my_colors_ct = c(cluster0_col,cluster1_col,cluster2_col,cluster3_col,cluster4_col,cluster5_col,cluster6_col,cluster7_col,cluster8_col,cluster9_col,cluster10_col,cluster11_col,cluster12_col,cluster13_col,cluster14_col,cluster15_col)#,cluster16_col,cluster17_col)

##### Load data after CellRanger #####
TFSEM_vs_Nat_E75.data <- Read10X(data.dir = "<path to data>/outs/count/filtered_feature_bc_matrix/")
TFSEM.data <- CreateSeuratObject(counts = TFSEM_vs_Nat_E75.data, project = "SEM2", min.cells = 3, min.features = 200)

#####Add sample names ######
TFSEM.data$sample<-as.numeric(str_extract(colnames(TFSEM.data),"\\d+"))
type<- TFSEM.data$sample
type<-ifelse(type==1,"Head_fold",type)
type<-ifelse(type==2,"Late_streak",type)
type<-ifelse(type==3,"Nat_E75",type)
type<-ifelse(type==4,"TFSEM_E75",type)
TFSEM.data$sample<-type
###### QC ###### 

TFSEM.data[["percent.mt"]] <- PercentageFeatureSet(TFSEM.data, pattern = "^mt-")
VlnPlot(TFSEM.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0,same.y.lims=F,log = T,split.by = "sample",group.by="sample")+theme(legend.position = "none")+xlab("")
plot1 <- FeatureScatter(TFSEM.data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(TFSEM.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

##### Filter ######
TFSEM.data_sub <- subset(TFSEM.data, subset = nFeature_RNA > 3000 & nFeature_RNA <8000 & percent.mt < 10)
rm(TFSEM.data)

##### Normalization ######
new_sub.list <- SplitObject(TFSEM.data_sub, split.by = "sample")
new_sub.list <- lapply(X = new_sub.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

##### Combined UMAP analysis, with anchoring ######
dims=1:15
anchors <- FindIntegrationAnchors(object.list = new_sub.list, dims = dims,anchor.features = 2000)
combined <- IntegrateData(anchorset = anchors, dims = dims)
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = dims)
combined <- FindNeighbors(combined, reduction = "pca", dims = dims)
combined <- FindClusters(combined, resolution = 0.5)

###### Visualization #####
# To upload the exact seurat object used in the paper: 
# combined <- readRDS("combined_SEM_Nat_E75.rds")

DimPlot(combined, reduction = "umap", group.by="seurat_clusters",split.by="sample",label=T,pt.size = 1,cols = my_colors_ct) & NoAxes() 
VlnPlot(combined, features = c("nCount_RNA"), ncol = 1,pt.size = 0,same.y.lims=F,log = T,group.by="seurat_clusters")+theme(legend.position = "none",plot.title = element_blank())+xlab("")
FeaturePlot(combined,"nCount_RNA")
FeaturePlot(combined,"nFeature_RNA")

##### Find markers ######
TFSEM_7_5_markers = FindAllMarkers(combined,logfc.threshold = 0.5,only.pos =T )
write.csv(TFSEM_7_5_markers,"TFSEM7.5_markers.csv")


##### EPC analysis #########
TGC_progenitors =c("Adm","Fosl1","Hand1","Trpm5","Maged2","Prl5a1")
Uncommited_EPC =c("Chsy1","Gjb3","Krt19","Lgals1","Cald1","Ctsl")
Intermediate_Chorion =c("Ascl2","Fgfr2","Cited1","Gjb3","Ndrg1","Irx2","Irx3")
Chorion_progenitors =c("Sox3","Dusp6","Nat8l","Bmp4","Sox2","Esrrb","Eomes")
general_placenta =c("Rhox5")
SpA_TGC =c("Ctla2a","Pecam1","Ramp3","Igfbp7","Nos3")
SpT_Gly =c("Dlx3","Car2","Ncam1","Pcdh12","Tpbpa")
p_TGC =c("Star","Serpinb9d","Hsd3b6","Rhox6","Cts7","Prl2c3")
Chorion =c("Irx4","Esx1","Id1","Id3","Phlda2","Klhl13")


combined$TGC_progenitors_mean <- rowMeans(FetchData(combined, vars = TGC_progenitors))
combined$Uncommited_EPC_mean <- rowMeans(FetchData(combined, vars = Uncommited_EPC))
combined$Intermediate_Chorion_mean <- rowMeans(FetchData(combined, vars = Intermediate_Chorion))
combined$Chorion_mean <- rowMeans(FetchData(combined, vars = Chorion))
combined$Chorion_progenitors_mean <- rowMeans(FetchData(combined, vars = Chorion_progenitors))
combined$General_placenta_mean <- rowMeans(FetchData(combined, vars = general_placenta))
combined$SpA_TGC_mean <- rowMeans(FetchData(combined, vars = SpA_TGC))
combined$p_TGC_mean <- rowMeans(FetchData(combined, vars = p_TGC))
combined$SpT_Gly_mean <- rowMeans(FetchData(combined, vars = SpT_Gly))

par(mfrow = c(3, 3))
EPC_groups_mean = c("Uncommited_EPC_mean","General_placenta_mean","Chorion_progenitors_mean","Intermediate_Chorion_mean",
                    "Chorion_mean","TGC_progenitors_mean","SpA_TGC_mean","p_TGC_mean","SpT_Gly_mean")

FeaturePlot(combined,features = EPC_groups_mean,cols=c("grey","red","brown3","brown"),min.cutoff = 0) & NoAxes() & NoLegend()

##### Dotplot analysis ########
Viseral_end = c("Ttr","Hnf4a","Lrp2","Afp","Apoa1") # cluster 0, 12, 13
Parietal_end = c("Lama1","Lamb1","Dab2","Sparc","Snai1") # cluster 11
DE = c("Cer1","Sox17","Hesx1","Foxa1","Lhx1","Sox7") # cluster 10
Notochord = c("Noto") # cluster 10
Early_PrE = c("Gata4","Gata6","Foxa2","Pdgfra") # cluster 11
Uncommited_EPC =c("Gjb3","Krt19","Lgals1") # clusters 3,6
Chorion_prog = c("Cdx2","Eomes","Tfap2c","Elf5","Esrrb","Bmp4","Nat8l") # clusters 1,14
Chorion = c("Gcm1","Ascl2","Plet1") # cluster 7
Intermediate_Chorion =c("Fgfr2","Ndrg1","Irx2","Irx3") # cluster 7
TGC_prog =c("Adm","Fosl1","Hand1","Trpm5","Maged2","Prl5a1") # cluster 3
Amnion = c("Postn","Hoxc8","Tdo2","Vim") # cluster 4
SpT_Gly =c("Dlx3","Ncam1","Pcdh12","Tpbpa") # cluster 6
Nascent_mesoderm = c("Mixl1","T","Fgf8","Mesp1") # cluster 2
Epiblast = c("Otx2","Sox2","Zic2","Zic3","Pou5f1") # cluster 5, 8
Blood_prog = c("Runx1","Gata2","Lmo2","Tal1") # cluster 9
Erythroids = c("Hbb-y","Flt3","Spi1") # cluster 15

dotblot_epiblast_features = c(Nascent_mesoderm,Epiblast,Notochord,Blood_prog,Erythroids)
dotblot_endoderm_features = c(Viseral_end,Parietal_end,DE,Early_PrE)
dotblot_epc_features = c(Uncommited_EPC,Chorion_prog,Chorion,Intermediate_Chorion,TGC_prog,Amnion,SpT_Gly)
dotblot_features = c(dotblot_endoderm_features,dotblot_epc_features,dotblot_epiblast_features)

combined$order <- as.character(combined$seurat_clusters)
combined$order[combined$seurat_clusters %in% c("0","12","13")] <- 
  paste0("A_", combined$seurat_clusters[combined$seurat_clusters %in% c("0","12","13")])
combined$order[combined$seurat_clusters %in% c("10","11")] <- 
  paste0("B_", combined$seurat_clusters[combined$seurat_clusters %in% c("10","11")])
combined$order[combined$seurat_clusters %in% c("14","1","3","4","6","7")] <- 
  paste0("C_", combined$seurat_clusters[combined$seurat_clusters %in% c("14","1","3","4","6","7")])
combined$order[combined$seurat_clusters %in% c("5","8","2","4")] <- 
  paste0("D_", combined$seurat_clusters[combined$seurat_clusters %in% c("5","8","2","4")])
combined$order[combined$seurat_clusters %in% c("9","15")] <- 
  paste0("E_", combined$seurat_clusters[combined$seurat_clusters %in% c("9","15")])

DotPlot(combined,features = dotblot_features,group.by = c("order"),cluster.idents = F,assay = "RNA",scale =T,scale.min =0) + scale_y_discrete(expand = expansion(mult = c(0.04, 0.04))) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size =8))

#saveRDS(combined,"combined_SEM_Nat_E75.rds")
