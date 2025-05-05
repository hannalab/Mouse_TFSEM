####################################
# Written by Dr. Noa Novershtern 
# https://github.com/hannalab/Mouse_TFSEM
# Run on R.4.3.1
#####################################
WORKING_DIR = "<enter your working dir here>" 
setwd(WORKING_DIR)
library(cowplot)
library(dplyr)
library(patchwork)
library(Seurat)
library(Matrix)
library(ggplot2)
library(stringr)
library(pheatmap)
library(data.table)
library(ggtern)
pluri_color = "#FFCC99"
te_color = "#339999"
pre_color = "#CC6666"
twoc_color = "#9933CC"
morula_color="blue"

###### Read Seurat objects #####
D_96h.data<-readRDS("D_96h.rds")
C_48h.data<-readRDS("C_48h.rds")
B_24h.data<-readRDS("B_24h.rds")
A_2i_LIF.data<-readRDS("A_2i_LIF.rds")

##### Define annotations ######
TE = c("Cdx2","Id2","Elf5","Lgals1","Krt8","Eomes","Cldn7","Tfap2c","Gata3","Krt18")
PrE = c("Sox17","Gata6","Pdgfra","Dab2","Apoa1","Foxa2","Gata4","Serpinh1","Amn","Srgn")
pluri = c("Sox2","Pou5f1","Nanog")


##### Calculate average counts, Used in the triangle analysis #####
calculate_avg_signature <- function(model,layer="counts")
{
  count_data <- GetAssayData(model, layer = layer)
  total_counts = sum(count_data)
  TE_matrix <- count_data[TE, ] 
  pluri_matrix <- count_data[pluri, ] 
  PrE_matrix <- count_data[PrE, ] 

  model  <- AddMetaData(model, metadata = colMeans(TE_matrix)*10E6/total_counts, col.name = "TE_avg")
  model  <- AddMetaData(model, metadata = colMeans(pluri_matrix)*10E6/total_counts, col.name = "pluri_avg")
  model  <- AddMetaData(model, metadata = colMeans(PrE_matrix)*10E6/total_counts, col.name = "PrE_avg")

  return(model)
}

A_2i_LIF.data <- calculate_avg_signature(A_2i_LIF.data, layer="counts")
B_24h.data <- calculate_avg_signature(B_24h.data, layer="counts")
C_48h.data <- calculate_avg_signature(C_48h.data, layer="counts")
D_96h.data <- calculate_avg_signature(D_96h.data, layer="counts")


##### Mark the cells that express the majority of the features #####
count_cells <- function(model,layer="counts")
{
  count_data <- GetAssayData(model, layer = layer)
  TE_matrix <- count_data[TE, ] 
  pluri_matrix <- count_data[pluri, ] 
  PrE_matrix <- count_data[PrE, ] 
  
  model  <- AddMetaData(model, metadata = colSums(TE_matrix>0)>dim(TE_matrix)[1]/2, col.name = "TE_cells") 
  model  <- AddMetaData(model, metadata = colSums(pluri_matrix>0)>2, col.name = "pluri_cells") 
  model  <- AddMetaData(model, metadata = colSums(PrE_matrix>0)>dim(PrE_matrix)[1]/2, col.name = "PrE_cells") 
  return(model)
}

A_2i_LIF.data <- count_cells(A_2i_LIF.data, layer="counts")
B_24h.data <- count_cells(B_24h.data, layer="counts")
C_48h.data <- count_cells(C_48h.data, layer="counts")
D_96h.data <- count_cells(D_96h.data, layer="counts")

##### Plot seurat cluters #####
cowplot::plot_grid(ncol = 2,
  DimPlot(A_2i_LIF.data, reduction = "umap", group.by="seurat_clusters",label=T,pt.size = 0.1)& NoAxes(),
  DimPlot(B_24h.data, reduction = "umap", group.by="seurat_clusters",label=T,pt.size = 0.1)& NoAxes(),
  DimPlot(C_48h.data, reduction = "umap", group.by="seurat_clusters",label=T,pt.size = 0.1)& NoAxes(),
  DimPlot(D_96h.data, reduction = "umap", group.by="seurat_clusters",label=T,pt.size = 0.1)& NoAxes()
)

##### Color the cells that express the majority of the features #####
cowplot::plot_grid(ncol = 4,
                   FeaturePlot(A_2i_LIF.data, features = c("pluri_cells"),max.cutoff=1,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",pluri_color),keep.scale = "all") & NoAxes() & NoLegend(),
                   FeaturePlot(B_24h.data, features = c("pluri_cells"),max.cutoff=1,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",pluri_color),keep.scale = "all") & NoAxes()& NoLegend(),
                   FeaturePlot(C_48h.data, features = c("pluri_cells"),max.cutoff=1,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",pluri_color),keep.scale = "all") & NoAxes()& NoLegend(),
                   FeaturePlot(D_96h.data, features = c("pluri_cells"),max.cutoff=1,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",pluri_color),keep.scale = "all") & NoAxes()& NoLegend(),
                   FeaturePlot(A_2i_LIF.data, features = c("TE_cells"),order=F,pt.size = 0.05)+ scale_color_gradientn(colors=c("lightgrey",te_color),limits=c(0,1), na.value = "lightgray") & NoAxes() & NoLegend(),
                   FeaturePlot(B_24h.data, features = c("TE_cells"),max.cutoff=1,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",te_color),keep.scale = "all") & NoAxes()& NoLegend(),
                   FeaturePlot(C_48h.data, features = c("TE_cells"),max.cutoff=1,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",te_color),keep.scale = "all") & NoAxes()& NoLegend(),
                   FeaturePlot(D_96h.data, features = c("TE_cells"),max.cutoff=1,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",te_color),keep.scale = "all") & NoAxes()& NoLegend(),
                   FeaturePlot(A_2i_LIF.data, features = c("PrE_cells"),order=F,pt.size = 0.05)+ scale_color_gradientn(colors=c("lightgrey",pre_color),limits=c(0,1), na.value = "lightgray") & NoAxes() & NoLegend(),
                   FeaturePlot(B_24h.data, features = c("PrE_cells"),max.cutoff=1,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",pre_color),keep.scale = "all") & NoAxes()& NoLegend(),
                   FeaturePlot(C_48h.data, features = c("PrE_cells"),max.cutoff=1,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",pre_color),keep.scale = "all") & NoAxes()& NoLegend(),
                   FeaturePlot(D_96h.data, features = c("PrE_cells"),max.cutoff=1,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",pre_color),keep.scale = "all") & NoAxes()& NoLegend()
)

#### Color the cells by type ####
cowplot::plot_grid(ncol = 4,
  DimPlot(A_2i_LIF.data, reduction = "umap", group.by="type",label=F,pt.size = 1,cols=c(twoc_color,pluri_color)) & NoAxes()& NoLegend(),
  DimPlot(B_24h.data, reduction = "umap", group.by="type",label=F,pt.size = 1,cols=c(twoc_color,pluri_color,pre_color)) & NoAxes()& NoLegend(),
  DimPlot(C_48h.data, reduction = "umap", group.by="type",label=F,pt.size = 1,cols=c(twoc_color,pluri_color,pre_color,te_color)) & NoAxes()& NoLegend(),
  DimPlot(D_96h.data, reduction = "umap", group.by="type",label=F,pt.size = 1,cols=c(pluri_color,pre_color,te_color)) & NoAxes()& NoLegend()
)


##### Triangle analysis using average expression #####
plot_triangle <- function(seurat_data)
{
  columns = c("pluri_avg","TE_avg","PrE_avg")
  df <- seurat_data@meta.data[, columns]
  df <- df / (rowSums(df[])+1e-5)
  df$type <- seurat_data$type
  ggtern(data = df, aes(x = pluri_avg, y = TE_avg, z = PrE_avg, color=type)) + 
  geom_point(size=1.5) +
  scale_color_manual(values = c("TE" = te_color, "EPI" = pluri_color, "PrE" = pre_color,"2CL"=twoc_color)) + 
  theme_minimal() +
  theme(legend.position = "none")  
}
plotA <- plot_triangle(A_2i_LIF.data)
plotB <- plot_triangle(B_24h.data)
plotC <- plot_triangle(C_48h.data)
plotD <- plot_triangle(D_96h.data)

table(A_2i_LIF.data$type)
table(B_24h.data$type)
table(C_48h.data$type)
table(D_96h.data$type)

plot_grid(plotA,plotB,plotC,plotD,ncol=2)

##### Pie chart analysis #####
plot_piechart <- function(seurat_data)
{
type_df <- as.data.frame(table(seurat_data$type))
percentage = table(seurat_data$type)*100/length(seurat_data$type)
colnames(type_df) <- c("type", "count")
ggplot(type_df, aes(x = "", y = count, fill = type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_manual(values = c("TE" = te_color, "EPI" = pluri_color, "PrE" = pre_color,"2CL"=twoc_color)) 
}
plot_grid(ncol=4,
  plot_piechart(A_2i_LIF.data) & NoLegend(),
  plot_piechart(B_24h.data)& NoLegend(),
  plot_piechart(C_48h.data)& NoLegend(),
  plot_piechart(D_96h.data)& NoLegend()
  )

##### Calculate cell type distribution ##### 
as.data.frame(table(A_2i_LIF.data$type)*100/length(A_2i_LIF.data$type))
as.data.frame(table(B_24h.data$type)*100/length(B_24h.data$type))
as.data.frame(table(C_48h.data$type)*100/length(C_48h.data$type))
as.data.frame(table(D_96h.data$type)*100/length(D_96h.data$type))


##### Dotplot analysis #####
combined = readRDS("ref_all_induction_samples.rds")
pluri_markers = c("Pou5f1","Sox2","Nanog")
pre_markers = c("Gata6","Sox17","Pdgfra","Dab2")
te_markers = c("Tfap2c","Cdx2","Id2","Elf5")
twoc_markers = c("Zscan4c","Zscan4d","Tcstv3")
morula_markers=c("Pdzd3","Oasl2","Sult5a1","Dnah2os","Ifi27")

dot_features = c(pluri_markers,twoc_markers,te_markers,pre_markers)
plot_grid(ncol=1,rel_heights = c(3,4,3,2.5),
          DotPlot(D_96h.data,features = dot_features,group.by = c("type"),cluster.idents = F,assay = "RNA",scale =F,scale.min =0) + scale_y_discrete(expand = expansion(mult = c(0.3, 0.2))) & NoLegend(),
          DotPlot(C_48h.data,features = dot_features,group.by = c("type"),cluster.idents = F,assay = "RNA",scale =F,scale.min =0) + scale_y_discrete(expand = expansion(mult = c(0.3, 0.2)))& NoLegend(),  
          DotPlot(B_24h.data,features = dot_features,group.by = c("type"),cluster.idents = F,assay = "RNA",scale =F,scale.min =0) + scale_y_discrete(expand = expansion(mult = c(0.3, 0.2)))& NoLegend(), 
          DotPlot(A_2i_LIF.data,features = dot_features,group.by = c("type"),cluster.idents = F,assay = "RNA",scale =F,scale.min =0) + scale_y_discrete(expand = expansion(mult = c(0.4, 0.4)))& NoLegend()
)          

##### Feature plot for gene markers #####
plot_grid(ncol=4,
  FeaturePlot(A_2i_LIF.data, features = pluri_markers,max.cutoff=3,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",pluri_color),keep.scale = "all") & NoAxes() & NoLegend(),
  FeaturePlot(B_24h.data, features = pluri_markers,max.cutoff=3,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",pluri_color),keep.scale = "all") & NoAxes() & NoLegend(),
  FeaturePlot(C_48h.data, features = pluri_markers,max.cutoff=3,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",pluri_color),keep.scale = "all") & NoAxes() & NoLegend(),
  FeaturePlot(D_96h.data, features = pluri_markers,max.cutoff=3,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",pluri_color),keep.scale = "all") & NoAxes() & NoLegend()
)

plot_grid(ncol=4,
          FeaturePlot(A_2i_LIF.data, features = pre_markers,max.cutoff=3,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",pre_color),keep.scale = "all") & NoAxes() & NoLegend(),
          FeaturePlot(B_24h.data, features = pre_markers,max.cutoff=3,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",pre_color),keep.scale = "all") & NoAxes() & NoLegend(),
          FeaturePlot(C_48h.data, features = pre_markers,max.cutoff=3,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",pre_color),keep.scale = "all") & NoAxes() & NoLegend(),
          FeaturePlot(D_96h.data, features = pre_markers,max.cutoff=3,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",pre_color),keep.scale = "all") & NoAxes() & NoLegend()
)
FeaturePlot(D_96h.data, features = "Ttr",max.cutoff=3,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",pre_color),keep.scale = "all") & NoAxes()# & NoLegend()

plot_grid(ncol=4,
          FeaturePlot(A_2i_LIF.data, features = te_markers,max.cutoff=3,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",te_color),keep.scale = "all") & NoAxes() & NoLegend(),
          FeaturePlot(B_24h.data, features = te_markers,max.cutoff=3,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",te_color),keep.scale = "all") & NoAxes() & NoLegend(),
          FeaturePlot(C_48h.data, features = te_markers,max.cutoff=3,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",te_color),keep.scale = "all") & NoAxes() & NoLegend(),
          FeaturePlot(D_96h.data, features = te_markers,max.cutoff=3,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",te_color),keep.scale = "all") & NoAxes() & NoLegend()
)

plot_grid(ncol=4,
          FeaturePlot(A_2i_LIF.data, features = twoc_markers,max.cutoff=3,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",twoc_color),keep.scale = "all") & NoAxes() & NoLegend(),
          FeaturePlot(B_24h.data, features = twoc_markers,max.cutoff=3,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",twoc_color),keep.scale = "all") & NoAxes() & NoLegend(),
          FeaturePlot(C_48h.data, features = twoc_markers,max.cutoff=3,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",twoc_color),keep.scale = "all") & NoAxes() & NoLegend(),
          FeaturePlot(D_96h.data, features = twoc_markers,max.cutoff=3,min.cutoff = 0,ncol=1,order=T,cols=c("lightgrey",twoc_color),keep.scale = "all") & NoAxes() & NoLegend()
)

plot_grid(ncol=4,
          FeaturePlot(A_2i_LIF.data, features = morula_markers,max.cutoff=3,min.cutoff = 0.5,ncol=1,order=T,cols=c("lightgrey",morula_color),keep.scale = "all") & NoAxes() & NoLegend(),
          FeaturePlot(B_24h.data, features = morula_markers,max.cutoff=3,min.cutoff = 0.5,ncol=1,order=T,cols=c("lightgrey",morula_color),keep.scale = "all") & NoAxes() & NoLegend(),
          FeaturePlot(C_48h.data, features = morula_markers,max.cutoff=3,min.cutoff = 0.5,ncol=1,order=T,cols=c("lightgrey",morula_color),keep.scale = "all") & NoAxes() & NoLegend(),
          FeaturePlot(D_96h.data, features = morula_markers,max.cutoff=3,min.cutoff = 0.5,ncol=1,order=T,cols=c("lightgrey",morula_color),keep.scale = "all") & NoAxes() & NoLegend()
)


##### Blend analysis #####
overlap_percentage <- function(model,gene1,gene2,layer="counts")
{
    count_data <- GetAssayData(model, layer = layer)
    cut <- table(count_data[gene1,]>0 & count_data[gene2,]>0)["TRUE"]
    union <- table(count_data[gene1,]>0 | count_data[gene2,]>0)["TRUE"]
    print(paste(gene1,gene2,cut,union,cut*100/union,sep=" "))
}

overlap_percentage(A_2i_LIF.data,"Gata4","Pou5f1")
overlap_percentage(A_2i_LIF.data,"Gata6","Pou5f1")
overlap_percentage(A_2i_LIF.data,"Cdx2","Pou5f1")

overlap_percentage(A_2i_LIF.data,"Nanog","Pou5f1")
overlap_percentage(A_2i_LIF.data,"Gata4","Nanog")
overlap_percentage(A_2i_LIF.data,"Gata6","Nanog")
overlap_percentage(A_2i_LIF.data,"Cdx2","Nanog")
overlap_percentage(A_2i_LIF.data,"Elf5","Nanog")

overlap_percentage(A_2i_LIF.data,"Gata4","Zscan4d")
overlap_percentage(A_2i_LIF.data,"Gata6","Zscan4d")
overlap_percentage(A_2i_LIF.data,"Pdgfra","Zscan4d")
overlap_percentage(A_2i_LIF.data,"Cdx2","Zscan4d")
overlap_percentage(A_2i_LIF.data,"Elf5","Zscan4d")

overlap_percentage(B_24h.data,"Gata4","Pou5f1")
overlap_percentage(B_24h.data,"Gata6","Pou5f1")
overlap_percentage(B_24h.data,"Cdx2","Pou5f1")

overlap_percentage(B_24h.data,"Nanog","Pou5f1")
overlap_percentage(B_24h.data,"Gata4","Nanog")
overlap_percentage(B_24h.data,"Gata6","Nanog")
overlap_percentage(B_24h.data,"Cdx2","Nanog")
overlap_percentage(B_24h.data,"Elf5","Nanog")

overlap_percentage(B_24h.data,"Gata4","Zscan4d")
overlap_percentage(B_24h.data,"Gata6","Zscan4d")
overlap_percentage(B_24h.data,"Pdgfra","Zscan4d")
overlap_percentage(B_24h.data,"Cdx2","Zscan4d")
overlap_percentage(B_24h.data,"Elf5","Zscan4d")

overlap_percentage(C_48h.data,"Gata4","Pou5f1")
overlap_percentage(C_48h.data,"Gata6","Pou5f1")
overlap_percentage(C_48h.data,"Cdx2","Pou5f1")

overlap_percentage(C_48h.data,"Nanog","Pou5f1")
overlap_percentage(C_48h.data,"Gata4","Nanog")
overlap_percentage(C_48h.data,"Gata6","Nanog")
overlap_percentage(C_48h.data,"Cdx2","Nanog")
overlap_percentage(C_48h.data,"Elf5","Nanog")

overlap_percentage(C_48h.data,"Gata4","Zscan4d")
overlap_percentage(C_48h.data,"Gata6","Zscan4d")
overlap_percentage(C_48h.data,"Pdgfra","Zscan4d")
overlap_percentage(C_48h.data,"Cdx2","Zscan4d")
overlap_percentage(C_48h.data,"Elf5","Zscan4d")

overlap_percentage(D_96h.data,"Gata4","Pou5f1")
overlap_percentage(D_96h.data,"Gata6","Pou5f1")
overlap_percentage(D_96h.data,"Cdx2","Pou5f1")

overlap_percentage(D_96h.data,"Nanog","Pou5f1")
overlap_percentage(D_96h.data,"Gata4","Nanog")
overlap_percentage(D_96h.data,"Gata6","Nanog")
overlap_percentage(D_96h.data,"Cdx2","Nanog")
overlap_percentage(D_96h.data,"Elf5","Nanog")

overlap_percentage(D_96h.data,"Gata4","Zscan4d")
overlap_percentage(D_96h.data,"Gata6","Zscan4d")
overlap_percentage(D_96h.data,"Pdgfra","Zscan4d")
overlap_percentage(D_96h.data,"Cdx2","Zscan4d")
overlap_percentage(D_96h.data,"Elf5","Zscan4d")

overlap_percentage(A_2i_LIF.data,"Gata4","Cdx2")
overlap_percentage(A_2i_LIF.data,"Gata6","Cdx2")
overlap_percentage(A_2i_LIF.data,"Pdgfra","Cdx2")

overlap_percentage(B_24h.data,"Gata4","Cdx2")
overlap_percentage(B_24h.data,"Gata6","Cdx2")
overlap_percentage(B_24h.data,"Pdgfra","Cdx2")

overlap_percentage(C_48h.data,"Gata4","Cdx2")
overlap_percentage(C_48h.data,"Gata6","Cdx2")
overlap_percentage(C_48h.data,"Pdgfra","Cdx2")

overlap_percentage(D_96h.data,"Gata4","Cdx2")
overlap_percentage(D_96h.data,"Gata6","Cdx2")
overlap_percentage(D_96h.data,"Pdgfra","Cdx2")

plot_grid(ncol=1,
        FeaturePlot(A_2i_LIF.data,features = c("Pou5f1","Gata4"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
        FeaturePlot(A_2i_LIF.data,features = c("Pou5f1","Gata6"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
        FeaturePlot(A_2i_LIF.data,features = c("Pou5f1","Cdx2"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
        FeaturePlot(A_2i_LIF.data,features = c("Pou5f1","Nanog"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
        FeaturePlot(A_2i_LIF.data,features = c("Nanog","Gata4"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
        FeaturePlot(A_2i_LIF.data,features = c("Nanog","Gata6"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
        FeaturePlot(A_2i_LIF.data,features = c("Nanog","Cdx2"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
        FeaturePlot(A_2i_LIF.data,features = c("Nanog","Elf5"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes()
)

plot_grid(ncol=1,
          FeaturePlot(B_24h.data,features = c("Pou5f1","Gata4"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(B_24h.data,features = c("Pou5f1","Gata6"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(B_24h.data,features = c("Pou5f1","Cdx2"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(B_24h.data,features = c("Pou5f1","Nanog"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(B_24h.data,features = c("Nanog","Gata4"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(B_24h.data,features = c("Nanog","Gata6"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(B_24h.data,features = c("Nanog","Cdx2"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(B_24h.data,features = c("Nanog","Elf5"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes()
)

plot_grid(ncol=1,
          FeaturePlot(C_48h.data,features = c("Pou5f1","Gata4"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(C_48h.data,features = c("Pou5f1","Gata6"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(C_48h.data,features = c("Pou5f1","Cdx2"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(C_48h.data,features = c("Pou5f1","Nanog"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(C_48h.data,features = c("Nanog","Gata4"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(C_48h.data,features = c("Nanog","Gata6"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(C_48h.data,features = c("Nanog","Cdx2"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(C_48h.data,features = c("Nanog","Elf5"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes()
)

plot_grid(ncol=1,
          FeaturePlot(D_96h.data,features = c("Pou5f1","Gata4"),blend = T,order=T,cols=c("lightgrey","red2","green3")) & NoAxes(),
          FeaturePlot(D_96h.data,features = c("Pou5f1","Gata6"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(D_96h.data,features = c("Pou5f1","Cdx2"),blend = T,order=T,cols=c("lightgrey","red2","green3")) & NoAxes(),
          FeaturePlot(D_96h.data,features = c("Pou5f1","Nanog"),blend = T,order=T,cols=c("lightgrey","red2","green3")) & NoAxes(),
          FeaturePlot(D_96h.data,features = c("Nanog","Gata4"),blend = T,order=T,cols=c("lightgrey","red2","green3")) & NoAxes(),
          FeaturePlot(D_96h.data,features = c("Nanog","Gata6"),blend = T,order=T,cols=c("lightgrey","red2","green3")) & NoAxes(),
          FeaturePlot(D_96h.data,features = c("Nanog","Cdx2"),blend = T,order=T,cols=c("lightgrey","red2","green3")) & NoAxes(),
          FeaturePlot(D_96h.data,features = c("Nanog","Elf5"),blend = T,order=T,cols=c("lightgrey","red2","green3")) & NoAxes()
)

plot_grid(ncol=1,
          FeaturePlot(A_2i_LIF.data,features = c("Cdx2","Gata4"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(A_2i_LIF.data,features = c("Cdx2","Gata6"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(A_2i_LIF.data,features = c("Cdx2","Pdgfra"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes()
          )
plot_grid(ncol=1,
          FeaturePlot(A_2i_LIF.data,features = c("Zscan4d","Gata4"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(A_2i_LIF.data,features = c("Zscan4d","Gata6"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(A_2i_LIF.data,features = c("Zscan4d","Pdgfra"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(A_2i_LIF.data,features = c("Zscan4d","Cdx2"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(A_2i_LIF.data,features = c("Zscan4d","Elf5"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes()
)

plot_grid(ncol=1,
          FeaturePlot(B_24h.data,features = c("Cdx2","Gata4"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(B_24h.data,features = c("Cdx2","Gata6"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(B_24h.data,features = c("Cdx2","Pdgfra"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes()
          )
plot_grid(ncol=1,
          FeaturePlot(B_24h.data,features = c("Zscan4d","Gata4"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(B_24h.data,features = c("Zscan4d","Gata6"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(B_24h.data,features = c("Zscan4d","Pdgfra"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(B_24h.data,features = c("Zscan4d","Cdx2"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(B_24h.data,features = c("Zscan4d","Elf5"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes()
)

plot_grid(ncol=1,
          FeaturePlot(C_48h.data,features = c("Cdx2","Gata4"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(C_48h.data,features = c("Cdx2","Gata6"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(C_48h.data,features = c("Cdx2","Pdgfra"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes()
          )
plot_grid(ncol=1,
          FeaturePlot(C_48h.data,features = c("Zscan4d","Gata4"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(C_48h.data,features = c("Zscan4d","Gata6"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(C_48h.data,features = c("Zscan4d","Pdgfra"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(C_48h.data,features = c("Zscan4d","Cdx2"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(C_48h.data,features = c("Zscan4d","Elf5"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes()
)

plot_grid(ncol=1,
          FeaturePlot(D_96h.data,features = c("Cdx2","Gata4"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(D_96h.data,features = c("Cdx2","Gata6"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(D_96h.data,features = c("Cdx2","Pdgfra"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes()
        )
plot_grid(ncol=1,
          FeaturePlot(D_96h.data,features = c("Zscan4d","Gata4"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(D_96h.data,features = c("Zscan4d","Gata6"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(D_96h.data,features = c("Zscan4d","Pdgfra"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(D_96h.data,features = c("Zscan4d","Cdx2"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes(),
          FeaturePlot(D_96h.data,features = c("Zscan4d","Elf5"),blend = T,order=T,cols=c("lightgrey","red2","green3"))& NoAxes()
)
