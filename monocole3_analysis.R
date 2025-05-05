####################################
# Written by Dr. Noa Novershtern 
# https://github.com/hannalab/Mouse_TFSEM
# Run on R.4.3.1
#####################################
#WORKING_DIR = "<enter your working dir here>" 
#setwd(WORKING_DIR)
library(Seurat)
library(Matrix)
library(monocle3)
set.seed(1234)

##### color definition ######
pluri_color = "#FFCC99"
te_color = "#339999"
pre_color = "#CC6666"
twoc_color = "#9933CC"

##### Load Seurat object rds #####
seurat.data<-readRDS("C_48h.rds")

##### Generate cds object, pre-process and cluster #####
cds <- as.cell_data_set(seurat.data,reductions = "umap",assay = "RNA",slot="counts")
cds <- preprocess_cds(cds,method='PCA')
cds <- reduce_dimension(cds,reduction_method = 'UMAP',preprocess_method = 'PCA')
cds <- cluster_cells(cds, resolution=1e-3,reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)

##### Set metadata #####
cds$seurat_clusters <- seurat.data$seurat_clusters
cds$type =  seurat.data$type
custom_colors <- c("EPI" = pluri_color,"PrE" = pre_color, "TE" = te_color, "2CL" = twoc_color)
color_vector = cds$type
color_vector[color_vector=="TE"] = te_color
color_vector[color_vector=="PrE"] = pre_color
color_vector[color_vector=="EPI"] = pluri_color
color_vector[color_vector=="2CL"] = twoc_color

### Plot UMAP, cells are colored by type ####
plot_cells(cds,
           color_cells_by = "type",
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_roots = F,
           label_branch_points=FALSE,cell_size = 1.5,graph_label_size = 3,
           ) + geom_point(colour = color_vector,fill="black")

### Plot UMAP with pseudotime trajectory, cells are colored by cluster ####
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=T,
           label_leaves=FALSE,
           label_roots = F,
           label_branch_points=FALSE,cell_size = 0.1,graph_label_size = 10,
) 

### Set pseudotime root (cluster number may change) ###
cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 2])) # EPI

### Plot UMAP with pseudotime trajectory, cells are colored by pseudotime ####
plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=F,
           label_branch_points=FALSE,
           label_roots = F,
           trajectory_graph_color = "lightgrey",cell_size = 1.5)


