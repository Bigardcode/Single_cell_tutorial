#######load libary**********************************************************

library(multtest)
library(remotes)
library(installr)
library(tidyverse)
library(ggplot2)
library(data.table)
library(Matrix)
library(cowplot)
library(devtools)
library(gplots)
library(sleepwalk)
library(SCINA)
library(dplyr)
library(Matrix)
library(e1071)
library(rlang)
library(patchwork)
library(SCINA)
library(Seurat)
library(ggplot2)
library(grid)
library(patchwork)
library(ggplot2)
library(patchwork)
library(stats)
library(grid)
library(graphics)
library(utils)
library(grDevices)
library(graphics)
library(utils)
library(gtable)





#----------------------------------------------------------------------------------

theme_complete_bw <- function(base_size = 24, base_family = "") 
{
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      axis.line =         element_blank(),
      axis.text.x =       element_text(size = base_size * 0.8 , lineheight = 0.9, colour = "black", vjust = 1),
      axis.text.y =       element_text(size = base_size * 0.8, lineheight = 0.9, colour = "black", hjust = 1),
      axis.ticks =        element_line(colour = "black"),
      axis.title.x =      element_text(size = base_size, vjust = 0.5),
      axis.title.y =      element_text(size = base_size, angle = 90, vjust = 0.5),
      axis.ticks.length = unit(0.15, "cm"),
      axis.ticks.margin = unit(0.1, "cm"),
      
      legend.background = element_rect(colour=NA), 
      legend.key =        element_rect(fill =NA, colour = "black", size = 0.25),
      legend.key.size =   unit(1.5, "lines"),
      legend.text =       element_text(size = base_size * 0.7),
      legend.title =      element_text(size = base_size * 0.8),
      legend.position =   "top",
      
      panel.background = element_rect(fill = "white", colour = NA), 
      panel.border =     element_rect(fill = NA, colour = "black", size=2), 
      panel.grid.major = element_line(colour = NA, size = 0.2), #"grey"
      panel.grid.minor = element_line(colour = NA, size = 0.5), #"grey"
      panel.margin =     unit(0.25, "lines"),
      
      strip.background = element_rect(fill = NA, colour = NA), 
      strip.text.x =     element_text(colour = "black", size = base_size * 0.8),
      strip.text.y =     element_text(colour = "black", size = base_size * 0.8, angle = +90),
      
      plot.background =  element_rect(colour = NA, fill = "white"),
      plot.title =       element_text(size = base_size*.8),
      plot.margin =      unit(c(1, 1, .5, .5), "lines"))
}

#---------------------------------------------------------------------------------------



######## Strat the pipeline #####################################################



#Step 1. Load datasets and Create a Seurat object*******************************

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Step 1. Create a Seurat object************************************************


setwd("D:/black board/sngle-cell/workflow_3")

library(Matrix)
counts <- readMM("matrix.mtx.gz")
barcodes <- read.table("barcodes.tsv.gz", stringsAsFactors=F)[,1]
features <- read.csv("features.tsv.gz", stringsAsFactors=F, sep="\t", header=F)
rownames(counts) <- make.unique(features[,2])
colnames(counts) <- barcodes


dim(counts) # report number of genes (rows) and number of cells (columns)


# Lets examine a few genes in the first thirty cells
counts[c("CD3D", "TCL1A", "MS4A1"), 1:30]

dense.size <- object.size(as.matrix(counts))
dense.size


sparse.size <- object.size(counts)
sparse.size

dense.size/sparse.size






#########CreateSeuratObject---------------------------------------------------

Zari <- CreateSeuratObject(counts, project="DS1")

head(Zari@meta.data, 6)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Step 2. Quality control********************************************************
#Pre-processing workflow


# Show QC metrics for the first 5 cells
head(Zari@meta.data, 5)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Zari[["percent.mt"]] <- PercentageFeatureSet(Zari, pattern = "^MT-")


Zari[["percent.rb"]] <- PercentageFeatureSet(Zari, pattern = "^RP[SL]")

head(Zari[[]])




###Identify cells that fail quality control using 3 metrics: 1) UMI counts, 2) 
#numberof detected genes, 3) percent mitochondrial genes

# Visualize QC metrics as a violin plot
VlnPlot(Zari, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

VlnPlot(Zari, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"), ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))




# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(Zari, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Zari, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(Zari, feature1 = "nCount_RNA", feature2 = "percent.rb")
plot1 + plot2 + plot3


#Let’s plot metadata only for cells that pass tentative QC:

Zari <- subset(Zari, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

Zari <- subset(Zari, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Step 3. Normalizing the data***************************************************

Zari <- NormalizeData(Zari)


#Identification of highly variable features (feature selection)

#selection.method = "vst"
Zari <- FindVariableFeatures(Zari, nfeatures = 3000)


# Identify the 10 most highly variable genes

# plot variable features with and without labels
top_features <- head (VariableFeatures(Zari), 20)
plot1 <- VariableFeaturePlot(Zari)
plot2 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE)
plot1 + plot2


plot1 <- VariableFeaturePlot(Zari, cols = c("red", "green"), pt.size = 1)
LabelPoints(plot = plot1, points = top10, colour="black", repel = TRUE, xnudge = 0, ynudge = 0)





#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Step 5. Data scaling*********************************************************


Zari <- ScaleData(Zari)

Zari <- ScaleData(Zari, vars.to.regress = c("nFeature_RNA", "percent.mt"))


#Alternative step 3-5: using SCTransform

Zari <- SCTransform(Zari, variable.features.n = 3000)


Zari <- SCTransform(Zari, vars.to.regress = c("nFeature_RNA", "percent.mt"),
                    variable.features.n = 3000)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Step 6. Perform Linear Dimensional Reduction***********************************

Zari <- RunPCA(Zari, npcs = 50)

ElbowPlot(Zari, ndims = ncol(Embeddings(Zari, "pca")))

Zari <- RunPCA(Zari, features = VariableFeatures(object = Zari))


# Examine and visualize PCA results a few different ways
print(Zari[["pca"]], dims = 1:5, nfeatures = 5)


VizDimLoadings(Zari, dims = 1:2, reduction = "pca")

DimPlot(Zari, reduction = "pca")


DimHeatmap(Zari, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(Zari, dims = 1:15, cells = 500, balanced = TRUE)



#Determining the “dimensionality” of the dataset



# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Zari <- JackStraw(Zari, num.replicate = 100)
Zari <- ScoreJackStraw(Zari, dims = 1:20)


JackStrawPlot(Zari, dims = 1:15)








#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Step 7. Non-linear dimension reduction for visualization***********************


Zari <- RunTSNE(Zari, dims = 1:20)
Zari <- RunUMAP(Zari, dims = 1:20)

#The results can be then visualized:
plot1 <- TSNEPlot(Zari)
plot2 <- UMAPPlot(Zari)
plot1 + plot2



plot1 <- FeaturePlot(Zari, c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","TFAP2A"),
                     ncol=3, reduction = "tsne")

plot2 <- FeaturePlot(Zari , c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","TFAP2A"),
                     ncol=3, reduction = "umap")

plot1 / plot2

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Step 8. Cluster the cells******************************************************

Zari <- FindNeighbors(Zari, dims = 1:20)

Zari <- FindClusters(Zari, resolution = 1)


plot1 <- DimPlot(Zari, reduction = "tsne", label = TRUE)
plot2 <- DimPlot(Zari, reduction = "umap", label = TRUE)
plot1 + plot2


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#https://rstudio-pubs-static.s3.amazonaws.com/682508_2e57986a58734fbd959aada7996e1b27.html
#Step 9. Annotate cell clusters*************************************************
#Assigning Cell Type Identity to Clusters

ct_markers <- c("MKI67","NES","DCX","FOXG1", # G2M, NPC, neuron, telencephalon
                "DLX2","DLX5","ISL1","SIX3","NKX2.1","SOX6","NR2F2", # ventral telencephalon related
                "EMX1","PAX6","GLI3","EOMES","NEUROD6", # dorsal telencephalon related
                "RSPO3","OTX2","LHX9","TFAP2A","RELN","HOXB2","HOXB5") # non-telencephalon related



DoHeatmap(Zari, features = ct_markers) + NoLegend()
"DLX2","DLX5","ISL1","SIX3","NKX2.1","SOX6","NR2F2", # ventral telencephalon related
"EMX1","PAX6","GLI3","EOMES","NEUROD6", # dorsal telencephalon related
"RSPO3","OTX2","LHX9","TFAP2A","RELN","HOXB2","HOXB5") # non-telencephalon related
DoHeatmap(Zari, features = ct_markers) + NoLegend()


cl_markers <- FindAllMarkers(Zari, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
library(dplyr)
cl_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)


cl_markers <- FindAllMarkers(Zari, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
library(dplyr)

cl_markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC) -> top4



cl_markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC) -> top4

head(top4)


library(presto)
cl_markers_presto <- wilcoxauc(seurat)
cl_markers_presto %>%
  filter(logFC > log(1.2) & pct_in > 20 & padj < 0.05) %>%
  group_by(group) %>%
  arrange(desc(logFC), .by_group=T) %>%
  top_n(n = 2, wt = logFC) %>%
  print(n = 40, width = Inf)






top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Zari, features = top10_cl_markers$gene) + NoLegend()



plot1 <- FeaturePlot(Zari, c("NEUROD2","NEUROD6"), ncol = 1)
plot2 <- VlnPlot(Zari, features = c("NEUROD2","NEUROD6"), pt.size = 0)
plot1 + plot2 + plot_layout(widths = c(1, 2))




library(voxhunt)
load_aba_data('ABA_data')
genes_use <- variable_genes('E13', 300)$gene
vox_map <- voxel_map(seurat, genes_use=genes_use)
plot_map(vox_map)




#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Step 9. Pseudo temporal cell ordering******************************************************



new_ident <- setNames(c("Dorsal telen. NPC",
                        "Midbrain-hindbrain boundary neuron",
                        "Dorsal telen. neuron",
                        "Dien. and midbrain excitatory neuron",
                        "MGE-like neuron","G2M dorsal telen. NPC",
                        "Dorsal telen. IP","Dien. and midbrain NPC",
                        "Dien. and midbrain IP and excitatory early neuron",
                        "G2M Dien. and midbrain NPC",
                        "G2M dorsal telen. NPC",
                        "Dien. and midbrain inhibitory neuron",
                        "Dien. and midbrain IP and early inhibitory neuron",
                        "Ventral telen. neuron",
                        "Unknown 1",
                        "Unknown 2"),
                      levels(Zari))
Zari <- RenameIdents(Zari, new_ident)


DimPlot(Zari, reduction = "umap", label = TRUE) + NoLegend()








seurat_dorsal <- subset(Zari, subset = RNA_snn_res.1 %in% c(0,2,5,6,10))
seurat_dorsal <- FindVariableFeatures(seurat_dorsal, nfeatures = 2000)



#different phase of cell cycle.

VariableFeatures(Zari) <- setdiff(VariableFeatures(Zari), unlist(cc.genes))

seurat_dorsal <- RunPCA(seurat_dorsal) %>% RunUMAP(dims = 1:20)
FeaturePlot(seurat_dorsal, c("MKI67","GLI3","EOMES","NEUROD6"), ncol = 4)




seurat_dorsal <- CellCycleScoring(seurat_dorsal,
                                  s.features = cc.genes$s.genes,
                                  g2m.features = cc.genes$g2m.genes,
                                  set.ident = TRUE)
seurat_dorsal <- ScaleData(seurat_dorsal, vars.to.regress = c("S.Score", "G2M.Score"))




seurat_dorsal <- RunPCA(seurat_dorsal) %>% RunUMAP(dims = 1:20)
FeaturePlot(seurat_dorsal, c("MKI67","GLI3","EOMES","NEUROD6"), ncol = 4)



library(destiny)
dm <- DiffusionMap(Embeddings(seurat_dorsal, "pca")[,1:20])
dpt <- DPT(dm)
seurat_dorsal$dpt <- rank(dpt$dpt)
FeaturePlot(seurat_dorsal, c("dpt","GLI3","EOMES","NEUROD6"), ncol=4)



fitted curve is usually a straightforward way.
library(ggplot2)
plot1 <- qplot(seurat_dorsal$dpt, as.numeric(seurat_dorsal@assays$RNA@data["GLI3",]),
               xlab="Dpt", ylab="Expression", main="GLI3") +
  geom_smooth(se = FALSE, method = "loess") + theme_bw()
plot2 <- qplot(seurat_dorsal$dpt, as.numeric(seurat_dorsal@assays$RNA@data["EOMES",]),
               xlab="Dpt", ylab="Expression", main="EOMES") +
  geom_smooth(se = FALSE, method = "loess") + theme_bw()
plot3 <- qplot(seurat_dorsal$dpt, as.numeric(seurat_dorsal@assays$RNA@data["NEUROD6",]),
               xlab="Dpt", ylab="Expression", main="NEUROD6") +
  geom_smooth(se = FALSE, method = "loess") + theme_bw()
plot1 + plot2 + plot3







#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Step 11. Save the result******************************************************




#saveRDS/readRDS to save/load every Seurat object separately,


saveRDS(seurat, file="DS1/seurat_obj_all.rds")
saveRDS(seurat_dorsal, file="DS1/seurat_obj_dorsal.rds")

seurat <- readRDS("DS1/seurat_obj_all.rds")
seurat_dorsal <- readRDS("DS1/seurat_obj_dorsal.rds")


#or use save/load to save multiple objects together

save(seurat, seurat_dorsal, file="DS1/seurat_objs.rdata")
load("DS1/seurat_objs.rdata")

















