library(dplyr)
library(Seurat)
library(patchwork)
pbmc <- as.matrix(read.table(file = "exprMatrix.tsv" ,row.names = 1 , header = TRUE))
pbmc <- CreateSeuratObject(counts = pbmc , project = "RA")

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:20)
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 2)

#UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap",label=TRUE)

#Visualize canonical markers
library(ggplot2)
features <- c("CD3D","CD8A","CD4","GNLY","NKG7","IL7R","CCR7","CD14","LYZ","FCGR3A","MS4A7","MS4A1","CD27","CD38","TNFRSF17","COL1A1","COL1A2","DCN")
# Function to create a violin plot
create_vln_plot <- function(feature) {
  VlnPlot(pbmc,features = c(feature), pt.size = 0, add.noise = FALSE) +
    theme(legend.position = 'none', 
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          text = element_text(family = "Arial"),
          plot.title = element_text(size = 15))
}
create_feature_plot <- function(feature) {
  FeaturePlot(pbmc,
              features = c(feature), 
              pt.size = 0, 
              label = TRUE) +
    theme(legend.position = "none",  # Removes color bar
          axis.title.x = element_blank(),  # Removes x axis label
          axis.title.y = element_blank(),  # Removes y axis label
          axis.text.x = element_blank(),   # Optional: removes x axis text
          axis.text.y = element_blank())   # Optional: removes y axis text
}
# Create combined violin plot
all_vln_plots <- lapply(features, create_vln_plot)
combined_vln <- wrap_plots(all_vln_plots, ncol = 3)  # Adjust ncol as needed
png("all_vln.png", res = 300, height = 3000, width = 3000)  # Adjust dimensions as needed
print(combined_vln)
dev.off()
# Create combined feature plot
all_feature_plots <- lapply(features, create_feature_plot)
combined_feature <- wrap_plots(all_feature_plots, ncol = 5)  # Adjust ncol as needed
png("all_feature.png", res = 300, height = 3000, width = 3000)  # Adjust dimensions as needed
print(combined_feature)
dev.off()


##Annotate cell type based on canonical marker expression
pbmc <- RenameIdents(pbmc, `0` = "CD4+T", `1` = "CD8+T", `2` = "Fibroblast",`3` = "Fibroblast", `4` = "Plasmablast", `5` = "B", `6` = "B",`7`="CD4+T",`8`="Monocyte",`9`="Fibroblast",`10`="Fibroblast",`11`="Fibroblast",`12`="Fibroblast",`13`="Fibroblast",`14`="Monocyte",`15`="CD4+T",`16`="Monocyte",`17`="B",`18`="CD8+T",`19`="CD8+T",`20`="CD4+T",`21`="Monocyte",`22`="B")
metadata<-read.table(file="celseq_meta.tsv.725591",sep="\t",header=1,row=1)
pbmc<-AddMetaData(pbmc,metadata)

#Visualize HGF and MET expression
orders=c("CD4+T", "CD8+T", "B", "Plasmablast", "Monocyte", "Fibroblast")
pbmc_RA=subset(pbmc, subset = disease != "OA")
png("~/feature/HGF.png", res = 300, height = 1500, width = 1500)
p=FeaturePlot(pbmc_RA,features=c("HGF"),label=TRUE)
print(p)
dev.off()
png("~/vln/HGF.png", res = 300, height = 1000, width = 1400)
p = VlnPlot(pbmc_RA,features = c("HGF")) +NoLegend() +scale_x_discrete(limits = orders)  +xlab("")
print(p)
dev.off()
png("~/feature/MET.png", res = 300, height = 1500, width = 1500)
p=FeaturePlot(pbmc_RA,features=c("MET"),label=TRUE)
print(p)
dev.off()
png("~/vln/MET.png", res = 300, height = 1000, width = 1500)
p = VlnPlot(pbmc_RA,features = c("MET")) +NoLegend() +scale_x_discrete(limits = orders)+xlab("")
print(p)
dev.off()
