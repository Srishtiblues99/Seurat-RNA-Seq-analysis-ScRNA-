install.packages('Seurat')
library(Seurat)
remove.packages("Matrix")
install.packages("Matrix")
library(Matrix)
devtools::install_github("thomasp85/patchwork")

install.packages("devtools")
library(devtools)
library(patchwork)
data = Read10X(data.dir = "C:/Users/shristi/Documents/GSE192723_RAW/GSM5763644_filtered_feature_bc_matrix_A8/filtered_feature_bc_matrix_A8")
osm = CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)
osm
data[1:50, 1:10]

osm[["percent.mt"]] = PercentageFeatureSet(osm, pattern = "^MT-")
head(osm@meta.data)
VlnPlot(osm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 = FeatureScatter(osm, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(osm, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2
osm = subset(osm, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
osm
osm = NormalizeData(osm)
osm = FindVariableFeatures(osm, selection.method = "vst", nfeatures = 2000)
top10 = head(VariableFeatures(osm), 10)
top10
plot1 = VariableFeaturePlot(osm)
plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
all.genes = rownames(osm)
osm = ScaleData(osm, features = all.genes)
osm@assays$RNA@scale.data[1:50, 1:5]
osm = RunPCA(osm, features = VariableFeatures(object = osm))
DimHeatmap(osm, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(osm)
osm = FindNeighbors(osm, dims = 1:10)
osm = FindClusters(osm, resolution = 0.5)
head(osm@meta.data)
osm = RunUMAP(osm, dims = 1:10)
DimPlot(osm, reduction = "umap")
DimPlot(osm, reduction = "umap", label = T)
osm.markers = FindAllMarkers(osm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(osm.markers)
if (packageVersion("devtools") < 1.6) {
  install.packages("devtools")
}
devtools::install_github("hadley/lazyeval")
devtools::install_github("hadley/dplyr")
library(dplyr)
a1 = osm.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
a1
genes = a1 %>% pull(gene)
genes
FeaturePlot(osm, features = genes[1:2])
FeaturePlot(osm, features = genes[1:2], cols = c("blue", "pink"))
