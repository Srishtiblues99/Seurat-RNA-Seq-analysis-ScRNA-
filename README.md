# Seurat-RNA-Seq-analysis-ScRNA data is generated using the 10X genomics platform. The cellranger pipeline outputs two types of feature-barcode matrices. Each element of the matrix is the number of UMIs (Short random molecular tags which are added to DNA fragments in library preparation process before PCR amplification. UMIs are used to identify input DNA molecules) associated with a feature (row) and a barcode (column).
Filtered feature-barcode matrix Contains only detected cellular barcodes. For Targeted Gene Expression samples, non-targeted genes are removed from the filtered matrix.
Each matrix is stored in the Market Exchange Format (MEX) for sparse matrices. It also contains gzipped TSV files with feature and barcode sequences corresponding to row and column indices respectively.
Start by reading in the data. The Read10X() function reads in the output of the cellranger pipeline from 10X, returning a unique molecular identified (UMI) count matrix. The values in this matrix represent the number of molecules for each feature (i.e. gene; row) that are detected in each cell (column).
Next use the count matrix to create a Seurat object. The object serves as a container that contains both data (like the count matrix) and analysis (like PCA, or clustering results) for a single-cell dataset.
scRNA-seq raw data includes reads with cell barcode (unique nucleic acid sequences, termed barcodes, which are used to label individual cells, so that they can be tracked through space and time in scRNA-seq) and UMIs. Before alignment of reads to genome, reads can be grouped using cell barcodes and the frequency of each read is estimated per cell per gene using UMIs. After alignment and frequency calculations, we have a gene expression table containing cells represented in the columns and genes represented in the rows. This is what is analyzed in scRNA-seq comlutational workflow.
Here is an image of ScRNA seq workflow- ![image](https://user-images.githubusercontent.com/114133070/198524323-3d33551d-b628-439f-a494-0f3a89c9587e.png)

# install.packages('Seurat')
library(Seurat)
#remove.packages("Matrix")
#install.packages("Matrix")
library(Matrix)
#devtools::install_github("thomasp85/patchwork")

#install.packages("devtools")
library(devtools)
library(patchwork)

# Load the PBMC dataset
data = Read10X(data.dir = "C:/Users/shristi/Documents/GSE192723_RAW/GSM5763644_filtered_feature_bc_matrix_A8/filtered_feature_bc_matrix_A8")

# Initialize the Seurat object with the raw (non-normalized data)
osm = CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)
osm
#Here the parameters used are:
counts: Un normalized data such as raw counts or TPMs
min.cells: Include features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.
min.features: Include cells where at least this many features are detected.
The . values in the matrix represent 0s (no molecules detected). Since most values in an scRNA-seq matrix are 0, Seurat uses a sparse-matrix representation whenever possible. This results in significant memory and speed savings for Drop-seq/inDrop/10x data.

data[1:50, 1:10]

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#Mitochondrial QC metrics is calculated with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features
#Set of all genes is used starting with MT- as a set of mitochondrial genes
osm[["percent.mt"]] = PercentageFeatureSet(osm, pattern = "^MT-")
head(osm@meta.data)

# Visualize QC metrics as a violin plot
VlnPlot(osm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
#for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 = FeatureScatter(osm, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(osm, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2
osm = subset(osm, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
osm
osm = NormalizeData(osm)
osm = FindVariableFeatures(osm, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 = head(VariableFeatures(osm), 10)
top10

# plot variable features with and without labels
plot1 = VariableFeaturePlot(osm)
plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# Now we calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others).
#Genes in downstream analysis helps to highlight biological signal in single-cell datasets. Seurat return 2,000 features per dataset.
#Apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function:
a.Shifts the expression of each gene, so that the mean expression across cells is 0
b.Scales the expression of each gene, so that the variance across cells is 1
# For highly-expressed genes do not dominate, downstream analyses is done in equal weight
all.genes = rownames(osm)
osm = ScaleData(osm, features = all.genes)
osm@assays$RNA@scale.data[1:50, 1:5]

# PCA on the scaled data is performed
osm = RunPCA(osm, features = VariableFeatures(object = osm))

#DimHeatmap() allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses
DimHeatmap(osm, dims = 1:15, cells = 500, balanced = TRUE)

#‘Elbow plot’: a ranking of principle components is based on the percentage of variance explained by ElbowPlot() function
ElbowPlot(osm)

# Computing nearest neighbour
osm = FindNeighbors(osm, dims = 1:10)

# Computing SNN
osm = FindClusters(osm, resolution = 0.5)
head(osm@meta.data)
osm = RunUMAP(osm, dims = 1:10)

#SNE and UMAP are non-linear dimensional reduction techniquest offered by Seurat, to visualize and explore datasets. The goal of these algorithms is to learn the underlying manifold of the data in order to place similar cells together in low-dimensional space. 
DimPlot(osm, reduction = "umap")
DimPlot(osm, reduction = "umap", label = T)
osm.markers = FindAllMarkers(osm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(osm.markers)
#if (packageVersion("devtools") < 1.6) {
  #install.packages("devtools")
}
#devtools::install_github("hadley/lazyeval")
#devtools::install_github("hadley/dplyr")
library(dplyr)
a1 = osm.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
a1
genes = a1 %>% pull(gene)
genes
FeaturePlot(osm, features = genes[1:2])
FeaturePlot(osm, features = genes[1:2], cols = c("blue", "pink"))
