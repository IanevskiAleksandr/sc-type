
# ScType: Fully-automated cell identification by means of cell type-specific markers extracted from single-cell RNA-seq data

**Article**: [https://www.biorxiv.org/content/10.1101/812131v2.full.pdf]

<p style="text-align:justify;"> <b>ScType</b> a computational method for automated selection of marker genes based merely on scRNA-seq data. The open-source portal (http://sctype.fimm.fi) provides an interactive web-implementation of the method.</p>

##
<br><br>

![alt text](https://github.com/IanevskiAleksandr/sc-type/blob/master/ScTypePlan.png)

<br><br>


### Cell type annotation example 

First let's load a PBMC 3k example dataset (see Seurat tutorial for more details on how to load the dataset using Seurat, https://satijalab.org/seurat/articles/pbmc3k_tutorial.html). The raw data can be found <a href='https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz' download>here</a>.

```R
# load libraries
lapply(c("dplyr","Seurat","patchwork","HGNChelper","geneSynonym"), library, character.only = T)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
```

Next, let's normalize and cluster the data.

```R
# normalize data
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# scale and run PCA
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Check number of PC components
ElbowPlot(pbmc)

# cluster and visualize
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
```


Now, let's automatically assign cell types using ScType

```R
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score.R")

```



<br><br>

<br><br>

For any questions please contact **Aleksandr Ianevski** [@IanevskiAleksandr](aleksandr.ianevski@helsinki.fi)

## Copyright and license

Code copyright 2021 ScType
