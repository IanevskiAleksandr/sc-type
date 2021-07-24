
# ScType: Fully-automated cell identification by means of cell type-specific markers extracted from single-cell RNA-seq data

**Article**: [https://www.biorxiv.org/content/10.1101/812131v2.full.pdf]

<p style="text-align:justify;"> <b>ScType</b> a computational method for automated selection of marker genes based merely on scRNA-seq data. The open-source portal (<a href="//sctype.app">http://sctype.app</a>) provides an interactive web-implementation of the method.</p>

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
![alt text](https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/fig1.png)
<br><br>
Now, let's automatically assign cell types using ScType. For that, we first load 2 additional ScType functions:

```R
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

```

Next, let's prepare gene sets from the input cell marker file. By default, we use our in-built cell marker DB, however, feel free to use your own data.
Just prepare an input XLSX file in the same format as our DB file.
<br>In addition, provide a tissue type your data belongs to:

```R
# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx";
tissue = "Immune system" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

```

<br>
Finally, let's assign cell types to each cluster:

```R
cL_resutls = sctype_score(scRNAseqData = pbmc[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative, 
                      marker_sensitivity = gs_list$marker_sensitivity, verbose=!0)
cL_resutls %>% group_by(cluster) %>% top_n(n = 1)                
```

<br>
We can also overlay the identified cell types on UMAP plot:

```R
pbmc@meta.data$customclassif = ""
for(j in unique(cL_resutls$cluster)){
  cl_type = cL_resutls[cL_resutls$cluster==j,]; cl_type = cl_type[order(cl_type$scores, decreasing = T), ]
  if(cl_type$scores[1]>0){
    pbmc@meta.data$customclassif[pbmc@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  } else {
    pbmc@meta.data$customclassif[pbmc@meta.data$seurat_clusters == j] = "Unknown"
  }
}

DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')        

```
![alt text](https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/fig2.png)


<br><br>
sessionInfo();
[1] patchwork_1.1.1         geneSynonym_1.7.21.5.23 HGNChelper_0.8.1        SeuratObject_4.0.2      Seurat_4.0.3           
[6] dplyr_1.0.6            


<br><br>

For any questions please contact **Aleksandr Ianevski** (aleksandr.ianevski@helsinki.fi)

## Copyright and license

Code copyright 2021 ScType
