# scMC: Integrating and comparing multiple single cell genomic datasets

## Capabilities
- scMC is an R toolkit for integrating and comparing multiple single cell genomic datasets from single cell RNA-seq and ATAC-seq experiments across different conditions, time points and tissues. 
- scMC exhibits superior performance in detecting context-shared and -specific biological signals, particularly noticeable for the datasets with imbalanced cell population compositions across interrelated biological conditions. 
- scMC learns a shared reduced dimensional embedding of cells that retains the biological variation while removing the technical variation. This shared embedding can enhance a variety of single cell analysis tasks, such as low-dimensional visualization, cell clustering and pseudotemporal trajectory inference. 

## Installation
To make it easy to run scMC in most common scRNA-seq data analysis pipelines, scMC is now implemented within Seurat V3 workflow. Please first **[install Seurat R pacakge (>= 3.2.1)](https://satijalab.org/seurat/install.html)** via ```install.packages('Seurat')```. For the standalone implementent of scMC and reproducing results from manuscript, please check out [previous release](http://doi.org/10.5281/zenodo.4395119).

scMC R package can then be easily installed from Github using devtools:  

```
devtools::install_github("amsszlh/scMC")
```
 
### Installation of other dependencies
- Install Leiden python pacakge for identifying cell clusters: ```pip install leidenalg```. Please check [here](https://github.com/vtraag/leidenalg) if you encounter any issue.


## Tutorials
The implementent of scMC is now seamlessly compatible with the workflow of Seurat V3 package. The runtime is also significantly reduced now. 

Please check out the full workflow

- [Demo of scMC workflow](https://htmlpreview.github.io/?https://github.com/amsszlh/scMC/blob/master/tutorial/demo_scMC_dermis.html)

We also wrote a Seurat Wrapper function `RunscMC` to run scMC directly on Seurat objects. You can run scMC within your Seurat V3 workflow. You'll only need to make two changes to your code.

- Run scMC with the `RunscMC()` function

- In downstream analyses, use the scMC embeddings instead of PCA.

For example, run scMC and then UMAP in two lines.

```
combined <- RunscMC(seuratObj.list)
combined <- RunUMAP(combined, reduction = "scMC")
```

For details, please check out

- [Demo of scMC Seurat Wrapper](https://htmlpreview.github.io/?https://github.com/amsszlh/scMC/blob/master/tutorial/demo_scMC_Seurat_Wrapper_dermis.html)


Here we also showcase scMC’s superior performance in detecting context-shared and -specific biological signals by applying it to a mouse skin scRNA-seq dataset and comparing it with other methods (Seurat, Harmony and LIGER)

- [Benchmarking of scMC against other methods](https://htmlpreview.github.io/?https://github.com/amsszlh/scMC/blob/master/tutorial/benchmark_against_other_methods.html)


## Help
If you have any problems, comments or suggestions, please contact us at Lihua Zhang (lihuaz1@uci.edu).


