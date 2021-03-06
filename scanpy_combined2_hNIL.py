# -*- coding: utf-8 -*-
"""scanpy_pool1_hNIL.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1u3uymUGrV0wjqOx3QrLtxuScjxjESo1C

# Set up the working directory and enviroment
"""

#from google.colab import drive
#drive.mount('/content/drive/')

import os
os.chdir("/content/drive/Shared drives/CARD/projects/iNDI/line_prioritization/projects_lirong/combined2_hNIL")

"""## Assign folders for storing input data, object, and output
1. 'data' folder: storing .mtx or h5 from 10xGenomic sequencing 
2. 'interim' folder: storing h5ad object after processing and its output files, e.g., top_markers.csv
3. 'figures' folder: automatically generated
"""

# Commented out IPython magic to ensure Python compatibility.
# %ls

#! mkdir scanpy_out
results_file = 'scanpy_out/combined2_hNIL.h5ad'
figure_file = '_combined2_hNIL.pdf'



"""## Installed the required packages"""

# required for umap clustering
# pip3 install leidenalg 
# pip3 install scanpy
# ! pip install python-igraph
# ! pip install louvain
import h5py
import numpy as np
import pandas as pd
import scanpy as sc

"""## Set up the displaying and print parameters"""

sc.settings.verbosity = 3 
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80)

"""# Data preprocessing
1. Use gene symbols for the variable names (variables-axis index)
3. Write a cache file for faster subsequent reading
4. De-duplicates

## Read 10xGenomics sc-RNA sequencing data
"""

adata = sc.read_10x_h5("filtered_feature_bc_matrix.h5")

adata.var_names_make_unique()

#adata

"""## Basic filtering
1. cell based filtering: remove cells with less than 200 genes
2. gene based filtering: remove genes expressing in less than 3 cells
"""

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

"""## Calculate the percentage of mitochondrial genes"""

mito_genes = adata.var_names.str.startswith('MT-')

adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1

adata.obs['n_counts'] = adata.X.sum(axis=1).A1

"""## Check sequencing quality
1. choose the threthold of gene numbers to remove, e.g., n_genes = 4500
2. choose the threthold of mitochondial genes to remove, e.g., percent_mito = 0.15
"""

sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, save=figure_file)

sc.pl.scatter(adata, x='n_counts', y='percent_mito', save="scatter1"+figure_file)
sc.pl.scatter(adata, x='n_counts', y='n_genes', save="scatter2"+figure_file)

adata = adata[adata.obs.n_genes < 4500, :]
adata = adata[adata.obs.percent_mito < 0.15, :]
#adata

"""## Scale and logarithmize the data
option: store the unnormalized data in .raw
"""

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw=adata

"""## Choosing highly-variable genes for further analysis
Subset is optional. If subset using adata.var.highly_variable, the adata will only contain these genes.
"""

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

sc.pl.highly_variable_genes(adata, save=figure_file)

adata = adata[:, adata.var.highly_variable]
#adata

"""## Further scale on cofounder "n_counts" and "percent_mito"
1. Regression out n_counts and percent_mito effect and scale again on adata which only contains highly-variable genes now.
2. Clip values exceeding standard deviation 10 (max_value=10)
"""

sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
sc.pp.scale(adata, max_value=10)

"""# Principal component analysis"""

sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True,save=figure_file)

adata.write_h5ad("scanpy_out/combined2_hNIL_before_neighbors.h5ad")

"""### Computing, embedding, and clustering the neighborhood graph

1. Computing the neighborhood graph of cells using the PCA representation of the data matrix.
2. Embedding the graph in 2 dimensions using UMAP.
3. Clustering the neighborhood graph using Leiden graph-clustering method
"""

# defaults are: n_neighbors= 20, n_pcs=50
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# default resolution=1.0
sc.tl.leiden(adata)

#adata

#sc.tl.leiden(adata, resolution=0.4, key_added = "leiden_0.4")
#sc.tl.leiden(adata, resolution=0.6, key_added = "leiden_0.6")

"""## Visualize the clusters"""

sc.pl.umap(adata, color=['leiden'], save=figure_file)

#sc.pl.umap(adata, color=['leiden_0.6'])

sc.pl.umap(adata, color=['MAPT','TMEM106B'], save=figure_file)

"""## Save the processed adata object saved as .h5ad
The object can be load using adata = sc.read(results_file)
"""

adata.write(results_file)

"""# Finding marker genes
one vs rest comparison using Mann-Whitney-U-test (recommend)
"""

# reduce the verbosity from 3 to 2 in the setting of logging output
sc.settings.verbosity = 2

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, save=figure_file)

"""## Visualize the differential expression of marker genes across clusters
1. Differential expression a set of marker genes of a specific cluster vs the rest.
2. Differential expression of a single gene across all clusters
"""

#sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)

#sc.pl.violin(adata, ['TMEM106B', 'MAPT', "DLK1"], groupby='leiden')

"""## Annotate cell types for all the clusters 
Based on domain knowledge and information of marker genes

Wait for disscussion for this case
"""

# assign cell type names to the clusters


# adata.rename_categories('leiden', new_cluster_names)

# it automatically generates a folder of 'figures' and save the figure inside 
# sc.pl.umap(adata, color='leiden', legend_loc='on data', title='', frameon=False, save='.pdf')

"""# Export a list of marker genes"""

result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names

top_marker_genes = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(50)

top_marker_genes.to_csv('scanpy_out/top_markers_pool2_hNIL.csv')

#top_marker_genes.head()

# get a list of top20 marker genes of cluster 1
#marker_genes_cluster1 = top_marker_genes["1_n"][:20].values

"""# Visulize gene differential expression 
maker genes or other gene of interest
"""

#ax = sc.pl.dotplot(adata, marker_genes_cluster1, groupby='leiden')

#ax = sc.pl.stacked_violin(adata, marker_genes_cluster1, groupby='leiden', rotation=90)

"""# Save the work
option: export a set of csv using 
adata.write_csvs(results_file[:-5])
"""

adata.write(results_file)

# Commented out IPython magic to ensure Python compatibility.
# %ls