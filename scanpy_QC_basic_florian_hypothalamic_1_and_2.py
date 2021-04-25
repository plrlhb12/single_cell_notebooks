# -*- coding: utf-8 -*-

#make sure the current directory is specific to an individual sample and haveing data directly in this directory
#make sure there is a subfolder of scanpy_out under the currrent directory
import os
os.chdir("/content/drive/Shared drives/CARD/projects/iNDI/line_prioritization/projects_lirong/Florian_data/data/filtered_cellranger_matrix/hypothalamic_sample_2")

import h5py
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3 
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80)

sample_name = 'hypothalamic_sample_2'
results_file = 'scanpy_out/'+sample_name+'.h5ad'
results_file2 = 'scanpy_out/'+sample_name+'_unnormalized.h5ad'
figure_file = '_'+sample_name+'.pdf'

## Read 10xGenomics sc-RNA sequencing data using mtx

#adata = sc.read_10x_h5("filtered_feature_bc_matrix.h5")
adata = sc.read_10x_mtx("../"+sample_name, var_names='gene_symbols')
adata.var_names_make_unique()

#check the most abundantly expressed genes
sc.pl.highest_expr_genes(adata, n_top=20, )

## Basic filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_genes(adata, min_counts=1)
## Calculate the percentage of mitochondrial genes

mito_genes = adata.var_names.str.startswith('MT-')

adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1

adata.obs['n_counts'] = adata.X.sum(axis=1).A1

## Check sequencing quality

sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, save='violin_'+figure_file)

sc.pl.scatter(adata, x='n_counts', y='percent_mito', save="scatter1_"+figure_file)
sc.pl.scatter(adata, x='n_counts', y='n_genes', save="scatter2_"+figure_file)

adata = adata[adata.obs.n_genes < 10000, :]
adata = adata[adata.obs.percent_mito < 0.15, :]
adata.write(results_file2)

#Scale and logarithmize the data
#note that adata.raw has been processed by normalization and log tranformation, but not go through cofounder correction yet
#adata.raw can be used for differentail expression analaysis

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw=adata

# Choosing highly-variable genes (HVG) for further analysis, here we do not subset the adata using HVG

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata, save=figure_file)
#if not subseeting, the regress_out function will generate error
adata = adata[:, adata.var.highly_variable]

# Further scale on cofounder "n_counts" and "percent_mito"

sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
sc.pp.scale(adata, max_value=10)

# Principal component analysis

sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True,save=figure_file)
sc.pl.pca(adata, save=figure_file)

# Computing, embedding, and clustering the neighborhood graph
# defaults are: n_neighbors= 20, n_pcs=50

sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
sc.tl.leiden(adata, resolution=0.2, key_added = "leiden_0.2")
sc.tl.leiden(adata, resolution=0.4, key_added = "leiden_0.4")

## Visualize the clusters

sc.pl.umap(adata, color=['leiden'], save="1"+figure_file)
sc.pl.umap(adata, color=['leiden_0.2'], save="0.2"+figure_file)
sc.pl.umap(adata, color=['leiden_0.4'], save="0.4"+figure_file)
sc.pl.umap(adata, color=['MAPT','TMEM106B'], save="gene"+figure_file)

# Finding marker genes using one vs rest comparison using Mann-Whitney-U-test (recommend)
#sc.settings.verbosity = 2

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, save=figure_file)

adata.write(results_file)
# Export a list of marker genes

result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names

top_marker_genes = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(50)

top_marker_genes.to_csv('scanpy_out/top_markers_'+sample_name+'.csv')

