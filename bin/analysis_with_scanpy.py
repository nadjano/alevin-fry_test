import scanpy as sc
import pandas as pd
import argparse


parser = argparse.ArgumentParser(description='Convert Alevin outputs to 10X .mtx.')
parser.add_argument('mtx_dir1', help = 'Alevin output directory')
parser.add_argument('mtx_dir2', help = 'Alevin output directory')
parser.add_argument('out_dir', help = 'Output directory for converted results')
# parser.add_argument('--cell_prefix', dest='cell_prefix', default='', help = 'Prefix to apply to cell barcodes')
args = parser.parse_args() 

mtx_dir=args.mtx_dir
out_dir=args.out_dir

adata = sc.read_10x_mtx(mtx_dir)

sc.pp.filter_cells(adata, min_genes=750)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save= 'type.png')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)
sc.pl.umap(adata)