"""Run with: python -m pbmc_annotation.pipeline"""
import argparse
import importlib.metadata
import json
import shutil
import time
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models
from .core import MARKERS, CELLTYPIST_BROAD, qc_flags, marker_annotation, broad_labels

SEED = 42


def run(data_dir, output):
    started = time.monotonic()
    output.mkdir(parents=True, exist_ok=True)
    # An interrupted rerun must not leave an old success summary for mixed artifacts.
    (output / 'summary.json').unlink(missing_ok=True)
    sc.settings.n_jobs = 2
    sc.settings.set_figure_params(dpi=110, dpi_save=160, frameon=False)
    adata = sc.read_h5ad(data_dir / 'pbmc3k_raw.h5ad')
    adata.var_names_make_unique()
    input_cells, input_genes = adata.shape
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True)
    flags = qc_flags(adata.obs)
    audit = adata.obs.join(flags)
    audit['qc_pass'] = ~flags.any(axis=1)
    passing = adata[audit.qc_pass].copy()
    sc.pp.filter_genes(passing, min_cells=3)
    sc.pp.scrublet(passing, expected_doublet_rate=0.03, random_state=SEED,
                   n_prin_comps=30, verbose=True)
    audit = audit.join(passing.obs[['doublet_score', 'predicted_doublet']])
    audit['retained'] = audit.qc_pass & ~audit.predicted_doublet.fillna(False).astype(bool)
    audit.to_csv(output / 'qc_cells.csv', index_label='barcode')
    sc.pl.scrublet_score_distribution(passing, show=False)
    plt.savefig(output / 'doublet_scores.png', bbox_inches='tight')
    plt.close('all')
    doublets = int(passing.obs.predicted_doublet.sum())
    doublet_threshold = float(passing.uns['scrublet']['threshold'])
    adata = passing[~passing.obs.predicted_doublet].copy()
    del passing
    sc.pp.filter_genes(adata, min_cells=3)
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
    embedding = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(embedding, max_value=10)
    sc.tl.pca(embedding, n_comps=40, svd_solver='arpack', random_state=SEED)
    adata.obsm['X_pca'] = embedding.obsm['X_pca'].copy()
    adata.uns['pca'] = embedding.uns['pca'].copy()
    del embedding
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=SEED)
    sc.tl.leiden(adata, resolution=0.5, random_state=SEED, flavor='igraph', n_iterations=2, directed=False)
    sc.tl.umap(adata, random_state=SEED)
    annotations, panel_means, scores = marker_annotation(adata)
    for column in ['marker_type', 'marker_top_type', 'marker_margin', 'marker_support', 'uncertain']:
        adata.obs[column] = adata.obs.leiden.astype(str).map(annotations[column])
    panel_means.to_csv(output / 'cluster_marker_scores.csv')
    scores.to_csv(output / 'cell_marker_scores.csv', index_label='barcode')
    model_path = models.get_model_path('Immune_All_Low.pkl')
    # No cluster voting: preserve a per-cell comparison independent of Leiden smoothing.
    predictions = celltypist.annotate(adata, model=model_path, majority_voting=False)
    adata.obs['celltypist_label'] = predictions.predicted_labels.predicted_labels.astype(str)
    adata.obs['celltypist_broad'] = broad_labels(adata.obs.celltypist_label)
    adata.obs['celltypist_confidence'] = predictions.probability_matrix.max(axis=1)
    adata.obs['broad_agreement'] = adata.obs.marker_type == adata.obs.celltypist_broad
    pd.Series(CELLTYPIST_BROAD, name='broad_type').to_csv(output / 'celltypist_label_mapping.csv', index_label='model_label')
    pd.crosstab(adata.obs.marker_type, adata.obs.celltypist_broad).to_csv(output / 'annotation_crosstab.csv')
    pd.crosstab(adata.obs.leiden, adata.obs.celltypist_label).to_csv(output / 'cluster_celltypist_labels.csv')
    summaries = adata.obs.groupby('leiden', observed=True).agg(
        cells=('marker_type', 'size'), marker_type=('marker_type', 'first'),
        marker_margin=('marker_margin', 'first'), marker_support=('marker_support', 'first'), uncertain=('uncertain', 'first'),
        celltypist_mode=('celltypist_broad', lambda x: x.value_counts().index[0]),
        broad_agreement=('broad_agreement', 'mean'),
        mean_model_confidence=('celltypist_confidence', 'mean'))
    summaries.to_csv(output / 'cluster_annotations.csv')
    adata.obs.to_csv(output / 'cell_annotations.csv', index_label='barcode')
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', use_raw=False, pts=True)
    sc.get.rank_genes_groups_df(adata, group=None).to_csv(output / 'ranked_markers.csv', index=False)
    sc.pl.umap(adata, color=['marker_type', 'leiden', 'celltypist_broad'], wspace=0.65, show=False)
    plt.savefig(output / 'umap.png', bbox_inches='tight')
    plt.close('all')
    panels = {k: [g for g in v if g in adata.var_names] for k, v in MARKERS.items()}
    plot = sc.pl.dotplot(adata, panels, groupby='leiden', standard_scale='var',
                        use_raw=False, return_fig=True, show=False)
    plot.savefig(output / 'marker_dotplot.png', dpi=160, bbox_inches='tight')
    plt.close('all')
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
    for ax, metric in zip(axes, ['total_counts', 'n_genes_by_counts', 'pct_counts_mt']):
        edges = np.histogram_bin_edges(audit[metric], bins=50)
        ax.hist(audit[metric], bins=edges, alpha=.6, label='Input')
        ax.hist(adata.obs[metric], bins=edges, alpha=.6, label='Retained')
        ax.set_xlabel(metric)
        ax.set_ylabel('Cells')
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(output / 'qc.png', dpi=160)
    plt.close(fig)
    params = {'seed': SEED, 'min_genes': 200, 'max_mt_percent': 10, 'mad_cutoff': 5,
              'min_cells_per_gene': 3, 'expected_doublet_rate': 0.03,
              'hvg': 2000, 'pca_components': 40, 'neighbor_pcs': 30,
              'neighbors': 15, 'leiden_resolution': 0.5, 'marker_margin_cutoff': 0.15, 'marker_support_cutoff': 0.70}
    summary = {'input_cells': input_cells, 'input_genes': input_genes,
               'qc_pass_cells': int(audit.qc_pass.sum()), 'qc_removed_cells': int((~audit.qc_pass).sum()),
               'predicted_doublets_removed': doublets, 'scrublet_threshold': doublet_threshold,
               'retained_cells': adata.n_obs, 'retained_genes': adata.n_vars,
               'clusters': int(adata.obs.leiden.nunique()),
               'uncertain_cells': int(adata.obs.uncertain.sum()),
               'uncertain_clusters': int(summaries.uncertain.sum()),
               'unknown_celltypist_cells': int((adata.obs.celltypist_broad == 'Unknown').sum()),
               'candidate_agreement_cells': int((adata.obs.marker_top_type.astype(str) == adata.obs.celltypist_broad.astype(str)).sum()),
               'broad_agreement_cells': int(adata.obs.broad_agreement.sum()),
               'broad_agreement_fraction': float(adata.obs.broad_agreement.mean()),
               'cell_types': adata.obs.marker_type.value_counts().to_dict(),
               'qc_flag_counts_overlapping': flags.sum().astype(int).to_dict(),
               'parameters': params,
               'versions': {p: importlib.metadata.version(p) for p in
                            ['scanpy', 'anndata', 'celltypist', 'numpy', 'pandas', 'scikit-learn', 'scrublet']}}
    adata.uns['pipeline_parameters'] = params
    adata.write_h5ad(output / 'pbmc3k_annotated.h5ad', compression='gzip')
    # Reopen the persisted artifact before reporting success.
    saved = sc.read_h5ad(output / 'pbmc3k_annotated.h5ad')
    assert saved.shape == adata.shape
    assert saved.obs_names.equals(adata.obs_names)
    assert np.isfinite(saved.obsm['X_umap']).all()
    shutil.copyfile(data_dir / 'manifest.json', output / 'data_manifest.json')
    summary['elapsed_seconds'] = round(time.monotonic() - started, 3)
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=Path('data/raw'))
    parser.add_argument('--output', type=Path, default=Path('results'))
    args = parser.parse_args()
    run(args.data_dir, args.output)


if __name__ == '__main__':
    main()
