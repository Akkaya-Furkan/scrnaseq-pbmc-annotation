"""Transparent QC, marker scoring and comparison utilities."""
import numpy as np
import pandas as pd
from scipy import sparse

MARKERS = {
    'T cells': ['CD3D', 'CD3E', 'TRAC'],
    'B cells': ['MS4A1', 'CD79A', 'CD79B'],
    'NK cells': ['NKG7', 'GNLY', 'KLRD1', 'PRF1'],
    'Monocytes': ['LST1', 'S100A8', 'S100A9', 'FCN1'],
    'Dendritic cells': ['FCER1A', 'CD1C', 'CLEC10A'],
    'Platelets': ['PPBP', 'PF4'],
}


def mad_outliers(values, n_mads=5):
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).all():
        raise ValueError('QC metrics must be finite')
    median = np.median(values)
    mad = np.median(np.abs(values - median))
    # Constant metrics carry no evidence for removing a cell.
    return np.zeros(values.size, dtype=bool) if mad == 0 else abs(values - median) > n_mads * mad


def qc_flags(obs):
    """Five MADs on log depth/detection plus explicit permissive hard limits."""
    return pd.DataFrame({
        'low_genes': obs.n_genes_by_counts < 200,
        'high_mt': obs.pct_counts_mt > 10,
        'depth_outlier': mad_outliers(np.log1p(obs.total_counts)),
        'gene_outlier': mad_outliers(np.log1p(obs.n_genes_by_counts)),
    }, index=obs.index)


def marker_annotation(adata, cluster_key='leiden', margin_threshold=0.15, support_threshold=0.70):
    """Mean per-gene cell z scores, then equal-weight panel/cluster averages.

    Only marker columns are densified; no CellTypist predictions enter scoring.
    Low top scores, small margins or weak within-cluster support yield uncertain labels.
    """
    available = {label: [g for g in genes if g in adata.var_names]
                 for label, genes in MARKERS.items()}
    if any(len(genes) < 2 for genes in available.values()):
        raise ValueError('At least two detected marker genes per panel are required')
    genes = list(dict.fromkeys(g for panel in available.values() for g in panel))
    matrix = adata[:, genes].X
    matrix = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)
    std = matrix.std(axis=0)
    z = (matrix - matrix.mean(axis=0)) / np.where(std == 0, 1, std)
    z = np.clip(z, -5, 5)
    scores = pd.DataFrame({label: z[:, [genes.index(g) for g in panel]].mean(axis=1)
                           for label, panel in available.items()}, index=adata.obs_names)
    means = scores.groupby(adata.obs[cluster_key], observed=True).mean()
    rows = []
    for cluster, row in means.iterrows():
        ordered = row.sort_values(ascending=False)
        best = ordered.index[0]
        margin = float(ordered.iloc[0] - ordered.iloc[1])
        members = adata.obs[cluster_key].astype(str) == str(cluster)
        support = float((scores.loc[members].idxmax(axis=1) == best).mean())
        uncertain = (margin < margin_threshold or ordered.iloc[0] <= 0
                     or support < support_threshold)
        rows.append({'cluster': str(cluster), 'marker_top_type': best,
                     'marker_type': f'Uncertain ({best})' if uncertain else best,
                     'marker_margin': margin, 'marker_support': support, 'uncertain': uncertain})
    return pd.DataFrame(rows).set_index('cluster'), means, scores


# Exact label dictionary: unmapped labels remain Unknown, never guessed by substring.
CELLTYPIST_BROAD = {
    'B cells': 'B cells', 'Naive B cells': 'B cells', 'Memory B cells': 'B cells',
    'Age-associated B cells': 'B cells', 'Plasma cells': 'B cells',
    'CD16+ NK cells': 'NK cells', 'CD16- NK cells': 'NK cells', 'NK cells': 'NK cells',
    'Monocytes': 'Monocytes', 'Classical monocytes': 'Monocytes',
    'Non-classical monocytes': 'Monocytes', 'Intermediate monocytes': 'Monocytes',
    'DC': 'Dendritic cells', 'DC1': 'Dendritic cells', 'DC2': 'Dendritic cells',
    'pDC': 'Dendritic cells', 'DC3': 'Dendritic cells',
    'Megakaryocytes/platelets': 'Platelets',
    'T cells': 'T cells', 'CD4+ T cells': 'T cells', 'CD8+ T cells': 'T cells',
    'Tcm/Naive helper T cells': 'T cells', 'Tem/Trm cytotoxic T cells': 'T cells',
    'Tem/Effector helper T cells': 'T cells', 'Tcm/Naive cytotoxic T cells': 'T cells',
    'Treg': 'T cells', 'MAIT cells': 'T cells', 'gdT cells': 'T cells',
    'Tem/Temra cytotoxic T cells': 'T cells',
    'Trm cytotoxic T cells': 'T cells', 'Trm helper T cells': 'T cells',
    'Proliferating T cells': 'T cells',
    'Regulatory T cells': 'T cells', 'Treg(diff)': 'T cells',
    'Follicular helper T cells': 'T cells', 'Type 1 helper T cells': 'T cells',
    'Type 17 helper T cells': 'T cells', 'gamma-delta T cells': 'T cells',
    'CRTAM+ gamma-delta T cells': 'T cells', 'Cycling T cells': 'T cells',
    'Cycling gamma-delta T cells': 'T cells', 'Memory CD4+ cytotoxic T cells': 'T cells',
    'Tem/Effector helper T cells PD1+': 'T cells', 'CD8a/a': 'T cells',
    'CD8a/b(entry)': 'T cells', 'T(agonist)': 'T cells',
    'Cycling B cells': 'B cells', 'Follicular B cells': 'B cells',
    'Germinal center B cells': 'B cells', 'Transitional B cells': 'B cells',
    'Plasmablasts': 'B cells', 'Proliferative germinal center B cells': 'B cells',
    'Cycling NK cells': 'NK cells', 'Transitional NK': 'NK cells',
    'Cycling DCs': 'Dendritic cells', 'Migratory DCs': 'Dendritic cells',
    'Transitional DC': 'Dendritic cells', 'Cycling monocytes': 'Monocytes',
}


def broad_labels(labels):
    return labels.map(CELLTYPIST_BROAD).fillna('Unknown')
