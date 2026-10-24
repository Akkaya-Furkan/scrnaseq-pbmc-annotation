import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix
from pbmc_annotation.core import MARKERS, broad_labels, mad_outliers, marker_annotation, qc_flags


def test_mad_constant_and_extreme():
    assert not mad_outliers([5, 5, 5]).any()
    assert mad_outliers([1, 2, 3, 4, 100]).tolist() == [False, False, False, False, True]
    with pytest.raises(ValueError, match='finite'):
        mad_outliers([1, np.nan])


def test_qc_boundaries_and_overlapping_flags():
    obs = pd.DataFrame({'n_genes_by_counts': [199, 200, 201, 200, 200],
                        'total_counts': [1000] * 5, 'pct_counts_mt': [11, 10, 10.1, 0, 0]})
    flags = qc_flags(obs)
    assert flags.low_genes.tolist() == [True, False, False, False, False]
    assert flags.high_mt.tolist() == [True, False, True, False, False]
    assert (~flags.any(axis=1)).sum() == 3


def test_marker_scoring_separates_panels_and_flags_ties():
    genes = [g for panel in MARKERS.values() for g in panel]
    matrix = np.zeros((7, len(genes)))
    for i, panel in enumerate(MARKERS.values()):
        matrix[i, [genes.index(g) for g in panel]] = 5
    matrix[6] = matrix[0] + matrix[1]
    # Equal-size panels for T/B make the mixed cluster's scores identical.
    a = ad.AnnData(csr_matrix(matrix), obs=pd.DataFrame({'leiden': list('0123456')},
                                                       index=[f'c{i}' for i in range(7)]),
                   var=pd.DataFrame(index=genes))
    labels, _, _ = marker_annotation(a)
    assert labels.loc['0', 'marker_type'] == 'T cells'
    assert labels.loc['1', 'marker_type'] == 'B cells'
    assert labels.loc['6', 'uncertain']
    assert labels.loc['6', 'marker_type'].startswith('Uncertain')
    with pytest.raises(ValueError, match='two detected'):
        marker_annotation(a[:, :2])


def test_unknown_model_label_is_explicit():
    mapped = broad_labels(pd.Series(['Naive B cells', 'pDC', 'Novel type']))
    assert mapped.tolist() == ['B cells', 'Dendritic cells', 'Unknown']


def test_mixed_cluster_is_uncertain_despite_positive_margin():
    genes = [g for panel in MARKERS.values() for g in panel]
    matrix = np.zeros((12, len(genes)))
    # Mixed cluster has four NK-like and two T-like cells; other lineages supply
    # a background. A positive average margin must not conceal heterogeneity.
    panels = list(MARKERS.values())
    for i in range(12):
        panel = panels[2] if i < 4 else panels[0] if i < 6 else panels[(i - 6) % 6]
        matrix[i, [genes.index(g) for g in panel]] = 5
    a = ad.AnnData(csr_matrix(matrix), obs=pd.DataFrame(
        {'leiden': ['mixed'] * 6 + ['background'] * 6},
        index=[f'c{i}' for i in range(12)]), var=pd.DataFrame(index=genes))
    labels, _, _ = marker_annotation(a)
    row = labels.loc['mixed']
    assert row.marker_top_type == 'NK cells'
    assert row.marker_margin > 0.15
    assert row.marker_support == pytest.approx(4 / 6)
    assert row.uncertain
