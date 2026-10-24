"""Validate the actual committed analysis, including count conservation."""
import json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from PIL import Image

RESULTS = Path(__file__).resolve().parents[1] / 'results'


def test_persisted_analysis_integrity():
    s = json.loads((RESULTS / 'summary.json').read_text())
    a = ad.read_h5ad(RESULTS / 'pbmc3k_annotated.h5ad')
    qc = pd.read_csv(RESULTS / 'qc_cells.csv', index_col=0)
    clusters = pd.read_csv(RESULTS / 'cluster_annotations.csv')
    assert a.shape == (s['retained_cells'], s['retained_genes'])
    assert len(qc) == s['input_cells']
    assert s['input_cells'] == s['qc_removed_cells'] + s['predicted_doublets_removed'] + a.n_obs
    assert qc.retained.sum() == a.n_obs
    assert set(qc.index[qc.retained]) == set(a.obs_names)
    assert a.obs_names.is_unique and a.var_names.is_unique
    assert clusters.cells.sum() == a.n_obs and len(clusters) == s['clusters']
    assert np.isfinite(a.obsm['X_umap']).all()
    assert a.layers['counts'].data.min() >= 0
    assert np.equal(a.layers['counts'].data, np.floor(a.layers['counts'].data)).all()
    assert int(a.obs.broad_agreement.sum()) == s['broad_agreement_cells']
    assert s['broad_agreement_fraction'] == s['broad_agreement_cells'] / a.n_obs
    cross = pd.read_csv(RESULTS / 'annotation_crosstab.csv', index_col=0)
    assert cross.to_numpy().sum() == a.n_obs
    assert s['uncertain_cells'] == int(a.obs.uncertain.sum())
    assert int((a.obs.marker_top_type.astype(str) == a.obs.celltypist_broad.astype(str)).sum()) == s['candidate_agreement_cells']
    assert ((clusters.marker_support < 0.7) <= clusters.uncertain).all()


def test_pngs_are_readable():
    for name in ['umap.png', 'marker_dotplot.png', 'qc.png', 'doublet_scores.png']:
        with Image.open(RESULTS / name) as image:
            assert image.width > 300 and image.height > 200
            image.verify()


def test_artifact_hash_inventory_matches_files():
    import hashlib
    records = json.loads((RESULTS / 'artifact_manifest.json').read_text())
    assert len(records) == 17
    assert len({r['file'] for r in records}) == len(records)
    assert len({r['sha256'] for r in records}) == len(records)
    for record in records:
        path = RESULTS / record['file']
        assert path.stat().st_size == record['bytes']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record['sha256']
