#!/usr/bin/env python
"""Download only the public Scanpy PBMC3k count matrix and CellTypist model."""
import argparse
import hashlib
import json
from pathlib import Path

import scanpy as sc
from celltypist import models


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=Path('data/raw'))
    args = parser.parse_args()
    args.data_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.datasetdir = args.data_dir
    adata = sc.datasets.pbmc3k()
    model_name = 'Immune_All_Low.pkl'
    models.download_models(model=model_name, force_update=False)
    model_path = Path(models.get_model_path(model_name))
    matrix_path = args.data_dir / 'pbmc3k_raw.h5ad'
    manifest = {
        'dataset': '10x Genomics PBMC3k, healthy donor, hg19',
        'dataset_url': 'https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k',
        'retrieval': 'scanpy.datasets.pbmc3k (public raw count H5AD mirror)',
        'cells': adata.n_obs, 'genes': adata.n_vars,
        'matrix_sha256': sha256(matrix_path),
        'model': model_name, 'model_sha256': sha256(model_path),
        'model_source': 'https://www.celltypist.org/models',
    }
    expected = {
        'matrix_sha256': '89a96f1beaa2dd83a687666d3f19a4513ac27a2a2d12581fcd77afed7ea653a1',
        'model_sha256': '290874d35dac039d4c9218c343fde4aac1077709b72a331ce7266f6828c36502',
    }
    for key, digest in expected.items():
        if manifest[key] != digest:
            raise RuntimeError(f'Upstream artifact changed: {key}; review before updating pin')
    (args.data_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
