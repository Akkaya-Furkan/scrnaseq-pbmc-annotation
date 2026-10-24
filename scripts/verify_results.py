#!/usr/bin/env python
"""Re-read deliverables and record a content hash inventory."""
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
results = root / 'results'
expected = [
    'pbmc3k_annotated.h5ad', 'umap.png', 'marker_dotplot.png', 'qc.png',
    'doublet_scores.png', 'summary.json', 'data_manifest.json', 'qc_cells.csv',
    'cluster_marker_scores.csv', 'cell_marker_scores.csv', 'celltypist_label_mapping.csv',
    'annotation_crosstab.csv', 'cluster_celltypist_labels.csv', 'cluster_annotations.csv',
    'cell_annotations.csv', 'ranked_markers.csv', 'runtime.txt',
]
missing = [name for name in expected if not (results / name).is_file()]
if missing:
    raise SystemExit(f'Missing outputs: {missing}')
records = []
for name in expected:
    path = results / name
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    records.append({'file': name, 'bytes': path.stat().st_size, 'sha256': digest})
assert len({r['sha256'] for r in records}) == len(expected)
(results / 'artifact_manifest.json').write_text(json.dumps(records, indent=2) + '\n')
print(f'{len(records)}/{len(expected)} expected outputs; {len(records)} unique SHA-256 hashes')
