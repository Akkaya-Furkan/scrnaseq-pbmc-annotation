# Single-cell RNA-seq — PBMC cell-type annotation

**Question.** Which cell types do the clusters in a public PBMC dataset correspond to?

## Data
Public 10x Genomics PBMC data, analysed with scanpy.

## Approach
1. QC: genes per cell, mitochondrial fraction, doublets
2. Normalisation, HVGs, PCA, neighbours, Leiden, UMAP
3. Marker-based annotation compared with automated annotation (CellTypist)

## Deliverables
Annotated AnnData, UMAP and marker dot plot, write-up.

## Limitations
- Annotation depends on marker choice; without reference labels, consistency rather than accuracy is reported.

## Status
Planned — part of a 30-project computational biology and scientific software portfolio.
