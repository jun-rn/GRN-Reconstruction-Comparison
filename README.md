# GRN Reconstruction from RNA-seq Expression Profiles

Comparison of three gene regulatory network inference methods — Relevance
Network, Mutual Information, and ARACNE — applied to human RNA-seq data.

## Background

Gene regulatory networks describe how transcription factors and regulatory
elements control gene expression. Different inference algorithms trade
sensitivity against specificity, and this project compares three of them on
the same dataset.

## Data

Expression data from GEO accession **GSE324044** — 64 human samples profiling
the role of cytotoxic CD8 T cells, neutrophils, and type 1 interferon
signalling in hyperinflammatory pathology in HIV-associated TB meningitis.

Data is not included in this repository. Download the raw TAR from GEO and
place it in the project root before running `R/01_preprocessing.R`.

## Pipeline

**01_preprocessing.R**
- Merges 64 per-sample RSEM files into a single TPM matrix
- Filters to genes with TPM > 1 in at least 20% of samples
- Maps Ensembl IDs to HGNC symbols via biomaRt
- Log2(TPM + 1) transform
- Retains the top 2,000 most variable genes

**02_grn_inference.R**
- Builds a mutual information matrix (Spearman estimator, `minet`)
- Relevance Network via CLR
- Mutual Information network via MRNET
- ARACNE with DPI tolerance eps = 0.05
- Top 5,000 edges retained per method

**03_visualization_stats.R**
- Renders each network with `igraph`
- Computes node/edge counts, density, average degree, clustering
  coefficient, and average path length

## Results

| Method | Nodes | Edges | Density | Avg Degree | Clustering | Avg Path Length |
|---|---|---|---|---|---|---|
| Relevance Network | 1560 | 5000 | 0.00411 | 6.41 | 0.5284 | 60.82 |
| Mutual Information | 1731 | 5000 | 0.00334 | 5.78 | 0.2113 | 1.32 |
| ARACNE | 1363 | 5000 | 0.00539 | 7.34 | 0.1691 | 2.33 |

### Relevance Network
![Relevance Network](results/relevance_network.jpg)

### Mutual Information
![Mutual Information](results/mutual_information.jpg)

### ARACNE
![ARACNE Network](results/aracne_network.jpg)

## Interpretation

Relevance Network showed the highest clustering coefficient (0.528) but a very
high average path length (60.82), indicating a fragmented network of isolated
clusters — consistent with correlation-based methods capturing indirect
co-expression.

Mutual Information recovered the most nodes (1,731) and a short average path
length (1.32), reflecting sensitivity to non-linear dependencies.

ARACNE produced the most stringent network: fewest nodes (1,363) but highest
density (0.00539) and average degree (7.34). Pruning indirect edges via the
Data Processing Inequality leaves the most confident regulatory relationships,
making it the best of the three for hub gene identification.

## Limitations

Edge selection used a fixed top-5,000 cutoff per method rather than a
statistically derived threshold (e.g. permutation-based null distribution).
Results are not validated against a known reference network, so the comparison
is topological rather than a benchmark of accuracy.

## Requirements

```r
install.packages(c("tidyverse", "igraph", "ggraph"))
BiocManager::install(c("biomaRt", "minet"))
```

## Running

```r
source("R/01_preprocessing.R")
source("R/02_grn_inference.R")
source("R/03_visualization_stats.R")
```

## License

MIT
