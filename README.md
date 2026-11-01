
## Example: Breast Cancer Replication Analysis

This example demonstrates how to apply the `csmGmm` package
to identify replicated genetic associations across two
breast cancer GWAS datasets.

### Datasets

The example uses two publicly available GWAS summary
statistics datasets:

| Dataset | Ancestry | GWAS Catalog ID |
|---|---|---|
| African females | African | GCST90551892 |
| iCOGS | European | GCST010098 |

Both datasets use genome build 37.

Download the summary statistics from the
[GWAS Catalog](https://www.ebi.ac.uk/gwas/).

Place the downloaded files in the `data/` directory
using the following filenames:

- `GCST010098.tsv`
- `GCST90551892.tsv.gz`

### Running the analysis

Install and load the required packages:

```r
install.packages(c("dplyr", "readr"))
```

Install `csmGmm` following the installation instructions
provided in this repository.

Run the example script:

```r
source("examples/example_replication_analysis.R")
```

### Analysis settings

The example uses the following settings:

- Number of datasets: 2
- Z-score truncation: 8.1
- Directional constraint: `sameDirAlt = TRUE`
- EM convergence tolerance: `1e-5`
- Local FDR threshold: 0.10

### Output

The analysis generates:

- `results/csmgmm_rep_results.csv`: SNPs identified as
  replicated by csmGmm.
- Manhattan plot: Visualization of the genomic
  distribution of replicated SNPs.
- Z-score plot: Visualization of the association
  statistics across datasets.

For additional details about the methodology, see the
csmGmm package documentation.
```
