### Combined Manhattan and Miami Plot

This example generates a combined Manhattan and Miami plot
to visualize replicated genetic associations identified
by csmGmm.

The figure consists of:

- **Top panel:** Manhattan plot showing replicated SNPs
  based on cumulative mean local FDR.
- **Bottom panel:** Miami plot comparing GWAS association
  statistics between African and European populations.
- **Legend:** Identifies the African, iCOGS, and csmGmm
  replicated SNPs.

The script requires the significant SNP results generated
by `example_replication_analysis.R`.

Run the replication analysis first:

```r
source("examples/example_replication_analysis.R")
```

Then generate the combined figure:

```r
source("examples/example_combined_manhattan_miami_plot.R")
```

The resulting figure will be saved to:

`results/manhattan_miami_combined_new.png`
