# Combined analysis

This workflow can be run in the project subdirectory `combined`.
In the workflow, results from the aggregated and condition-specific
analyses are combined and visualized.

## Input files

Create a subdirectory, `derived_data`.

```
mkdir derived_data
```

Copy the following files into `derived_data`, or create symbolic links.
Note that `<sample>` represents `DMSO`, `iBET151`, `Prostrat`, or `SAHA`.

```
hiv_predict.rds
hiv_chromvar.rds
linked.gene.csv
linked.TF.csv
chromvar.<sample>.rds
hiv.rna.<sample>.rds
motif.obj.<sample>.rds
seurat.<sample>.rds
```

## Running the workflow

To perform a dry-run, type:

```
snakemake -np
```

To visualize the workflow, type:

```
snakemake --forceall --dag | dot -Tpdf > dag.pdf
```

To run the workflow, type:

```
snakemake -cores <number_of_cores>
```


