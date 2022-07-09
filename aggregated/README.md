# Aggregated analyses

This workflow can be run in the project subdirectory `aggregated`.

## Input files

Create a subdirectory, `derived_data`.

```
mkdir derived_data
```

Copy the following files into `derived_data`, or create symbolic links.

```
hiv_predict.rds
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
snakemake --cores <number_of_cores>
```

## Output files

The workflow will generate the following files in `figures`.

```
umap.pdf
vln.pdf
umap_no_hiv.pdf
vln_no_hiv.pdf
```

Also, the following files will be generated in `derived_data`.

```
no_hiv.rds
```
