# Condition-specific analysis

This workflow can be run in the project subdirectory `conditions`.

## Input files

The following files must be present in the `source_data` directory
under the project root directory.
Note that the relative path is `../source_data`. 
Also, note that `<sample>` represents `DMSO`, `iBET151`, `Prostrat`, or `SAHA`.

```
GRCh38_pHR_<sample>_filtered_feature_bc_matrix.h5
GRCh38_pHR_<sample>_atac_fragments.tsv.gz
GRCh38_pHR_<sample>_atac_fragments.tsv.gz.tbi
genes.gtf
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
