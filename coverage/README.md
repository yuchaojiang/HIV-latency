# Generating coverage plots

This workflow can be run in the project subdirectory `coverage`.

## Input files

The following files must be present in the `source_data` directory
under the project root directory.
Note that the relative path is `../source_data`. 
Also, note that `<sample>` represents `DMSO`, `iBET151`, `Prostrat`, or `SAHA`.

```
GRCh38_pHR_<sample>_gex_possorted_bam.bam
GRCh38_pHR_<sample>_atac_fragments.tsv.gz
GRCh38_pHR_<sample>_atac_fragments.tsv.gz.tbi
genes.gtf
chrNameLength.txt
```

Create a subdirectory, `derived_data`.

```
mkdir derived_data
```

Download the following input files into `derived_data`, or create
symbolic linkes.

```
hiv_predict.rds
transcripts.rds
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

