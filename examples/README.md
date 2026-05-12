# Example files

This folder contains small example input and output files for pgBoost.

Run from the main `pgBoost` folder:

```bash
Rscript pgBoost.R \
  --data_file examples/input/data.tsv \
  --training_file resources/training_data.tsv.gz \
  --predictor_file examples/input/predictors.txt \
  --drop_duplicates_file examples/input/drop_duplicates.txt \
  --outfile examples/output/pgboost_predictions.tsv
```

The files in `input/` are intentionally small and are meant only to show the expected format.

| File | Description |
| ---- | ----------- |
| `input/data.tsv` | Candidate SNP-gene links and pgBoost features. |
| `input/predictors.txt` | Names of columns to use as model predictors. |
| `input/drop_duplicates.txt` | Names of columns used to drop duplicate training instances. |
| `input/signac_peak_gene_links.tsv` | Example Signac peak-gene links. |
| `input/scent_peak_gene_links.tsv` | Example SCENT peak-gene links. |
| `input/cicero_peak_gene_links.tsv` | Example Cicero peak-gene links. |
| `input/macs2_peaks.bed` | Example peak file. |
| `output/pgboost_predictions.tsv` | Example pgBoost output. |
