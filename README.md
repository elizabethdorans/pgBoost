# pgBoost.R

The provided script pgBoost.R takes as input a data set of candidate gene regulatory links x link attributes and generates consensus linking scores using gradient boosting (in a leave-one-chromosome-out framework). This repo will be updated as pgBoost is developed into a software package with additional functionality.

<div align="center">
<img src="https://github.com/elizabethdorans/pgBoost/blob/main/figs/overview_method.jpg" width=75%>
</div>

### Arguments

Rscript pgBoost.R --data_file data.tsv --training_file training.tsv --predictor_file predictors.txt --drop_duplicates_file drop_duplicates.txt --LOO_colname "chrom" --outfile outfile.tsv

| Argument | Description |
| -------- | ----------- |
| __--data_file__ | A tab-separated data frame of candidate links (rows) x linking attributes (columns). Must contain all columns specified in __predictor_file__ _and_ the column specified in __LOO_colname__. Must also contain one or more columns which uniquely index candidate links (e.g. "SNP", "peak", "gene"). Additional columns will be ignored. |
| __--training_file__ | A tab-separated data frame of training links (rows) x training link attributes (columns). Must contain one column named "positive" which provides a binary indicator (0/1) of whether a link is a positive (1) or negative (0) training instance. Must contain columns which uniquely identify candidate links (e.g. "SNP", "peak", "gene") matching those in __data_file__. |
| __--predictor_file__ | A line-delimited text file containing names of columns to use as predictors. Must match columns in __data_file__. |
| __--drop_duplicates_file__ | (OPTIONAL) A line-delimited text file containing the names of columns used to drop duplicate instances during training. If more than one training instance contains the same values across these columns (and the same classification as positive/negative), only one instance will be retained. Must match columns in __data_file__. If not supplied, pgBoost will not check for duplicate training instances. |
| __--LOO_colname__ | Name of the column used to group links for the leave-one-chromosome-out framework, e.g. "chr" (default), "CHR", "chrom", "chromosome". |
| __--outfile__ | Name of tab-delimited file where pgBoost predictions will be saved (default "pgboost_predictions.tsv"). |

### Output

The pgBoost output file will include:
- All index columns uniquely identifying candidate links (e.g. "SNP", "peak", "gene").
- __pgBoost__: a column containing pgBoost predictions.
