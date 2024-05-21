# pgBoost.R

Download and run the script pgBoost.R on a user-supplied data frame of linking scores and distance-based features to train a gradient boosting model and generate predictions using a leave-one-chromosome-out framework.

## Arguments

- data_file: tab-separated file in which each row represents a candidate link. Should contain the following columns (but can include additional columns, which will be ignored):
  1. Index columns (e.g. "SNP", "peak", "gene"): columns which identify a candidate link. Column names _must_ be identical to the index columns in the training_file. There must be at least one index column, and each candidate link must be uniquely identified by the combination of index columns.
  2. Predictors: columns which contain features to be included in the model. Column names _must_ be identical to those provided in the predictor_file (see below).
  3. Leave-one-out column: _one_ column denoting the chromosome of the focal link (see below).
 
- training_file: tab-separated file in which each row represents a link in the training data. Should contain the following columns:
  1. Index columns (e.g. "SNP", "peak", "gene"): columns which identify a candidate link. Column names _must_ be identical to the index columns in the data_file. There must be at least one index column, and each candidate link must be uniquely identified by the combination of index columns.
     - The training data file can include candidate links that are not in the data file (these links will be discarded).
  2. "positive": a binary column (0/1) indicating whether a link is a positive (1) or negative (0) training instance.
 
- predictor_file: text file with one predictor name per line (corresponding to the name of a column in data_file).

- drop_duplicates_file: text file specifying columns across which to drop duplicate instances during training (one column name per line). Any duplicate training instances containing the exact same value in these columns as another instance (and also identically classified as positive/negative) will be dropped.

- LOO_colname: name of the column denoting the chromosome of the focal link - e.g. "chr" (default), "CHR", "chrom", "chromosome" - used to conduct the leave-one-chromosome-out training procedure.

- outfile: name of file where pgBoost predictions will be saved as a tab-delimited file with all index columns plus a "pgBoost" column.
