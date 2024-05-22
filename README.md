# pgBoost.R

The provided script pgBoost.R takes as input a data set of candidate regulatory links x linking scores and generates consensus linking scores using gradient boosting (in a leave-one-chromosome-out framework).

Example command: 

Rscript pgBoostR

## Arguments

Attempt | #1 | #2 | #3 | #4 | #5 | #6 | #7 | #8 | #9 | #10 | #11
--- | --- | --- | --- |--- |--- |--- |--- |--- |--- |--- |---
Seconds | 301 | 283 | 290 | 286 | 289 | 285 | 287 | 287 | 272 | 276 | 269

- data_file (tab-separated): candidate links (rows) x features (columns)
  - Should contain the following columns (additional columns will be ignored):
    - Index columns (_at least_ 1): this column / combination of columns (e.g. "SNP", "peak", "gene") should uniquely identify each candidate link (row names will be ignored). **_Must_ match columns in training_file.**
    - Predictors (_at least_ 1): features to be included in the model. **_Must_ match those provided in predictor_file (see below).**
    - Leave-one-out column (_exactly_ 1): column used to group links for the leave-one-chromosome-out approach (see below).
 
- training_file (tab-separated): training links (rows) x training classification (column).
  - Should contain the following columns:
    - Index columns (see above): this column / combination of columns (e.g. "SNP", "peak", "gene") should uniquely identify each candidate link (row names will be ignored). **_Must_ match columns in data_file.** Links not included in data_file will be ignored.
  - "positive": a binary indicator column (0/1) indicating whether a link is a positive (1) or negative (0) training instance.
 

 
- predictor_file: text file with one predictor name per line (corresponding to the name of a column in data_file).

- drop_duplicates_file: text file specifying columns across which to drop duplicate instances during training (one column name per line). Any duplicate training instances containing the exact same value in these columns as another instance (and also identically classified as positive/negative) will be dropped.

- LOO_colname: name of the column denoting the chromosome of the focal link - e.g. "chr" (default), "CHR", "chrom", "chromosome" - used to conduct the leave-one-chromosome-out training procedure.

- outfile: name of file where pgBoost predictions will be saved as a tab-delimited file with all index columns plus a "pgBoost" column.
