suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(xgboost))

parser <- ArgumentParser()

parser$add_argument("--data_file",
    required = TRUE,
    help = "tsv file containing candidate link data")
parser$add_argument("--training_file",
    required = TRUE,
    help = "tsv file containing training data")
parser$add_argument("--predictor_file",
    required = TRUE,
    help = "file containing names of model predictors")
parser$add_argument("--drop_duplicates_file",
    default = NULL,
    help = "file containing names of columns across which duplicates should be dropped")
parser$add_argument("--LOO_colname", default = "chr",
    help = "name of column used to separate data for leave-one-out training procedure")
parser$add_argument("--outfile", default = "pgboost_predictions.tsv",
    help = "name of output file for pgBoost predictions")
parser$add_argument("--seed", default = 511, type = "integer",
    help = "value to pass to set.seed()")
parser$add_argument("--nthread", default = 1, type = "integer",
    help = "number of threads to use for xgboost")

args <- parser$parse_args()

data_file <- args$data_file
training_file <- args$training_file
predictor_file <- args$predictor_file
drop_duplicates_file <- args$drop_duplicates_file
LOO_colname <- args$LOO_colname
outfile <- args$outfile
seed <- args$seed
nthread <- args$nthread

set.seed(seed)

writeLines("Arguments:\n")
writeLines(sprintf("data_file: %s", data_file))
writeLines(sprintf("training_file: %s", training_file))
writeLines(sprintf("predictor_file: %s", predictor_file))
writeLines(sprintf("drop_duplicates_file: %s", ifelse(is.null(drop_duplicates_file), "None", drop_duplicates_file)))
writeLines(sprintf("LOO_colname: %s", LOO_colname))
writeLines(sprintf("outfile: %s", outfile))
writeLines(sprintf("seed: %s", seed))
writeLines(sprintf("nthread: %s\n", nthread))

# Read in candidate link data and training data
writeLines("Reading candidate link data...")
full_data <- read.csv(data_file, sep = "\t", check.names = FALSE)
writeLines("Reading training data...\n")
training_data <- read.csv(training_file, sep = "\t", check.names = FALSE)

# Read in names of predictors to use in model
predictors <- read.table(predictor_file)$V1
writeLines(sprintf("Using %s features for model training and prediction.\n", length(predictors)))

if (!all(predictors %in% colnames(full_data))) {
    stop("Not all specified predictors are supplied in data_file.", call. = FALSE)
}

if (!(LOO_colname %in% colnames(full_data))) {
    stop("LOO column name not in data_file.", call. = FALSE)
}

if (!("positive" %in% colnames(training_data))) {
    stop("training_file must contain a column named positive.", call. = FALSE)
}

# Merge features and training data
index_cols <- colnames(training_data)[colnames(training_data) %in% colnames(full_data)]
writeLines("Index columns (shared between candidate link and training data):")
print(index_cols)
writeLines("Merging candidate link and training data on index columns...\n")

data <- merge(full_data, training_data, by = index_cols)

if (nrow(data) == 0) {
    stop("No overlap between candidate links and training data.", call. = FALSE)
}

# Subset to genes with a positive link in the merged set
positive_genes <- unique(data[data$positive == 1, ]$gene)
data <- data[data$gene %in% positive_genes, ]

if (nrow(data) == 0) {
    stop("No candidate genes with positive training links remain after filtering.", call. = FALSE)
}

# Order data
positive_set <- data[data$positive == 1, ]
negative_set <- data[data$positive == 0, ]
data <- rbind(positive_set, negative_set)

# Drop duplicate instances (same values across specified columns and training data classification)
if (!is.null(drop_duplicates_file)) {
    writeLines("Dropping duplicates across specified columns!\n")
    drop_duplicates_colnames <- read.table(drop_duplicates_file)$V1
    duplicate_key <- data[c(drop_duplicates_colnames, "positive")]
    data <- data[!duplicated(duplicate_key), ]
}

writeLines(sprintf("Number of training instances: %s", nrow(data)))
writeLines(sprintf("Number of positives: %s", sum(data$positive)))
writeLines(sprintf("Number of negatives: %s\n", nrow(data) - sum(data$positive)))

# Create data matrix and training labels
annotation_feature_tabb <- data[predictors]
labels <- data$positive

# Generate LOCO predictions (for candidate links with each value in the specified LOO column)
chromosomes <- unique(full_data[, LOO_colname])
predictions <- data.frame()

for (chromosome in chromosomes) {
    writeLines(sprintf("Chromosome %s:", chromosome))

    # Subset features and training labels to off-chromosome candidate links
    training_idx <- which(data[, LOO_colname] != chromosome)
    training_features <- annotation_feature_tabb[training_idx, , drop = FALSE]
    training_labels <- labels[training_idx]

    writeLines(sprintf("Number of training instances: %s", length(training_labels)))
    writeLines(sprintf("Number of positives: %s", sum(training_labels)))
    writeLines(sprintf("Number of negatives: %s", length(training_labels) - sum(training_labels)))

    # Train pgBoost model on off-chromosome data
    bstSparse <- xgboost(data = as.matrix(training_features),
                         label = training_labels,
                         max_depth = c(10, 15, 25),
                         learning_rate = 0.05,
                         gamma = 10,
                         min_child_weight = c(6, 8, 10),
                         nthread = nthread,
                         scale_pos_weight = 1,
                         subsample = c(0.6, 0.8, 1),
                         nrounds = 1000,
                         objective = "binary:logistic",
                         eval_metric = "auc",
                         verbose = 0
    )

    # Generate predictions for candidate links on the focal chromosome
    prediction_idx <- which(full_data[, LOO_colname] == chromosome)
    prediction_features <- full_data[prediction_idx, predictors, drop = FALSE]
    pgBoost_probability <- predict(bstSparse, as.matrix(prediction_features))
    chrom_predictions <- cbind(full_data[prediction_idx, index_cols, drop = FALSE], pgBoost_probability)

    # Add predictions from focal chromosome to full set of predictions
    predictions <- rbind(predictions, chrom_predictions)
    writeLines(sprintf("Number of scored links: %s\n", nrow(chrom_predictions)))
}

if (nrow(predictions) > 1) {
    predictions$pgBoost_percentile <- (rank(predictions$pgBoost_probability, ties.method = "min") - 1) /
        (nrow(predictions) - 1)
} else {
    predictions$pgBoost_percentile <- 0
}

# Save predictions to specified outfile name
write.table(predictions, outfile, sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(sprintf("Total number of scored links: %s", nrow(predictions)))
writeLines(sprintf("Saved predictions to %s.", outfile))
