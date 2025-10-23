suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(xgboost))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))

parser <- ArgumentParser()

parser$add_argument("--data_file",
    help="tsv file containing data")
parser$add_argument("--training_file",
    help="tsv file containing training data (see description for details)")
parser$add_argument("--predictor_file",
    help="file containing names of model predictors (should correspond to columns in data_file")
parser$add_argument("--drop_duplicates_file",
    help="file containing names of predictors across which duplicates should be dropped")
parser$add_argument("--distance_colname",
    default = "snp_gene_distance", 
    help="Name of column containing distance (SNP-gene, peak-gene, etc.)")
parser$add_argument("--promoter_distance_to_TSS",
    default = 1000,
    help="minimum candidate link distance to use for training")
parser$add_argument("--output_promoter_links",
    action = "store_true",
    help="supply to include promoter-gene links in output file (assigned probability 1)")
parser$add_argument("--LOO_colname", default = "chr",
    help="name of column (chr/chrom/chromosome) used to separate data for leave-one-out training procedure")
parser$add_argument("--gene_colname", default = "gene",
    help="name of column (gene/Gene/GeneSymbol/etc.) that denotes gene name")
parser$add_argument("--element_colname", default = "gene",
    help="name of column (snp/peak/ElementName/etc.) that denotes element name")
parser$add_argument("--outfile", default = "pgboost_predictions.tsv",
    help="name of output file for pgBoost predictions")
parser$add_argument("--save_model",
    action = "store_true",
    help="supply to save a .model file (in the same folder as <outfile>")
parser$add_argument("--seed", default = 511,
    help="value to pass to set.seed()")

args <- parser$parse_args()

data_file = args$data_file
training_file = args$training_file
output_promoter_links = args$output_promoter_links
predictor_file = args$predictor_file
drop_duplicates_file = args$drop_duplicates_file
distance_colname = args$distance_colname
promoter_distance_to_TSS = as.numeric(args$promoter_distance_to_TSS)
LOO_colname = args$LOO_colname
gene_colname = args$gene_colname
element_colname = args$element_colname
outfile = args$outfile
save_model = args$save_model
seed = args$seed

set.seed(seed)

# Ensure that required arguments are supplied
if (is.null(data_file)) {
    writeLines("File with candidate link data not supplied! Halting.")
    quit()
    }

if (is.null(training_file)) {
    writeLines("File with training data not supplied! Halting.")
    quit()
    }

if (is.null(predictor_file)) {
    writeLines("File with predictors not supplied! Halting.")
    quit()
    }

writeLines("\nArguments:\n")
writeLines(sprintf("data_file: %s", data_file))
writeLines(sprintf("training_file: %s", training_file))
writeLines(sprintf("predictor_file: %s", predictor_file))
writeLines(sprintf("drop_duplicates_file: %s", max(drop_duplicates_file, "None")))
writeLines(sprintf("gene_colname: %s", gene_colname))
writeLines(sprintf("element_colname: %s", element_colname))
writeLines(sprintf("distance_colname: %s", distance_colname))
writeLines(sprintf("LOO_colname: %s", LOO_colname))
writeLines(sprintf("outfile: %s\n", outfile))

# Read in candidate link data and training data
writeLines("Reading candidate link data...")
full_data = read.csv(data_file, sep = "\t", comment = "#")
writeLines("Reading training data...\n")
training_data = read.csv(training_file, sep = "\t")
index_cols = colnames(training_data)[colnames(training_data) %in% colnames(full_data)]

# Subset to links passing the supplied distance threshold
if (!is.null(promoter_distance_to_TSS)) {
    promoter_links = full_data[full_data[distance_colname] <= promoter_distance_to_TSS, index_cols]
    writeLines(sprintf("Dropping all links involving a peak < %s bp to any TSS. From %s candidate links...", promoter_distance_to_TSS, nrow(full_data)))
    promoter_peaks = promoter_links[,element_colname]
    full_data = full_data[!full_data[,element_colname] %in% promoter_peaks,]
    writeLines(sprintf("...to %s candidate links.", nrow(full_data)))
}

# Merge features and training data
writeLines("Index columns (shared between candidate link and training data):")
print(index_cols)
writeLines("Merging candidate link and training data on index columns...\n")
data = merge(full_data, training_data, by = index_cols)

# Subset to genes with a positive link in the merged set
positive_genes = unique(data[data$positive == 1,gene_colname])
data = data[data[,gene_colname] %in% positive_genes,]
# Order data
positive_set = data[data$positive == 1,]
negative_set = data[data$positive == 0,]
data = rbind(positive_set, negative_set)

# Drop duplicate instances (same values across specified columns and training data classification)
if (!is.null(drop_duplicates_file)) {
    writeLines("Dropping duplicates across specified columns!\n")
    drop_duplicates_colnames = read.table(drop_duplicates_file)$V1
    data = distinct_at(data, c(drop_duplicates_colnames, "positive"), .keep_all = TRUE)
    }

if (nrow(data) == 0) {
    writeLines("No overlap between candidate links and training data! Halting.")
    quit()
    }

# Read in names of predictors to use in model (should correspond to column names)
predictors = read.table(predictor_file)$V1
writeLines(sprintf("Using %s features for model training and prediction.\n", length(predictors)))

if (!all(predictors %in% colnames(data))) {
    writeLines("Not all specified predictors are supplied! Halting.")
    quit()
    }

# Create data matrix and training labels
annotation_feature_tabb = data[predictors]
labels = data$positive

writeLines("Data matrix head:")
head(annotation_feature_tabb)

# Generate LOCO predictions (for candidate links with each value in the specified LOO column)
if (!(LOO_colname %in% colnames(data))) {
    writeLines("LOO column name not in data! Halting.")
    quit()
    }

chromosomes = unique(full_data[,LOO_colname])
predictions = data.frame()
for (chromosome in chromosomes) {
    writeLines(sprintf("\nChromosome %s:", chromosome))
    # Subset features and training labels to off-chromosome candidate links
    training_idx = which(data[,LOO_colname] != chromosome)
    training_data = annotation_feature_tabb[training_idx, ]
    training_labels = labels[training_idx]
    
    # Train pgBoost model on off-chromosome data
    NUMITER_XGBoost=1000
    bstSparse <-  xgboost(data = as.matrix(training_data),
                          label = training_labels,
                          max_depth = c(10, 15, 25),
                          learning_rate = 0.05,
                          gamma = 10,
                          min_child_weight = c(6, 8, 10),
                          nthread = 2,
                          scale_pos_weight = 1,
                          subsample = c(0.6, 0.8, 1),
                          nrounds = NUMITER_XGBoost,
                          objective = "binary:logistic",
                          eval_metric = "auc",
                          verbose = 0
    )
    
    # Generate predictions for candidate links on the focal chromosome
    prediction_idx = which(full_data[,LOO_colname] == chromosome)
    prediction_features = full_data[prediction_idx, predictors]
    pgBoost = predict(bstSparse, as.matrix(prediction_features))
    chrom_predictions = cbind(full_data[prediction_idx, index_cols], pgBoost)
    
    # Add predictions from focal chromosome to full set of predictions
    predictions = rbind(predictions, chrom_predictions)
    writeLines(sprintf("Number of scored links: %s\n", nrow(chrom_predictions)))
    print(head(predictions))
}

# Add promoter links to output
if (output_promoter_links == TRUE) {
    promoter_links$pgBoost = 1
    writeLines(sprintf("Number of promoter links (assigned score of 1): %s\n", nrow(promoter_links)))
    head(promoter_links)
    head(predictions)
    predictions = rbind(predictions, promoter_links)
}

# Create directory for output
out_dir = dirname(outfile)
if (!dir.exists(out_dir)) {
    print(sprintf("Creating directory"))
    dir.create(out_dir, recursive = TRUE)
}

# Save trained model to outfile directory
if (save_model == TRUE) {
    model_outfile = sprintf("%s/pgboost_trained_model.model", out_dir)
    xgb.save(bstSparse, model_outfile)
    writeLines(sprintf("Saved trained model to %s.", model_outfile))
}

# Save predictions to specified outfile name
write.table(predictions, outfile, sep = "\t", row.names = F, quote = F)
writeLines(sprintf("Total number of scored links: %s", nrow(predictions)))
writeLines(sprintf("Saved predictions to %s.", outfile))
