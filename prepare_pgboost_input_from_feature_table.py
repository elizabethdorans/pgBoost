import pandas as pd
import numpy as np
import pybedtools
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
from glob import glob
import argparse

### Arguments

parser = argparse.ArgumentParser()

parser.add_argument("--feature_table",
                    help = "path to feature table file")
parser.add_argument("--training_file",
                    help = "File with training data")
parser.add_argument("--training_output_file",
                    help = "path to file to output training data")

args = parser.parse_args()

feature_table = args.feature_table
training_file = args.training_file
training_output_file = args.training_output_file

print(args)
print()

### PART 1: PREPARE FEATURES

# Read in constituent scores
element_gene_links = pd.read_csv(feature_table, sep = '\t', comment = '#')
    
### PART 2: PREPARE TRAINING

# Read in training data
print("Reading in training data, subsetting to 1kG SNPs and specified gene universe...")
training_data = pd.read_csv(training_file, sep = "\t")
training_data = training_data[training_data['gene'].isin(element_gene_links['GeneSymbol'].values)]
training_data[['chr', 'start', 'end']] = training_data['snp'].str.split("-", expand = True)
training_data[["start", "end"]] = training_data[["start", "end"]].astype(int)

# Merge training data and features on peak-gene pairs
print("Merging features and training data...")
element_gene_links["chr_gene"] = element_gene_links["ElementChr"] + "-" + element_gene_links["GeneSymbol"]
training_data["chr_gene"] = training_data["chr"] + "-" + training_data["gene"]

training_X_features = pybedtools.BedTool.from_dataframe(element_gene_links[['chr_gene', 'ElementStart', 'ElementEnd', 'ElementChr', 'ElementName', 'GeneSymbol']]
                                 ).intersect(
    pybedtools.BedTool.from_dataframe(training_data[["chr_gene", "start", "end", "positive"]]), wa = True, wb = True
).to_dataframe()

training_X_features = training_X_features.rename(
    columns = {"score": "ElementName", "strand": "GeneSymbol", "blockCount": "positive"}
)[["ElementName", "GeneSymbol", "positive"]]

training_X_features = training_X_features.groupby(["ElementName", "GeneSymbol"]).agg({"positive": max}).reset_index()

### SAVING OUTPUT

print("Outputting training data...")
training_X_features.to_csv(training_output_file, sep = "\t", index = False)
print("Training data written to %s!" % training_output_file)
print()