#!/usr/bin/env Rscript

# ------------------ Load Libraries ------------------
suppressMessages({
  library(DESeq2)
  library(optparse)
})

# ------------------ CLI Arguments ------------------
option_list <- list(
  make_option(c("-c", "--counts"), 
              type    = "character", 
              help    = "Path to Bracken counts CSV (samples as columns)"),
  make_option(c("-m", "--metadata"), 
              type    = "character", 
              help    = "Sample metadata TSV file (rows = samples)"),
  make_option(c("--condition"), 
              type    = "character", 
              help    = "Column name in metadata for condition (e.g. treatment group)"),
  make_option(c("--batch"), 
              type    = "character", 
              default = NULL, 
              help    = "Optional column name in metadata for batch effect"),
  make_option(c("--min-counts"), 
              type    = "numeric", 
              default = 10, 
              help    = "Discard taxa with total counts ≤ this threshold (default = 10)"),
  make_option(c("--sample-names"), 
              type    = "character", 
              default = "sample-id", 
              help    = "Column from metadata file with the Samples ID. The names should be the same as the ones that appears in '--counts' csv (default: SampleID)"),
  make_option(c("-o", "--output"), 
              type    = "character", 
              default = "deseq2_results", 
              help    = "Output file prefix")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Treat empty string as NULL for --batch
if (!is.null(opt$batch) && opt$batch == "") {
  opt$batch <- NULL
}

# ------------------ Load Data ------------------
counts   <- read.csv(opt$counts, row.names = 1, check.names = FALSE)

print("Head Counts: ")
print(head(counts))

print("Opt Values")
print(opt)

# Check if sample-names was provided and is valid
if (is.null(opt[["sample-names"]]) || opt[["sample-names"]] == "") {
  stop("❌ The --sample-names argument must be a valid column name.")
}

# Read the metadata file without setting row names to check the structure
metadata_preview <- read.table(opt$metadata, header = TRUE, sep = "\t", check.names = FALSE)
print("Column Names For Metadata: ")
print(colnames(metadata_preview))

if (!(opt[["sample-names"]] %in% colnames(metadata_preview))) {
  stop(paste("❌ Column", shQuote(opt[["sample-names"]]), "not found in metadata file."))
}

# Now load metadata with the specified column as row names
metadata <- read.table(opt$metadata, header = TRUE, sep = "\t",
                       row.names = opt[["sample-names"]], check.names = FALSE)

print("Head Metadata: ")
print(head(metadata))

# ------------------ Match Samples ------------------
common_samples <- intersect(colnames(counts), rownames(metadata))
if (length(common_samples) < 2) {
  stop("There must be at least two samples in common between counts and metadata.")
}

counts   <- counts[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# ------------------ Check Condition Column ------------------
if (!(opt$condition %in% colnames(metadata))) {
  stop(paste("Condition column", opt$condition, "not found in metadata."))
}
# Ensure the condition column is a factor
metadata[[opt$condition]] <- factor(metadata[[opt$condition]])
cond_table <- table(metadata[[opt$condition]])
if (length(cond_table) < 2) {
  stop("The condition column must have at least two levels.")
}
if (any(cond_table < 2)) {
  stop("Each level of the condition must have at least two samples.")
}

# ------------------ Check Batch Column (if Provided) ------------------
if (!is.null(opt$batch)) {
  if (!(opt$batch %in% colnames(metadata))) {
    stop(paste("Batch column", opt$batch, "not found in metadata."))
  }
  # Ensure batch is a factor as well
  metadata[[opt$batch]] <- factor(metadata[[opt$batch]])
}

# ------------------ Build Design Formula ------------------
if (!is.null(opt$batch)) {
  design_formula <- as.formula(paste("~", opt$batch, "+", opt$condition))
} else {
  design_formula <- as.formula(paste("~", opt$condition))
}

# ------------------ Print Sample Counts for Logging ------------------
message("Sample counts by condition:")
print(cond_table)

if (!is.null(opt$batch)) {
  message("Sample counts by batch:")
  print(table(metadata[[opt$batch]]))
}

# ------------------ DESeq2 Analysis ------------------
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData   = metadata,
  design    = design_formula
)

# Filter out taxa (rows) with very low counts across all samples
keep_rows <- rowSums(counts(dds)) > opt[["min-counts"]]
dds <- dds[keep_rows, ]

# Run DESeq
dds <- DESeq(dds)

# ------------------ Extract Results ------------------
condition_levels <- levels(metadata[[opt$condition]])
res <- results(
  dds, 
  contrast = c(opt$condition, condition_levels[2], condition_levels[1])
)

# ------------------ Construct Output File Names ------------------
prefix <- opt$output

csv_name <- paste0(
  prefix, 
  "_deseq2_", 
  opt$condition, 
  if (!is.null(opt$batch)) paste0("_batch_", opt$batch) else "", 
  ".csv"
)

ma_name <- paste0(opt$output, "_MAplot.pdf")

# ------------------ Save Output ------------------
write.csv(as.data.frame(res), file = csv_name)

pdf(ma_name)
plotMA(res, main = "DESeq2 MA-Plot", ylim = c(-4, 4))
dev.off()

# ------------------ Completion Message ------------------
message("✅ DESeq2 analysis complete. Files saved as:")
message("   • ", csv_name)
message("   • ", ma_name)
