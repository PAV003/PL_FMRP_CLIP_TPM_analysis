# ==============================================
# R Script: featureCounts Pipeline with TPM
# ==============================================

# ==== 0. Load packages ====
library(Rsubread)
library(GenomicRanges)
library(S4Vectors)
library(rtracklayer)
library(Rsamtools)
library(BiocManager)

list(bam_files)

# ==== 1. Set directories and files ====
setwd("/Users/pedro/Desktop/PL fmrp clip all")
bam_files <- list.files(pattern = "\\.bam$", full.names = TRUE)
gtf_file  <- "/Users/pedro/Desktop/NDEL pan/Mus_musculus.GRCm39.113.gtf"

# ==== 2. Import GTF and extract features ====
cat("Importing GTF...\n")
gtf <- import(gtf_file)

exons     <- gtf[gtf$type == "exon"]
fiveUTR   <- gtf[gtf$type == "five_prime_utr"]
threeUTR  <- gtf[gtf$type == "three_prime_utr"]

# ==== 3. Generate introns ====
if (file.exists("introns.rds")) {
  cat("Loading saved introns...\n")
  introns <- readRDS("introns.rds")
} else {
  cat("Generating introns per gene...\n")
  genes <- split(exons, exons$gene_id)
  introns_list <- vector("list", length(genes))
  
  for (i in seq_along(genes)) {
    exon_gr <- genes[[i]]
    gene_range <- range(exon_gr)
    exons_reduced <- reduce(exon_gr)
    intron_ranges <- setdiff(gene_range, exons_reduced)
    
    if (length(intron_ranges) > 0) {
      intron_ranges <- GRanges(intron_ranges)
      mcols(intron_ranges)$gene_id <- exon_gr$gene_id[1]
      mcols(intron_ranges)$type <- "intron"
      introns_list[[i]] <- intron_ranges
    }
    if (i %% 1000 == 0) cat("Processed gene", i, "of", length(genes), "\n")
  }
  introns <- do.call(c, introns_list)
  saveRDS(introns, file = "introns.rds")
}
cat("Introns ready:", length(introns), "\n")

# ==== 4. Generate promoters (2kb upstream of TSS) ====
if (file.exists("promoters.rds")) {
  cat("Loading saved promoters...\n")
  promoters <- readRDS("promoters.rds")
} else {
  cat("Generating promoters...\n")
  
  # Merge exons per gene and reduce
  genes <- split(exons, exons$gene_id)
  gene_ranges_list <- lapply(genes, reduce)
  gene_ranges <- do.call(c, gene_ranges_list)
  
  # Create promoters 2kb upstream
  promoters <- promoters(gene_ranges, upstream = 2000, downstream = 0)
  
  # Assign gene_id and type
  mcols(promoters)$gene_id <- rep(names(genes), times = elementNROWS(promoters))
  mcols(promoters)$type <- "promoter"
  
  # Save
  saveRDS(promoters, file = "promoters.rds")
}
cat("Promoters ready:", length(promoters), "\n")


# ==== 5. Flatten all features & ensure gene_id ====
flatten_and_fix <- function(gr) {
  if (is(gr, "GRangesList")) gr <- unlist(gr)
  if (is.null(mcols(gr)$gene_id)) mcols(gr)$gene_id <- names(gr)
  return(gr)
}

exons     <- flatten_and_fix(exons); mcols(exons)$type <- "exon"
introns   <- flatten_and_fix(introns)
fiveUTR   <- flatten_and_fix(fiveUTR); mcols(fiveUTR)$type <- "fiveUTR"
threeUTR  <- flatten_and_fix(threeUTR); mcols(threeUTR)$type <- "threeUTR"
promoters <- flatten_and_fix(promoters); mcols(promoters)$type <- "promoter"

all_features <- c(exons, introns, fiveUTR, threeUTR, promoters)
cat("Total raw features:", length(all_features), "\n")

# ==== 6. Create SAF for featureCounts ====
fc_annot <- data.frame(
  GeneID = as.character(mcols(all_features)$gene_id),
  Chr    = as.character(seqnames(all_features)),
  Start  = as.integer(start(all_features)),
  End    = as.integer(end(all_features)),
  Strand = as.character(strand(all_features)),
  Type   = as.character(mcols(all_features)$type), # keep track of feature type
  stringsAsFactors = FALSE
)

fc_annot$Strand <- ifelse(fc_annot$Strand %in% c("+","-"), fc_annot$Strand, "*")
fc_annot <- fc_annot[!duplicated(fc_annot), ]

# Match BAM chromosome names
bam_header <- scanBamHeader(bam_files[1])[[1]]$targets
if (!grepl("^chr", names(bam_header)[1])) {
  fc_annot$Chr <- gsub("^chr", "", fc_annot$Chr)
} else {
  fc_annot$Chr <- ifelse(grepl("^chr", fc_annot$Chr), fc_annot$Chr, paste0("chr", fc_annot$Chr))
}
fc_annot <- fc_annot[fc_annot$Chr %in% names(bam_header), ]

write.table(fc_annot, "all_features_clean.saf", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
cat("SAF saved as 'all_features_clean.saf'\n")

# ==== 7. Run featureCounts ====
all_counts <- featureCounts(
  files = bam_files,
  annot.ext = "all_features_clean.saf",
  isGTFAnnotationFile = FALSE,
  useMetaFeatures = TRUE,
  countMultiMappingReads = FALSE,
  strandSpecific = 1,
  nthreads = 4,
  isPairedEnd = FALSE
)

# ==== 8. Save counts ====
write.csv(all_counts$counts, "all_counts_featureCounts.csv", row.names = TRUE)
write.csv(all_counts$stat, "featureCounts_summary.csv", row.names = TRUE)
cat("Counts and summary saved.\n")

# ==== 9. Compute feature lengths per gene (robust) ====
# Keep only features with a gene_id
valid_features <- all_features[!is.na(mcols(all_features)$GeneID) & mcols(all_features)$GeneID != ""]

# Split by gene and sum widths
feature_lengths <- sapply(split(width(valid_features), mcols(valid_features)$GeneID), sum)

# Keep only genes present in counts
common_genes <- intersect(rownames(all_counts$counts), names(feature_lengths))
counts_subset <- all_counts$counts[common_genes, ]
feature_lengths_subset <- feature_lengths[common_genes]

# ==== 10. Calculate TPM ====
# 1. Ensure all_features has GeneID
mcols(all_features)$GeneID <- as.character(mcols(all_features)$gene_id)

# 2. Filter only features with valid GeneID
valid_features <- all_features[!is.na(mcols(all_features)$GeneID) & mcols(all_features)$GeneID != ""]

# 3. Split by GeneID
features_by_gene <- split(width(valid_features), mcols(valid_features)$GeneID)

# 4. Sum widths per gene to get total feature length
feature_lengths <- sapply(features_by_gene, sum)

# 5. Keep only genes present in counts
common_genes <- intersect(rownames(all_counts$counts), names(feature_lengths))
counts_subset <- all_counts$counts[common_genes, ]
feature_lengths_subset <- feature_lengths[common_genes]

# 6. Check
str(feature_lengths_subset)
head(feature_lengths_subset)

# 1. Convert lengths to kilobases
gene_kb <- feature_lengths_subset / 1000  # kb

# 2. Calculate Reads Per Kilobase (RPK)
rpk <- counts_subset / gene_kb

# 3. Compute scaling factors per sample (sum of RPKs)
scaling_factors <- colSums(rpk)

# 4. Calculate TPM
tpm <- sweep(rpk, 2, scaling_factors, FUN = "/") * 1e6

# 5. Save TPM table
write.csv(tpm, "all_counts_TPMalltrials.csv", row.names = TRUE)

# 6. Quick sanity check: total TPM per sample should be ~1e6
colSums(tpm)



# List of feature types you want exon 5utr 3utr promoter intron
feature_types <- c("exon", "intron", "fiveUTR", "threeUTR", "promoter")

# Initialize list to store TPM per type
tpm_by_type <- list()

for (ft in feature_types) {
  
  cat("Processing", ft, "...\n")
  
  # 1. Subset features of this type
  features_ft <- all_features[mcols(all_features)$type == ft]
  
  # 2. Ensure GeneID exists
  mcols(features_ft)$GeneID <- as.character(mcols(features_ft)$gene_id)
  
  # 3. Filter valid GeneID
  valid_ft <- features_ft[!is.na(mcols(features_ft)$GeneID) & mcols(features_ft)$GeneID != ""]
  
  # 4. Split widths by gene
  widths_by_gene <- split(width(valid_ft), mcols(valid_ft)$GeneID)
  
  # 5. Sum widths per gene
  feature_lengths_ft <- sapply(widths_by_gene, sum)
  
  # 6. Keep only genes present in counts
  common_genes_ft <- intersect(rownames(counts_subset), names(feature_lengths_ft))
  counts_ft <- counts_subset[common_genes_ft, ]
  lengths_ft <- feature_lengths_ft[common_genes_ft]
  
  # 7. Convert lengths to kb
  gene_kb <- lengths_ft / 1000
  
  # 8. Compute RPK
  rpk <- counts_ft / gene_kb
  
  # 9. Compute TPM
  scaling_factors <- colSums(rpk)
  tpm_ft <- sweep(rpk, 2, scaling_factors, FUN = "/") * 1e6
  
  # 10. Store in list
  tpm_by_type[[ft]] <- tpm_ft
  
  # Optional: quick sanity check
  cat("Total TPM per sample for", ft, ":", colSums(tpm_ft), "\n")
}

# Save all TPMs
for (ft in names(tpm_by_type)) {
  write.csv(tpm_by_type[[ft]], paste0("TPM_", ft, ".csv"), row.names = TRUE)
}

# put it all in one sheet List of feature types
feature_types <- c("exon", "intron", "fiveUTR", "threeUTR", "promoter")

# Initialize list to store TPM tables
tpm_list <- list()

for (ft in feature_types) {
  
  cat("Processing", ft, "...\n")
  
  # Subset features of this type
  features_ft <- all_features[mcols(all_features)$type == ft]
  mcols(features_ft)$GeneID <- as.character(mcols(features_ft)$gene_id)
  valid_ft <- features_ft[!is.na(mcols(features_ft)$GeneID) & mcols(features_ft)$GeneID != ""]
  
  # Split widths by gene
  widths_by_gene <- split(width(valid_ft), mcols(valid_ft)$GeneID)
  
  # Sum widths per gene
  feature_lengths_ft <- sapply(widths_by_gene, sum)
  
  # Keep only genes present in counts
  common_genes_ft <- intersect(rownames(counts_subset), names(feature_lengths_ft))
  counts_ft <- counts_subset[common_genes_ft, ]
  lengths_ft <- feature_lengths_ft[common_genes_ft]
  
  # Convert lengths to kb
  gene_kb <- lengths_ft / 1000
  
  # Compute RPK
  rpk <- counts_ft / gene_kb
  
  # Compute TPM
  scaling_factors <- colSums(rpk)
  tpm_ft <- sweep(rpk, 2, scaling_factors, FUN = "/") * 1e6
  
  # Rename rows to include feature type
  rownames(tpm_ft) <- paste0(rownames(tpm_ft), "_", ft)
  
  # Store
  tpm_list[[ft]] <- tpm_ft
}

# Combine all types into one big TPM matrix
all_tpm <- do.call(rbind, tpm_list)

# Save as CSV
write.csv(all_tpm, "TPM_all_feature_typesAlltrials.csv", row.names = TRUE)

cat("All TPMs saved as 'TPM_all_feature_types.csv'\n")

#TPM with gene symbol 
# Install biomaRt if needed
if (!requireNamespace("biomaRt", quietly = TRUE)) install.packages("biomaRt")
library(biomaRt)

# Extract unique Ensembl IDs from TPM rownames
ensembl_ids <- unique(sub("_(exon|intron|fiveUTR|threeUTR|promoter)$", "", rownames(all_tpm)))

# Connect to Ensembl
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Get gene symbols
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = mart
)

# Make a named vector for easy lookup
gene_symbol_lookup <- setNames(gene_map$mgi_symbol, gene_map$ensembl_gene_id)

# Add GeneSymbol column to TPM table
gene_ids <- sub("_(exon|intron|fiveUTR|threeUTR|promoter)$", "", rownames(all_tpm))
all_tpm_with_symbols <- cbind(GeneSymbol = gene_symbol_lookup[gene_ids], all_tpm)

# Save as CSV
write.csv(all_tpm_with_symbols, "TPM_all_feature_types_with_symbolsalltrials.csv", row.names = TRUE)

cat("TPMs with gene symbols saved as 'TPM_all_feature_types_with_symbols.csv'\n")

# fix missing feature type Extract gene ID and feature type from rownames
gene_ids <- sub("_(exon|intron|fiveUTR|threeUTR|promoter)$", "", rownames(all_tpm))
feature_type <- sub(".*_(exon|intron|fiveUTR|threeUTR|promoter)$", "\\1", rownames(all_tpm))

# Map gene symbols
gene_symbols <- gene_symbol_lookup[gene_ids]

# Build new table with GeneID, GeneSymbol, FeatureType, then TPMs
all_tpm_with_info <- cbind(
  GeneID = gene_ids,
  GeneSymbol = gene_symbols,
  FeatureType = feature_type,
  all_tpm
)

# Save as CSV
write.csv(all_tpm_with_info, "TPM_with_symbols_and_featureTypealltrials.csv", row.names = FALSE)

cat("TPMs with GeneID, GeneSymbol, and FeatureType saved as 'TPM_with_symbols_and_featureType.csv'\n")

# ==== Load feature-level TPM table ====
library(dplyr)
feature_tpm <- read.csv("TPM_with_symbols_and_featureTypealltrials.csv", stringsAsFactors = FALSE)

# Check columns
head(feature_tpm)
# Expect: GeneID, GeneSymbol, FeatureType, sample1, sample2, ...

# ==== Sum TPMs per gene ====
gene_tpm <- feature_tpm %>%
  group_by(GeneID, GeneSymbol) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")

# Optional: reorder columns with GeneSymbol first
gene_tpm <- gene_tpm %>%
  select(GeneSymbol, GeneID, everything())

# ==== Save gene-level TPM table ====
write.csv(gene_tpm, "TPM_gene_level_alltrials.csv", row.names = FALSE)

cat("Gene-level TPM saved as 'TPM_gene_level_alltrials.csv'\n")
