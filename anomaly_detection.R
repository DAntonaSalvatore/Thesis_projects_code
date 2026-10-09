library(data.table)
library(isotree)
library(ggplot2)
library(dplyr)
library(ggvenn)
library(UpSetR)
library(pheatmap)
library(grid)
library(biomaRt)
library(enrichR)



setwd("C:\\Users\\Salvatore\\Desktop\\Dottorato\\studio_DL\\Anomaly_Detection")

# ---------------------------------------------------------------------
# 1. Read and merge the five datasets
# ---------------------------------------------------------------------

my_files <- list.files()
my_files <- my_files[grep("_both_sexes.tsv.bgz$", my_files)]
my_files <- my_files[-4]
merged_dataset_cleaned <- NULL
vocabulary <- read.delim("variants_col_selected.txt")

for (file in my_files) {
  
  dataset <- fread(file)
  
  if (is.null(dataset) || nrow(dataset) == 0) {
    cat("Skipping empty file:", file, "\n")
    next
  }
  
  dataset_cleaned <- dataset[, .(variant, minor_AF, low_confidence_variant, beta, se, pval)]
  dataset_cleaned$pval <- -log10(dataset_cleaned$pval)
  setDT(dataset_cleaned)
  
  # Z-score
  dataset_cleaned[, z_score := beta / se]

  # Rename columns with disease suffix
  clean_name <- sub("_both_sexes.tsv.bgz$", "", file)
  setnames(dataset_cleaned,
           old = names(dataset_cleaned)[-1],
           new = paste0(names(dataset_cleaned)[-1], "_", clean_name))
  
  # Cumulative merge
  if (is.null(merged_dataset_cleaned)) {
    merged_dataset_cleaned <- dataset_cleaned
  } else {
    merged_dataset_cleaned <- merge(merged_dataset_cleaned, dataset_cleaned,
                                    by = "variant", all = TRUE)
  }
  cat("Merged:", file, "\n")
}

# ---------------------------------------------------------------------
# 2. Compute max -log10(p) across diseases (for reference only)
# ---------------------------------------------------------------------

p_cols <- c("pval_ALS", "pval_Alzheimer", "pval_DEMENTIA", "pval_Parkinson")
merged_dataset_cleaned[, max_logp := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = p_cols]

# ---------------------------------------------------------------------
# 3. Keep only necessary columns and remove rows with NA/Inf
# ---------------------------------------------------------------------

dataset_cleaned <- merged_dataset_cleaned
keep_cols <- grep("^(z_score_|max_logp|variant)", names(dataset_cleaned), value = TRUE)
dataset_cleaned <- dataset_cleaned[, ..keep_cols]
dataset_cleaned <- dataset_cleaned[apply(is.finite(as.matrix(dataset_cleaned[, -1, with = FALSE])), 1, all)]

# ---------------------------------------------------------------------
# 4. Build Z-score matrix
# ---------------------------------------------------------------------

z_cols   <- grep("^z_score_", names(dataset_cleaned), value = TRUE)
z_mat    <- as.matrix(dataset_cleaned[, ..z_cols])

# ---------------------------------------------------------------------
# 5. Isolation Forest on Z-score matrix (discovery)
# ---------------------------------------------------------------------

set.seed(42)
iso_z <- isolation.forest(z_mat, ntrees = 100, sample_size = 256, nthreads = 4)
scores_z <- predict(iso_z, z_mat)
dataset_cleaned$anomaly_z_score <- scores_z



# ---------------------------------------------------------------------
# 7. Compute classical outlier scores on the FULL pre-clumping data
# ---------------------------------------------------------------------

# Mahalanobis
maha_dist <- mahalanobis(
  z_mat,
  colMeans(z_mat),
  cov(z_mat) + diag(1e-6, ncol(z_mat))
)
dataset_cleaned$anomaly_maha <- maha_dist

# ---------------------------------------------------------------------
# PCA reconstruction error with anomaly‑aware k selection
# ---------------------------------------------------------------------

pca <- prcomp(z_mat, center = TRUE, scale. = TRUE)
var_explained <- pca$sdev^2 / sum(pca$sdev^2)
cum_var <- cumsum(var_explained)

# Option A: Lower threshold (90%) + cap at ncol-1
target_variance <- 0.90
optimal_k <- which(cum_var >= target_variance)[1]
max_k <- ncol(z_mat) - 1
if (is.na(optimal_k) || optimal_k > max_k) optimal_k <- max_k
if (optimal_k < 1) optimal_k <- 1

cat("PCA anomaly detection: using k =", optimal_k, 
    "(explained variance =", round(cum_var[optimal_k]*100, 1), "%)\n")

# Alternative: if you prefer Option 2 (hard cap 2), uncomment:
# optimal_k <- min(2, ncol(z_mat) - 1)

# Compute reconstruction error with chosen k
z_scaled <- scale(z_mat, center = pca$center, scale = pca$scale)
reconstructed_scaled <- z_scaled %*% pca$rotation[, 1:optimal_k] %*% t(pca$rotation[, 1:optimal_k])
reconstruction_error <- rowSums((z_scaled - reconstructed_scaled)^2)
dataset_cleaned$pca <- reconstruction_error

# Univariate: max absolute Z-score across diseases
dataset_cleaned$anomaly_univ <- apply(z_mat, 1, function(x) max(abs(x)))


# ---------------------------------------------------------------------
# Store all method scores together
# ---------------------------------------------------------------------

method_scores <- list(
  IF_Z = dataset_cleaned$anomaly_z_score,
  Mahalanobis = dataset_cleaned$anomaly_maha,
  PCA = dataset_cleaned$pca,
  Univariate = dataset_cleaned$anomaly_univ
)


# ---------------------------------------------------------------------
# 8. Define top 1% for each method (on FULL pre-clumping data)
# ---------------------------------------------------------------------

top1pct <- function(scores) {
  N <- length(scores)
  rank(-scores, ties.method = "first") <= ceiling(0.01 * N)
}

dataset_cleaned$top_z       <- top1pct(dataset_cleaned$anomaly_z_score)
dataset_cleaned$top_maha    <- top1pct(dataset_cleaned$anomaly_maha)
dataset_cleaned$top_pca     <- top1pct(dataset_cleaned$pca)
dataset_cleaned$top_univ    <- top1pct(dataset_cleaned$anomaly_univ)

# ---------------------------------------------------------------------
# 9. Function to prune a set of SNPs using LD pruning
# ---------------------------------------------------------------------

prune_snps <- function(variant_ids, vocabulary, method_name) {
  
  # Prepare SNP list
  dt <- data.table(variant = variant_ids)
  dt <- merge(vocabulary, dt, by = "variant")
  
  prune_input <- data.table(SNP = unique(dt$rsid))
  
  setwd("C:\\Users\\Salvatore\\Desktop\\Plink 1.9")
  
  write.table(
    prune_input,
    file = paste0("prune_", method_name, ".txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  # PLINK pruning
  plink_cmd <- sprintf(
    paste(
      "plink",
      "-memory 30000",
      "-bfile g1000p3_EUR",
      "--extract prune_%s.txt",
      "--indep-pairwise 50 5 0.2",
      "--out prune_%s"
    ),
    method_name,
    method_name
  )
  
  system(plink_cmd)
  
  # Read pruned SNPs
  pruned_file <- paste0("prune_", method_name, ".prune.in")
  
  if (file.exists(pruned_file)) {
    
    pruned_rsids <- fread(pruned_file, header = FALSE)$V1
    
    pruned_variants <- vocabulary[
      rsid %in% pruned_rsids,
      variant
    ]
    
    return(pruned_variants)
    
  } else {
    
    warning(paste("No pruned SNPs for", method_name))
    
    return(character(0))
  }
}

# ---------------------------------------------------------------------
# 10. Clump each method's top 1% separately (using the function above)
# ---------------------------------------------------------------------

# Prepare vocabulary as data.table for merging
setDT(vocabulary)

# For IF Z (discovery)
top_z_variants <- dataset_cleaned$variant[dataset_cleaned$top_z]
top_z_scores   <- dataset_cleaned$anomaly_z_score[dataset_cleaned$top_z]
clumped_z <- prune_snps(top_z_variants, vocabulary, "IF_Z")

# For Mahalanobis
top_maha_variants <- dataset_cleaned$variant[dataset_cleaned$top_maha]
top_maha_scores   <- dataset_cleaned$anomaly_maha[dataset_cleaned$top_maha]
clumped_maha <- prune_snps(top_maha_variants, vocabulary, "Mahalanobis")

# For PCA
top_pca_variants <- dataset_cleaned$variant[dataset_cleaned$top_pca]
top_pca_scores   <- dataset_cleaned$pca[dataset_cleaned$top_pca]
clumped_pca <- prune_snps(top_pca_variants, vocabulary, "PCA")

# For Univariate
top_univ_variants <- dataset_cleaned$variant[dataset_cleaned$top_univ]
top_univ_scores   <- dataset_cleaned$anomaly_univ[dataset_cleaned$top_univ]
clumped_univ <- prune_snps(top_univ_variants, vocabulary, "Univariate")

# ---------------------------------------------------------------------
# 11. Create final datasets (only clumped SNPs that were originally top 1%)
# ---------------------------------------------------------------------

# Discovery set (IF Z)
discovery_final <- dataset_cleaned[variant %in% clumped_z & top_z == TRUE, ]

# For each classical method, create a separate final set (optional)
maha_final <- dataset_cleaned[variant %in% clumped_maha & top_maha == TRUE, ]
pca_final  <- dataset_cleaned[variant %in% clumped_pca & top_pca == TRUE, ]
univ_final <- dataset_cleaned[variant %in% clumped_univ & top_univ == TRUE, ]



# ---------------------------------------------------------------------
# 12. Compare overlaps AFTER clumping
# ---------------------------------------------------------------------

jaccard <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(0)
  intersect_len <- length(intersect(a, b))
  union_len <- length(union(a, b))
  if (union_len == 0) return(0)
  intersect_len / union_len
}

# Get variant lists from each clumped set
z_vars    <- discovery_final$variant
maha_vars <- maha_final$variant
pca_vars  <- pca_final$variant
univ_vars <- univ_final$variant

# ---------------------------------------------------------------------
# Pairwise overlaps across all methods
# ---------------------------------------------------------------------

method_sets <- list(
  IF_Z = z_vars,
  Mahalanobis = maha_vars,
  PCA = pca_vars,
  Univariate = univ_vars
)

method_names <- names(method_sets)

jaccard_matrix <- matrix(
  NA,
  nrow = length(method_sets),
  ncol = length(method_sets),
  dimnames = list(method_names, method_names)
)

for (i in seq_along(method_sets)) {
  
  for (j in seq_along(method_sets)) {
    
    jaccard_matrix[i, j] <- jaccard(
      method_sets[[i]],
      method_sets[[j]]
    )
  }
}

print(round(jaccard_matrix, 3))

write.table(jaccard_matrix, file = "Overlapping_clumped_SNPs_for_each_method.txt", quote = F, row.names = F, sep = "\t")


# Rename methods for nicer labels
rownames(jaccard_matrix) <- c(
  "Isolation Forest",
  "Mahalanobis",
  "PCA",
  "Univariate"
)

colnames(jaccard_matrix) <- c(
  "Isolation Forest",
  "Mahalanobis",
  "PCA",
  "Univariate"
)

pheatmap(
  jaccard_matrix,
  display_numbers = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_row = 0,
  treeheight_col = 0,
  border_color = NA,
  angle_col = 45,
  filename = "jaccard_heatmap.png"   # saves directly
)

#grid.newpage()
#grid.draw(p$gtable)


# ---------------------------------------------------------------------
# 13. Common SNPs across all four methods (clumped sets)
# ---------------------------------------------------------------------

common_all <- Reduce(intersect, list(z_vars, maha_vars, pca_vars, univ_vars))
cat("\nCommon SNPs across all four methods (post-clumping):\n")
print(common_all)
cat("Number of common SNPs:", length(common_all), "\n")

# ---------------------------------------------------------------------
# Soft consensus across methods (3 or more models agree)
# ---------------------------------------------------------------------

# Rebuild using clumped variant lists
all_clumped_variants <- unique(c(z_vars, maha_vars, pca_vars, univ_vars))

consensus_clumped <- data.table(variant = all_clumped_variants)
consensus_clumped[, vote_z    := variant %in% z_vars]
consensus_clumped[, vote_maha := variant %in% maha_vars]
consensus_clumped[, vote_pca  := variant %in% pca_vars]
consensus_clumped[, vote_univ := variant %in% univ_vars]

consensus_clumped[, consensus_count := vote_z + vote_maha + vote_pca + vote_univ]

# Soft consensus on clumped SNPs
consensus_clumped_3 <- consensus_clumped[consensus_count >= 3]
cat("Clumped SNPs supported by >=3 methods:", nrow(consensus_clumped_3), "\n")





# ---------------------------------------------------------------------
# 14. INTERNAL VALIDATION: Leave-One-Disease-Out (LODO)
# ---------------------------------------------------------------------

diseases <- c("ALS", "Alzheimer", "DEMENTIA", "Parkinson")

# --------------------------------------------------
# Method-specific scoring
# --------------------------------------------------

compute_method_scores <- function(method_name, mat) {
  
  if (method_name == "IF_Z") {
    
    model <- isolation.forest(
      mat,
      ntrees = 100,
      sample_size = 256,
      nthreads = 4
    )
    
    return(predict(model, mat))
  }
  
  if (method_name == "Mahalanobis") {
    
    return(
      mahalanobis(
        mat,
        colMeans(mat),
        cov(mat) + diag(1e-6, ncol(mat))
      )
    )
  }
  
  if (method_name == "PCA") {
    
    pca_tmp <- prcomp(mat, center = TRUE, scale. = TRUE)
    
    var_exp <- pca_tmp$sdev^2 / sum(pca_tmp$sdev^2)
    cum_var <- cumsum(var_exp)
    
    k <- which(cum_var >= 0.90)[1]
    
    max_k <- ncol(mat) - 1
    
    if (is.na(k) || k > max_k) k <- max_k
    if (k < 1) k <- 1
    
    scaled <- scale(
      mat,
      center = pca_tmp$center,
      scale = pca_tmp$scale
    )
    
    recon <- scaled %*%
      pca_tmp$rotation[, 1:k] %*%
      t(pca_tmp$rotation[, 1:k])
    
    return(rowSums((scaled - recon)^2))
  }
  
  if (method_name == "Univariate") {
    
    return(apply(mat, 1, function(x) max(abs(x))))
  }
}

# --------------------------------------------------
# Generalized LODO
# --------------------------------------------------

get_lodo_anomalies <- function(method_name,
                               holdout_disease,
                               top_percent = 0.01) {
  
  train_cols <- grep(
    paste0(
      "z_score_(",
      paste(setdiff(diseases, holdout_disease),
            collapse = "|"),
      ")"
    ),
    names(dataset_cleaned),
    value = TRUE
  )
  
  train_mat <- as.matrix(dataset_cleaned[, ..train_cols])
  
  scores <- compute_method_scores(
    method_name,
    train_mat
  )
  
  top_flag <- rank(
    -scores,
    ties.method = "first"
  ) <= ceiling(top_percent * nrow(dataset_cleaned))
  
  top_variants <- dataset_cleaned$variant[top_flag]
  top_scores <- scores[top_flag]
  
  clumped <- prune_snps(
    top_variants,
    vocabulary,
    paste0(method_name, "_LODO_", holdout_disease)
  )
  
  return(clumped)
}

# --------------------------------------------------
# Final SNP sets for each method
# --------------------------------------------------

final_sets <- list(
  IF_Z = discovery_final$variant,
  Mahalanobis = maha_final$variant,
  PCA = pca_final$variant,
  Univariate = univ_final$variant
)

# --------------------------------------------------
# Run LODO for ALL methods
# --------------------------------------------------

lodo_results <- list()

for (method in names(final_sets)) {
  
  cat("\nRunning LODO for:", method, "\n")
  
  lodo_sets <- lapply(
    diseases,
    function(d)
      get_lodo_anomalies(method, d)
  )
  
  names(lodo_sets) <- diseases
  
  overlaps <- sapply(
    lodo_sets,
    function(x)
      jaccard(final_sets[[method]], x)
  )
  
  lodo_results[[method]] <- data.frame(
    Method = method,
    Disease_Left_Out = diseases,
    Jaccard = overlaps
  )
}

lodo_table <- do.call(rbind, lodo_results)

print(lodo_table)

write.table(lodo_table, file = "lodo_table_stability_among_methods.txt", quote = F, row.names = F, sep = "\t")
# --------------------------------------------------
# Mean LODO stability
# --------------------------------------------------

lodo_summary <- lodo_table %>%
  group_by(Method) %>%
  summarise(
    Mean_Jaccard = round(mean(Jaccard), 3),
    SD_Jaccard = round(sd(Jaccard), 3),
    Min_Jaccard = round(min(Jaccard), 3),
    Max_Jaccard = round(max(Jaccard), 3)
  )

print(lodo_summary)

# ---------------------------------------------------------------------
# 15. RANDOM BASELINE FOR EACH METHOD
# ---------------------------------------------------------------------

set.seed(42)

random_baseline <- lapply(names(final_sets), function(method) {
  
  observed <- final_sets[[method]]
  
  rand_overlap <- replicate(100, {
    
    rand <- sample(
      dataset_cleaned$variant,
      length(observed)
    )
    
    jaccard(observed, rand)
  })
  
  data.frame(
    Method = method,
    Random_Jaccard = mean(rand_overlap)
  )
})

random_baseline <- do.call(rbind, random_baseline)

print(random_baseline)

write.table(lodo_summary, file = "LODO_summary.txt", quote = F, row.names = F, sep = "\t")
write.table(random_baseline, file = "random_baseline.txt", quote = F, row.names = F, sep = "\t")

# ---------------------------------------------------------------------
# 16. PERMUTATION TESTS FOR ALL METHODS
# ---------------------------------------------------------------------

#run_permutation_test <- function(method_name,
#                             observed_scores,
#                                 z_mat,
#                                 top_percent = 0.01,
#                                 n_perm = 200) {
#  
#  top_n <- ceiling(top_percent * nrow(z_mat))
#  
# top_idx <- order(
#    observed_scores,
#    decreasing = TRUE
#  )[1:top_n]
#  
#  z_top <- z_mat[top_idx, ]
#  
#  obs_top <- observed_scores[top_idx]
#  
#  perm_scores <- matrix(
#    NA,
#    nrow = top_n,
#    ncol = n_perm
#  )
#  
#  for (i in 1:n_perm) {
#    
#    cat(method_name, "- permutation:", i, "\n")
#    
#    perm_mat <- t(apply(z_top, 1, sample))
#    
#    perm_scores[, i] <- compute_method_scores(
#      method_name,
#      perm_mat
#    )[1:top_n]
#  }
#  
#  p_vals <- rowSums(
#    perm_scores >= obs_top
#  ) / n_perm
#  
#  p_adj <- p.adjust(
#    p_vals,
#    method = "fdr"
#  )
#  
#  validated_idx <- top_idx[p_adj < 0.05]
#  
#  return(list(
#    summary = data.frame(
#      Method = method_name,
#      Significant_SNPs = length(validated_idx),
#      Mean_Adjusted_P = mean(p_adj)
#    ),
#    
#    validated_variants = dataset_cleaned$variant[validated_idx]
#  ))
#}

# --------------------------------------------------
# Run permutation tests
# --------------------------------------------------

#perm_results <- list()

#for (method in names(method_scores)) {
#  
#  cat("\nRunning permutation test for:", method, "\n")
#  
#  perm_results[[method]] <- run_permutation_test(
#    method,
#    method_scores[[method]],
#    z_mat
#  )
#}

# --------------------------------------------------
# Summary table
# --------------------------------------------------

#perm_summary <- do.call(
#  rbind,
#  lapply(perm_results, function(x) x$summary)
#)

#print(perm_summary)


# --------------------------------------------------
# IF validated SNPs only (for downstream compatibility)
# --------------------------------------------------

#validated_if <- perm_results$IF_Z$validated_variants

#validated_if_scores <- dataset_cleaned[
#variant %in% validated_if
#]$anomaly_z_score

#validated_if_clumped <- clump_snps(
#validated_if,
#validated_if_scores,
#vocabulary,
#"IF_validated"
#)

#validated_snps <- dataset_cleaned[
#  variant %in% validated_if_clumped
#]

#dataset_cleaned <- merge(
#  dataset_cleaned,
#  vocabulary,
#  by = "variant"
#)

#validated_snps2.0 <- merge(
#  validated_snps,
#  vocabulary,
#  by = "variant"
#)

#global_medians <- apply(z_mat, 2, median)
#global_mads <- apply(z_mat, 2, mad)

#deviation_matrix <- sweep(z_mat, 2, global_medians, "-")
#deviation_matrix <- sweep(deviation_matrix, 2, global_mads, "/")

#top50_idx <- order(dataset_cleaned$anomaly_z_score,
#                   decreasing = TRUE)[1:50]

#heatmap_data <- deviation_matrix[top50_idx, ]
#rownames(heatmap_data) <- dataset_cleaned$rsid[top50_idx]

# Remove "z score " from labels
#clean_labels <- gsub("_", " ", gsub("DEMENTIA", "Dementia", colnames(heatmap_data)))
#clean_labels <- gsub("z score ", "", clean_labels)

#p <- pheatmap(
#  heatmap_data,
#  
#  cluster_rows = TRUE,
#  cluster_cols = TRUE,
#  
#  fontsize_col = 12,
# 
#  treeheight_row = 0,
#  treeheight_col = 0,
#  
#  scale = "row",
#  
#  labels_col = clean_labels,
#  angle_col = 45,
#  
#  legend = TRUE,
#  silent = TRUE  # prevents immediate plotting
#)

# Draw the plot
#grid.newpage()
#grid.draw(p$gtable)

consensus_clumped_4 <- merge(consensus_clumped_3, vocabulary, by = "variant")

# BIOMART
snp_mart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")
common_all2.0 <- consensus_clumped_4



rsids <- unique(common_all2.0$rsid)

annot <- getBM(
  attributes = c(
    "refsnp_id",
    "chr_name",
    "chrom_start",
    "allele",
    "associated_gene"
  ),
  
  filters = "snp_filter",
  
  values = rsids,
  
  mart = snp_mart
)

gene_mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)



gene_info <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "hgnc_symbol",
    "gene_biotype"
  ),
  
  filters = "ensembl_gene_id",
  
  values = unique(annot$associated_gene),
  
  mart = gene_mart
)

final_annot <- merge(
  annot,
  gene_info,
  by.x = "associated_gene",
  by.y = "ensembl_gene_id",
  all.x = TRUE
)
# Motifbreak R
snps <- common_all2.0$rsid



print(lodo_summary)
print(random_baseline)
#print(perm_summary)



# ---------------------------------------------------------------------
# 17. Write final outputs
# ---------------------------------------------------------------------

setwd("C:\\Users\\Salvatore\\Desktop\\Dottorato\\studio_DL\\Anomaly_Detection")

#write.table(validated_snps2.0,
#            "Validated_Anomalous_SNPs.txt",
#            row.names = FALSE,
#            sep = "\t",
#            quote = FALSE)

write.table(discovery_final, file = "Zscore_anomalies_discovery_clumped.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(common_all, file = "common_high_confidence_SNPs.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(consensus_clumped_3, 
            file = "SNPs_3plus_models_consensus_with_counts.txt", 
            row.names = FALSE, 
            quote = FALSE, 
            sep = "\t")

write.table(
  retro_confronto_geni_split3.0,
  "validated_snps_gene_annotation.txt",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  lodo_table,
  "LODO_results_all_methods.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  lodo_summary,
  "LODO_summary_all_methods.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  random_baseline,
  "Random_baseline_all_methods.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

#write.table(
#  perm_summary,
#  "Permutation_summary_all_methods.txt",
#  sep = "\t",
#  row.names = FALSE,
#  quote = FALSE
#)


Genes <- unique(final_annot$associated_gene)
Genes_split <- unlist(strsplit(Genes, ",\\s*"))
Genes_split <- unique(Genes_split)



dbs <- listEnrichrDbs()
# View(dbs)
selected_dbs <- c("DisGeNET", "Jensen_DISEASES_Curated_2025", "GWAS_Catalog_2023", "ClinVar_2025", "OMIM_Disease", "Human_Phenotype_Ontology")

enrich_results <- enrichr(
  Genes_split,
  selected_dbs
)

sig_results <- lapply(names(enrich_results), function(db) {
  
  res <- enrich_results[[db]]
  
  # Keep only significant rows
  res_sig <- res[res$Adjusted.P.value < 0.05, , drop = FALSE]
  
  # Skip empty results
  if (nrow(res_sig) == 0) {
    return(NULL)
  }
  
  # Add database name
  res_sig$Database <- db
  
  return(res_sig)
})
sig_table <- do.call(rbind, sig_results)
sig_table <- sig_table[
  order(sig_table$Adjusted.P.value),
]


sig_table_clean <- sig_table[, c(
  "Database",
  "Term",
  "Adjusted.P.value",
  "P.value",
  "Overlap",
  "Odds.Ratio",
  "Combined.Score",
  "Genes"
)]

diseases <- c("ALS", "Alzheimer", "DEMENTIA", "Parkinson", "amyotrophic lateral sclerosis", "dementia")


disease_related_genes <- sig_table_clean[
  sapply(sig_table_clean$Term, function(x)
    any(grepl(paste(diseases, collapse = "|"), x, ignore.case = TRUE))
  ) &
    sig_table_clean$Adjusted.P.value < 0.05,
]




write.table(
  sig_table_clean,
  "Significant_Disease_Enrichment_all_4_methods.txt",
  row.names = FALSE,
  quote = F,
  sep = "\t"
)


write.table(
  disease_related_genes,
  "Significant_Disease_Enrichment_at_least_3_methods.txt",
  row.names = FALSE,
  quote = F,
  sep = "\t"
)

retro_confronto_geni <- disease_related_genes$Genes
retro_confronto_geni_split <- unique(unlist(strsplit(retro_confronto_geni, ";")))
retro_confronto_geni_df <- data.frame(associated_gene = retro_confronto_geni_split)

retro_confronto_geni_split2.0 <- merge(retro_confronto_geni_df, final_annot, by = "associated_gene")
names(retro_confronto_geni_split2.0)[names(retro_confronto_geni_split2.0) == "refsnp_id"] <- "rsid"
retro_confronto_geni_split3.0 <- merge(retro_confronto_geni_split2.0, vocabulary, by ="rsid")
dataset_cleaned2.0 <- dataset_cleaned[,c(1:5)]

pre_heatmap <- merge(dataset_cleaned2.0, retro_confronto_geni_split3.0, by = "variant")
pre_heatmap2.0 <- pre_heatmap[,c(1:5)]


# Build matrix with variants as row names
heatmap_mat <- as.matrix(pre_heatmap2.0[, .(z_score_ALS, z_score_Alzheimer, 
                                            z_score_DEMENTIA, z_score_Parkinson)])

pre_heatmap3.0 <- merge(pre_heatmap2.0, vocabulary, by = "variant")

rownames(heatmap_mat) <- pre_heatmap3.0$rsid
heatmap_mat <- unique(heatmap_mat)

# Cleaner column labels
colnames(heatmap_mat) <- c("ALS", "Alzheimer", "Dementia", "Parkinson")

# Cluster all rows
row_dist <- dist(scale(heatmap_mat))   # Euclidean after scaling
row_hc <- hclust(row_dist, method = "ward.D2")

# Cut the dendrogram into 3 groups
row_groups <- cutree(row_hc, k = 3)



# Single heatmap for all SNPs (no splitting)
pheatmap(heatmap_mat,
         scale = "row",
         cluster_rows = TRUE,      # still cluster rows for meaningful order
         cluster_cols = TRUE,      # still cluster columns (diseases)
         treeheight_row = 0,       # removes row dendrogram
         treeheight_col = 0,       # removes column dendrogram
         angle_col = 45,           # tilts column names 45 degrees
         main = "All top anomalous",
         fontsize_row = 8,
         filename = "heatmap_all_clusters.png",
         width = 10,
         height = 12
)



all_snps <- unique(c(z_vars, maha_vars, pca_vars, univ_vars))

# Create a data frame with binary columns
upset_data <- data.frame(
  SNP = all_snps,
  Isolation_Forest = as.integer(all_snps %in% z_vars),
  Mahalanobis = as.integer(all_snps %in% maha_vars),
  PCA = as.integer(all_snps %in% pca_vars),
  Univariate = as.integer(all_snps %in% univ_vars)
)

# Rename column 2 to "Isolation Forest" (with a space)
names(upset_data)[2] <- "Isolation Forest"

# Remove the SNP column
upset_data_plot <- upset_data[, -1]

# Create the plot with bars ordered from highest to lowest frequency
UpSetR::upset(
  upset_data_plot,
  sets = c("Isolation Forest", "Mahalanobis", "PCA", "Univariate"),
  order.by = "freq",     # order intersections by frequency (size)
  decreasing = TRUE      # highest first
)

png("upset_plot.png", width = 1200, height = 800, res = 150)

UpSetR::upset(
  upset_data_plot,
  sets = c("Isolation Forest", "Mahalanobis", "PCA", "Univariate"),
  order.by = "freq",
  decreasing = TRUE
)

dev.off()

# Validation with other daasets
setwd("C:\\Users\\Salvatore\\Desktop\\Dottorato\\studio_DL\\Anomaly_Detection\\DROPPED_Validation_Anomaly_detection")
#heatmap_mat <- cbind(
#  variant = rownames(heatmap_mat),
#  as.data.frame(heatmap_mat)
#)

#repeat the following steps for each phenotype
external <- fread("ALS_GRCh37_preprocessed.tsv")
external <- external[,c(1,3)]
external2.0 <- merge(heatmap_mat, external, by ="variant")
external3.0 <- external2.0[,c(1,2,6)]
external3.0$sign_concordance <- ifelse(
  sign(external3.0$z_score_ALS) == sign(external3.0$beta),
  "T",
  "F"
)
table(external3.0$sign_concordance)

write.table(external3.0, file = "SNP_concordi_segno_beta_ALS.txt", quote = F, row.names = F, sep = "\t")
rm(external2.0)
rm(external3.0)
write.table(heatmap_mat, file = "heatmap_mat_info.txt", quote = F, row.names = F, sep = "\t")


# SNPs validated on external datasets
SNPs_validated <- as.data.frame(c("rs11103954","rs137853920","rs141068435","rs145259784","rs146685642","rs147523347","rs149258721","rs17881850","rs182380523","rs190707033","rs199836710","rs20554","rs2230863","rs35277675","rs35381732","rs36220450","rs41276550","rs58973334","rs5952","rs5970","rs61752932","rs61969075","rs9456735","rs971647","rs9893223"))
colnames(SNPs_validated) <- c("variant")

#new heatmap on validated SNPs
heatmap_mat2.0 <- merge(heatmap_mat, SNPs_validated, by ="variant")

heatmap_matrix <- as.matrix(
  heatmap_mat2.0[, c(
    "z_score_ALS",
    "z_score_Alzheimer",
    "z_score_DEMENTIA",
    "z_score_Parkinson"
  )]
)


colnames(heatmap_matrix) <- c("ALS", "Alzheimer", "Dementia", "Parkinson")
rownames(heatmap_matrix) <- heatmap_mat2.0$variant
pheatmap(
  heatmap_matrix,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_row = 0,
  treeheight_col = 0,
  angle_col = 45,
  main = "Validated SNPs",
  fontsize_row = 8,
  filename = "heatmap_validated.png",
  width = 10,
  height = 12
)
                 
