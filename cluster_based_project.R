# START OF SIMULATED DATA SECTION

##
## Linear regression with 10 significant SNPs
##

library(doParallel)
library(foreach)
library(data.table)
library(dplyr)
library(QRank)
library(quantreg)
library(dplyr)
library(ggplot2)
options(scipen = 999)
setwd(getwd())
load("simulation.RData")

# ============================================
# ADDED: Rank-based Inverse Normal Transformation function
# ============================================
rank_INT <- function(x) {
  # Perform rank-based inverse normal transformation
  # x: numeric vector
  n <- length(x)
  # Rank with average ties
  r <- rank(x, ties.method = "average")
  # Convert to quantiles (Blom's approximation)
  q <- (r - 0.375) / (n + 0.25)
  # Inverse normal transform
  transformed <- qnorm(q)
  return(transformed)
}

Num_varianti <- 250 # in origine 250
ult_iterazione = 100 # in origine 100
n_causal_snps <- 10  # Number of causal SNPs to select

# Initialize confusion matrix storage
models <- c("hom", "cauchy", "locscale", "local", "dom")

# Create a list to store all performance data across iterations
all_performance_data <- list()

# ADDED: Option to turn INT on/off
USE_INT <- TRUE  # Set to FALSE to run without normalization

for (j in 1:ult_iterazione) {
  print(j)
  # Reset confusion results for each j iteration
  confusion_results <- lapply(models, function(m) {
    list(TP = 0, FP = 0, FN = 0, TN = 0)
  })
  names(confusion_results) <- models

  # Randomly select columns
  random_cols <- sample(ncol(geno), Num_varianti)
  print(paste("Iteration:", j))
  geno_subset <- geno[1:10000, random_cols] # in origine 10000

  # phenotype simulation function (modified for multiple causal SNPs)
  simu <- function(gt_matrix, covar, model = "hom", beta = rep(0.1, ncol(gt_matrix)), beta.dom = rep(0.2, ncol(gt_matrix)), gamma.cov = 0.5, local.tau = 0.7) {
    n.sample = nrow(gt_matrix)
    base_effect <- gamma.cov * covar

    # Calculate total genetic effect as sum of effects from all causal SNPs
    genetic_effect <- rep(0, n.sample)
    for (s in 1:ncol(gt_matrix)) {
      genetic_effect <- genetic_effect + beta[s] * gt_matrix[, s]
    }

    if (model == "hom") {
      y = base_effect + genetic_effect + rnorm(n.sample)
    } else if (model == "cauchy") {
      y = base_effect + genetic_effect + rcauchy(n.sample)
    } else if (model == "locscale") {
      u = seq(0, 1, length = 50)
      betax = function(u, u0) pmax(0, 5 * (u - u0) / (1 - u0))
      uu = runif(n.sample)
      slopex = betax(uu, local.tau)
      # Apply slope modifier to each causal SNP's effect
      modified_genetic_effect <- rep(0, n.sample)
      for (s in 1:ncol(gt_matrix)) {
        modified_genetic_effect <- modified_genetic_effect + beta[s] * slopex * gt_matrix[, s]
      }
      y = base_effect + modified_genetic_effect + (1 + 0.5 * modified_genetic_effect) * rnorm(n.sample)
    } else if (model == "local") {
      # Step 1: baseline phenotype without genetic effect
      y0 <- base_effect + rnorm(n.sample)

      # Step 2: empirical quantiles of y0
      q_low  <- quantile(y0, probs = 0.6)
      q_high <- quantile(y0, probs = 0.7)

      # Step 3: indicator of local region in the *outcome distribution*
      slopex <- ifelse(y0 > q_low & y0 < q_high, 1, 0)

      # Step 4: apply genetic effect ONLY in that region
      modified_genetic_effect <- rep(0, n.sample)
      for (s in 1:ncol(gt_matrix)) {
        modified_genetic_effect <- modified_genetic_effect +
          beta[s] * slopex * gt_matrix[, s]
      }

      # Step 5: final phenotype
      y <- y0 + modified_genetic_effect

    } else if (model == "dom") {
      # Add dominance effects for each causal SNP
      dominance_effect <- rep(0, n.sample)
      for (s in 1:ncol(gt_matrix)) {
        gt.dom = ifelse(gt_matrix[, s] == 1, 1, 0)
        dominance_effect <- dominance_effect + beta.dom[s] * gt.dom
      }
      y = base_effect + genetic_effect + dominance_effect + rnorm(n.sample)
    }
    return(y)
  }

  # Generate covariate
  COVAR1 <- sample(0:1, nrow(geno_subset), replace = TRUE)

  cat(">> Replicate number:", j, "\n")

  # beta erano 0.2
  beta_values <- rep(0.15, n_causal_snps)
  beta.dom_values <- rep(0.15, n_causal_snps)

  # prima era local.tau = 0.7
  local.tau = 0.85

  n.seed = j
  set.seed(seed = n.seed)

  # Select 5 causal variants (without replacement)
  idx.casual = sample(1:ncol(geno_subset), size = n_causal_snps, replace = FALSE)
  gt_matrix = geno_subset[, idx.casual]
  causal_ids <- colnames(geno_subset)[idx.casual]
  cat("Causal SNPs:", paste(causal_ids, collapse = ", "), "\n")

  # SPEED OPTIMIZATION: Generate all phenotypes at once
  set.seed(n.seed)
  phenotypes <- list(
    y.hom = simu(gt_matrix = gt_matrix, covar = COVAR1, model = "hom", beta = beta_values, beta.dom = beta.dom_values),
    y.cauchy = simu(gt_matrix = gt_matrix, covar = COVAR1, model = "cauchy", beta = beta_values, beta.dom = beta.dom_values),
    y.locscale = simu(gt_matrix = gt_matrix, covar = COVAR1, model = "locscale", beta = beta_values, beta.dom = beta.dom_values, local.tau = local.tau),
    y.local = simu(gt_matrix = gt_matrix, covar = COVAR1, model = "local", beta = beta_values, beta.dom = beta.dom_values, local.tau = local.tau),
    y.dom = simu(gt_matrix = gt_matrix, covar = COVAR1, model = "dom", beta = beta_values, beta.dom = beta.dom_values)
  )

  # ============================================
  # ADDED: Apply INT transformation if USE_INT is TRUE
  # ============================================
  if (USE_INT) {
    cat("Applying rank-based inverse normal transformation (INT) to phenotypes\n")
    phenotypes <- lapply(phenotypes, function(y) {
      # Check for extreme values that might cause issues
      if (any(is.infinite(y))) {
        y[is.infinite(y)] <- sign(y[is.infinite(y)]) * max(abs(y[is.finite(y)]), na.rm = TRUE) * 10
      }
      return(rank_INT(y))
    })
  }

  # Create IID
  protoIID <- seq(1:nrow(geno_subset))
  IID <- paste("HG", protoIID, sep = "")
  geno_current <- cbind(IID, geno_subset)

  # SPEED OPTIMIZATION: Create all phenotype data frames efficiently
  create_pheno_df <- function(y_values, model_name) {
    df <- data.frame(IID = IID, PHENOTYPE = y_values, COVAR1 = COVAR1)
    colnames(df)[2] <- "PHENOTYPE"
    return(df)
  }

  pheno_list <- list(
    y.hom = create_pheno_df(phenotypes$y.hom, "y.hom"),
    y.cauchy = create_pheno_df(phenotypes$y.cauchy, "y.cauchy"),
    y.locscale = create_pheno_df(phenotypes$y.locscale, "y.locscale"),
    y.local = create_pheno_df(phenotypes$y.local, "y.local"),
    y.dom = create_pheno_df(phenotypes$y.dom, "y.dom")
  )

  modelli <- c("y.hom", "y.cauchy", "y.locscale", "y.local", "y.dom")

  geno.mat <- geno_current[,-1]

  # SPEED OPTIMIZATION: Set up parallel cluster once per iteration
  ncores <- detectCores()
  main_cluster <- makeCluster(ncores)
  registerDoParallel(main_cluster)

  # Linear regression function
  linear_regression_test <- function(phe, cov, snp) {
    # Create data frame for regression
    df <- data.frame(
      PHENOTYPE = phe$PHENOTYPE,
      COVAR1 = cov$COVAR1,
      GENOTYPE = snp$GENOTYPE
    )

    # Fit linear regression model
    lm_model <- lm(PHENOTYPE ~ COVAR1 + GENOTYPE, data = df)

    # Extract p-value for GENOTYPE
    summary_lm <- summary(lm_model)
    p_value <- summary_lm$coefficients["GENOTYPE", "Pr(>|t|)"]
    beta <- summary_lm$coefficients["GENOTYPE", "Estimate"]

    return(list(p_value = p_value, beta = beta))
  }

  for (m in modelli) {
    print(paste("Model:", m))
    pheno <- pheno_list[[m]]

    # ADDED: Diagnostic plot to see phenotype distribution (optional)
    # Uncomment to see before/after INT effects
    # if (USE_INT && j == 1) {
    #   png(paste0("phenotype_distribution_", m, "_seed", n.seed, ".png"))
    #   par(mfrow=c(1,2))
    #   hist(pheno$PHENOTYPE, main=paste("Original", m), breaks=50)
    #   qqnorm(pheno$PHENOTYPE, main=paste("QQ-plot", m))
    #   qqline(pheno$PHENOTYPE)
    #   dev.off()
    # }

    # Linear regression analysis
    is.effect.estimated = TRUE

    # SPEED OPTIMIZATION: Use data.table for faster operations in parallel
    clusterExport(main_cluster, c("pheno"), envir = environment())

    # Parallel loop for linear regression
    df.lr.parallel <- foreach(
      idx.geno = 1:ncol(geno.mat),
      .combine = function(...) {
        # SPEED OPTIMIZATION: Use rbindlist for faster combining
        data.table::rbindlist(list(...), fill = TRUE)
      },
      .packages = c("dplyr", "data.table"),
      .errorhandling = "remove"
    ) %dopar% {
      # Create data frame for current SNP
      df.reg <- data.frame(
        PHENOTYPE = pheno$PHENOTYPE,
        COVAR1 = pheno$COVAR1,
        GENOTYPE = geno.mat[, idx.geno],
        check.names = FALSE
      )

      # Perform linear regression test
      result <- linear_regression_test(
        phe = df.reg %>% select(PHENOTYPE),
        cov = df.reg %>% select(COVAR1),
        snp = df.reg %>% select(GENOTYPE)
      )

      # Create result row
      result_row <- data.frame(
        P_LR = result$p_value,
        BETA_LR = result$beta,
        stringsAsFactors = FALSE
      )

      # Return the row
      result_row
    }

    # Convert back to data.frame if needed
    if ("data.table" %in% class(df.lr.parallel)) {
      df.lr <- as.data.frame(df.lr.parallel)
    } else {
      df.lr <- df.lr.parallel
    }

    colnames(df.lr) = c("P_LR", "BETA_LR")
    df.lr = df.lr %>% mutate(ID = colnames(geno.mat), .before = 1)

    # Apply FDR correction to p-values
    df.lr$P_LR_corrected <- p.adjust(df.lr$P_LR, method = "fdr")

    # Update confusion matrix using corrected p-values
    model_name <- sub("y.", "", m)

    significant <- df.lr$P_LR_corrected < 0.05
    is_causal <- df.lr$ID %in% causal_ids

    confusion_results[[model_name]]$TP <- confusion_results[[model_name]]$TP + sum(significant & is_causal, na.rm = TRUE)
    confusion_results[[model_name]]$FP <- confusion_results[[model_name]]$FP + sum(significant & !is_causal, na.rm = TRUE)
    confusion_results[[model_name]]$FN <- confusion_results[[model_name]]$FN + sum(!significant & is_causal, na.rm = TRUE)
    confusion_results[[model_name]]$TN <- confusion_results[[model_name]]$TN + sum(!significant & !is_causal, na.rm = TRUE)

    # Save results
    int_suffix <- ifelse(USE_INT, "_INT", "")
    data.table::fwrite(df.lr, paste0(m, "_linear_regression", int_suffix, "_seed_", n.seed, "_num_varianti_", Num_varianti, ".tsv"),
                       sep = "\t", row.names = FALSE, quote = FALSE)

  }  # End model loop

  # Performance table for this iteration
  perf_table <- data.frame(
    Model = character(),
    TP = integer(),
    FP = integer(),
    FN = integer(),
    TN = integer(),
    Power = numeric(),
    Type_I_Error = numeric(),
    Type_II_Error = numeric(),
    Precision = numeric(),
    F1_Score = numeric(),
    stringsAsFactors = FALSE
  )

  for (model in models) {
    cm <- confusion_results[[model]]

    # Compute core metrics
    power <- cm$TP / (cm$TP + cm$FN + 1e-8)
    Type_I_Error <- cm$FP / (cm$FP + cm$TN + 1e-8)
    Type_II_Error <- cm$FN / (cm$TP + cm$FN + 1e-8)
    precision <- cm$TP / (cm$TP + cm$FP + 1e-8)
    f1_score <- 2 * (precision * power) / (precision + power + 1e-8)

    # Append row to the performance table
    perf_table <- rbind(perf_table, data.frame(
      Model = model,
      TP = cm$TP,
      FP = cm$FP,
      FN = cm$FN,
      TN = cm$TN,
      Power = round(power, 4),
      Type_I_Error = round(Type_I_Error, 4),
      Type_II_Error = round(Type_II_Error, 4),
      Precision = round(precision, 4),
      F1_Score = round(f1_score, 4),
      Iteration = j,
      INT = USE_INT  # ADDED: Track whether INT was used
    ))
  }

  # Store performance data for later aggregation
  all_performance_data[[length(all_performance_data) + 1]] <- perf_table

  # View table
  print(perf_table)

  # Save to CSV or TXT
  int_suffix <- ifelse(USE_INT, "_INT", "")
  filename <- paste0("performance_metrics_linear_regression", int_suffix, "_seed", j, ".txt")
  write.table(perf_table, filename, row.names = FALSE, quote = FALSE, sep = "\t")

  # SPEED OPTIMIZATION: Stop cluster after each iteration
  stopCluster(main_cluster)
  registerDoSEQ()
}  # End iteration loop

# After all iterations are complete, create boxplots
if (length(all_performance_data) > 0) {
  # Combine all performance data
  all_perf_combined <- do.call(rbind, all_performance_data)

  # Convert metrics to percentage for plotting
  all_perf_combined <- all_perf_combined %>%
    mutate(
      Power_pct = Power * 100,
      F1_Score_pct = F1_Score * 100,
      Type_I_Error_pct = Type_I_Error * 100
    )

  # Define consistent color palette for the five models
  model_colors <- c(
    "hom" = "#1f77b4",       # blue
    "cauchy" = "#ff7f0e",    # orange
    "locscale" = "#2ca02c",  # green
    "local" = "#d62728",     # red
    "dom" = "#9467bd"        # purple
  )

  # ======== 1️⃣ POWER BOXPLOT ========
  int_label <- ifelse(USE_INT, "with INT", "without INT")
  p1 <- ggplot(
    all_perf_combined,
    aes(
      x = Model,
      y = Power_pct,
      fill = Model
    )
  ) +
    geom_boxplot(outlier.shape = 19, outlier.size = 1.5, alpha = 0.7) +
    geom_jitter(
      position = position_jitter(width = 0.2),
      size = 0.8,
      alpha = 0.5
    ) +
    scale_fill_manual(values = model_colors) +
    labs(
      title = paste("Distribution of Power (Linear Regression", int_label, ")"),
      x = "Model",
      y = "Power (%)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  int_suffix <- ifelse(USE_INT, "_INT", "")
  ggsave(
    paste0("boxplot_power_linear_regression", int_suffix, ".png"),
    p1,
    width = 10,
    height = 8,
    dpi = 300
  )

  # ======== 2️⃣ F1 SCORE BOXPLOT ========
  p2 <- ggplot(
    all_perf_combined,
    aes(
      x = Model,
      y = F1_Score_pct,
      fill = Model
    )
  ) +
    geom_boxplot(outlier.shape = 19, outlier.size = 1.5, alpha = 0.7) +
    geom_jitter(
      position = position_jitter(width = 0.2),
      size = 0.8,
      alpha = 0.5
    ) +
    scale_fill_manual(values = model_colors) +
    labs(
      title = paste("Distribution of F1 Score (Linear Regression", int_label, ")"),
      x = "Model",
      y = "F1 Score (%)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  ggsave(
    paste0("boxplot_f1score_linear_regression", int_suffix, ".png"),
    p2,
    width = 10,
    height = 8,
    dpi = 300
  )

  # ======== 3️⃣ TYPE I ERROR BOXPLOT ========
  p3 <- ggplot(
    all_perf_combined,
    aes(
      x = Model,
      y = Type_I_Error_pct,
      fill = Model
    )
  ) +
    geom_boxplot(outlier.shape = 19, outlier.size = 1.5, alpha = 0.7) +
    geom_jitter(
      position = position_jitter(width = 0.2),
      size = 0.8,
      alpha = 0.5
    ) +
    scale_fill_manual(values = model_colors) +
    labs(
      title = paste("Distribution of Type I Error (Linear Regression", int_label, ")"),
      x = "Model",
      y = "Type I Error (%)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  ggsave(
    paste0("boxplot_typeIerror_linear_regression", int_suffix, ".png"),
    p3,
    width = 10,
    height = 8,
    dpi = 300
  )

  # Save the combined data for reference
  int_suffix <- ifelse(USE_INT, "_INT", "")
  write.table(all_perf_combined, paste0("all_performance_data_linear_regression", int_suffix, "_combined.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # Print summary statistics
  cat("\n📊 Summary statistics across all iterations (Linear Regression", ifelse(USE_INT, "with INT", "without INT"), "):\n")
  summary_stats <- all_perf_combined %>%
    group_by(Model) %>%
    summarise(
      Mean_Power = mean(Power_pct, na.rm = TRUE),
      SD_Power = sd(Power_pct, na.rm = TRUE),
      Mean_F1 = mean(F1_Score_pct, na.rm = TRUE),
      SD_F1 = sd(F1_Score_pct, na.rm = TRUE),
      Mean_TypeI = mean(Type_I_Error_pct, na.rm = TRUE),
      SD_TypeI = sd(Type_I_Error_pct, na.rm = TRUE),
      .groups = "drop"
    )

  print(summary_stats)

  # Save summary statistics
  write.table(summary_stats, paste0("performance_summary_statistics_linear_regression", int_suffix, ".txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

##
## Quantile regression with 10 significant SNPs


library(data.table)
library(QRank)
library(quantreg)
library(reshape2)
library(viridis)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)
library(doParallel)
library(foreach)
options(scipen = 999)

# ---------------------------
# WORKDIR & DATA
# ---------------------------
setwd(getwd())

# load genotype / simulation data (assumes object 'geno' inside)
load("simulation.RData")

print("dati caricati")

Num_varianti <- 250
n_causal_snps <- 10  # Added: Number of causal SNPs to select
n_replicates <- 100  # Do the analysis X times

# Models
models <- c("hom", "cauchy", "locscale", "local", "dom")

# Create a list to store all performance data across iterations
all_performance_data <- list()  # Added: to store performance data for boxplots

print("entrando nella funzione")
# ---------------------------
# DEFINE THE REPLICATE FUNCTION
# ---------------------------
run_replication <- function(replicate_num, Num_varianti_local = Num_varianti) {
  tryCatch({
    cat(">> Replicate number:", replicate_num, "\n")
    # Set seed
    n.seed <- replicate_num
    print(n.seed)
    set.seed(seed = n.seed)

    # Initialize confusion matrix for this replicate
    replicate_confusion_results <- lapply(models, function(m) list(TP = 0, FP = 0, FN = 0, TN = 0))
    names(replicate_confusion_results) <- models

    # Randomly select columns for this replicate
    random_cols <- sample(ncol(geno), Num_varianti_local)
    geno_subset <- geno[1:10000, random_cols]
    geno_replicate <- geno_subset

    # Generate covariate
    COVAR1 <- sample(0:1, nrow(geno_replicate), replace = TRUE)

    # Causal variants - CHANGED: select 10 causal SNPs instead of 1
    idx.casual <- sample(1:ncol(geno_replicate), size = n_causal_snps, replace = FALSE)
    gt_matrix <- geno_replicate[, idx.casual, drop = FALSE]  # CHANGED: to matrix
    causal_ids <- colnames(geno_replicate)[idx.casual]  # CHANGED: to vector of IDs

    beta_values <- rep(0.15, n_causal_snps)
    beta.dom_values <- rep(0.15, n_causal_snps)

    local.tau = 0.85

    # Simulation function - UPDATED to use multiple causal SNPs
    simu <- function(gt_matrix, covar, model = "hom", beta = rep(0.1, ncol(gt_matrix)), beta.dom = rep(0.2, ncol(gt_matrix)), gamma.cov = 0.5, local.tau = 0.7) {
      n.sample = nrow(gt_matrix)
      base_effect <- gamma.cov * covar

      # Calculate total genetic effect as sum of effects from all causal SNPs
      genetic_effect <- rep(0, n.sample)
      for (s in 1:ncol(gt_matrix)) {
        genetic_effect <- genetic_effect + beta[s] * gt_matrix[, s]
      }

      if (model == "hom") {
        y = base_effect + genetic_effect + rnorm(n.sample)
      } else if (model == "cauchy") {
        y = base_effect + genetic_effect + rcauchy(n.sample)
      } else if (model == "locscale") {
        u = seq(0, 1, length = 50)
        betax = function(u, u0) pmax(0, 5 * (u - u0) / (1 - u0))
        uu = runif(n.sample)
        slopex = betax(uu, local.tau)
        # Apply slope modifier to each causal SNP's effect
        modified_genetic_effect <- rep(0, n.sample)
        for (s in 1:ncol(gt_matrix)) {
          modified_genetic_effect <- modified_genetic_effect + beta[s] * slopex * gt_matrix[, s]
        }
        y = base_effect + modified_genetic_effect + (1 + 0.5 * modified_genetic_effect) * rnorm(n.sample)
      } else if (model == "local") {
        # Step 1: baseline phenotype without genetic effect
        y0 <- base_effect + rnorm(n.sample)

        # Step 2: empirical quantiles of y0
        q_low  <- quantile(y0, probs = 0.6)
        q_high <- quantile(y0, probs = 0.7)

        # Step 3: indicator of local region in the *outcome distribution*
        slopex <- ifelse(y0 > q_low & y0 < q_high, 1, 0)

        # Step 4: apply genetic effect ONLY in that region
        modified_genetic_effect <- rep(0, n.sample)
        for (s in 1:ncol(gt_matrix)) {
          modified_genetic_effect <- modified_genetic_effect +
            beta[s] * slopex * gt_matrix[, s]
        }

        # Step 5: final phenotype
        y <- y0 + modified_genetic_effect

      } else if (model == "dom") {
        # Add dominance effects for each causal SNP
        dominance_effect <- rep(0, n.sample)
        for (s in 1:ncol(gt_matrix)) {
          gt.dom = ifelse(gt_matrix[, s] == 1, 1, 0)
          dominance_effect <- dominance_effect + beta.dom[s] * gt.dom
        }
        y = base_effect + genetic_effect + dominance_effect + rnorm(n.sample)
      }
      return(y)
    }

    # Generate phenotypes - UPDATED to use gt_matrix
    set.seed(n.seed)
    print("genera i fenotipi")
    phenotypes <- list(
      y.hom = simu(gt_matrix = gt_matrix, covar = COVAR1, model = "hom", beta = beta_values, beta.dom = beta.dom_values),
      y.cauchy = simu(gt_matrix = gt_matrix, covar = COVAR1, model = "cauchy", beta = beta_values, beta.dom = beta.dom_values),
      y.locscale = simu(gt_matrix = gt_matrix, covar = COVAR1, model = "locscale", beta = beta_values, beta.dom = beta.dom_values, local.tau = local.tau),
      y.local = simu(gt_matrix = gt_matrix, covar = COVAR1, model = "local", beta = beta_values, beta.dom = beta.dom_values, local.tau = local.tau),
      y.dom = simu(gt_matrix = gt_matrix, covar = COVAR1, model = "dom", beta = beta_values, beta.dom = beta.dom_values)
    )

    protoIID <- seq_len(nrow(geno_replicate))
    IID <- paste0("HG", protoIID)
    geno_replicate <- cbind(IID, geno_replicate)

    # Create phenotype data frames more efficiently
    create_pheno_df <- function(y_values, model_name) {
      df <- data.frame(IID = IID, PHENOTYPE = y_values, COVAR1 = COVAR1)
      colnames(df)[2] <- "PHENOTYPE"
      return(df)
    }

    y.hom <- create_pheno_df(phenotypes$y.hom, "y.hom")
    y.cauchy <- create_pheno_df(phenotypes$y.cauchy, "y.cauchy")
    y.locscale <- create_pheno_df(phenotypes$y.locscale, "y.locscale")
    y.local <- create_pheno_df(phenotypes$y.local, "y.local")
    y.dom <- create_pheno_df(phenotypes$y.dom, "y.dom")

    modelli <- c("y.hom", "y.cauchy", "y.locscale", "y.local", "y.dom")

    # ACAT function with bounds to avoid precision issues (from SCRIPT 1)
    acat <- function(pvals) {
      pvals = pvals[!is.na(pvals)]
      if (length(pvals) == 0) return(NA)
      pvals <- pmax(pmin(pvals, 1 - 1e-15), 1e-15)
      cauchy = tanpi(0.5 - pvals)
      stats = mean(cauchy)
      p = pcauchy(q = stats, lower.tail = FALSE)
      return(p)
    }

    # CORRECTED QRank test function (from SCRIPT 1)
    qrank.test <- function(phe, cov, snp, fit.residual, tau) {
      ltau = length(tau)
      x = as.matrix(snp)
      y = as.matrix(phe)
      n <- nrow(y)
      zz = cbind(rep(1, nrow(y)), as.matrix(cov))

      # SPEED OPTIMIZATION: Precompute VN matrix more efficiently
      VN = outer(tau, tau, function(t1, t2) pmin(t1, t2) - t1 * t2)

      # Get residuals for x
      xstar = lm(x ~ zz - 1)$residual

      SN = NULL
      for(i in 1:ltau) {
        start_idx <- (i-1)*n + 1
        end_idx <- i*n
        resid_vector <- fit.residual$Residual[start_idx:end_idx]

        dual = 1 * (resid_vector > 0)
        ranks = dual - (1 - tau[i])
        Sn = as.matrix(t(xstar) %*% ranks)
        SN = c(SN, Sn)
      }

      # Calculate VN2
      VN2 = VN * as.numeric(t(xstar) %*% xstar)

      # Calculate p-values
      pvalue = pchisq(SN^2/diag(VN2), 1, lower.tail = FALSE)
      pvalue.cauchy = acat(pvalue)
      names(pvalue) = tau

      result = list(tau = tau, quantile.pvalue = pvalue, cauchy.pvalue = pvalue.cauchy)
      return(result)
    }

    # Analysis
    nq <- 9  # Fixed number of bins
    qntl <- seq(1, nq) / (nq + 1)
    cat("Running analysis with", nq, "quantiles for replicate", replicate_num, "\n")

    geno.mat <- geno_replicate[, -1, drop = FALSE]

    # SPEED OPTIMIZATION: Set up parallel cluster
    main_cluster <- makeCluster(n_cores - 1)
    registerDoParallel(main_cluster)

    # SPEED OPTIMIZATION: Pre-export common objects to cluster
    clusterExport(main_cluster, c("qrank.test", "acat"), envir = environment())

    for (m in modelli) {
      cat("  Processing model:", m, "\n")

      pheno <- get(m)
      pheno <- as.data.frame(pheno)
      pheno$PHENOTYPE <- as.numeric(pheno$PHENOTYPE)
      pheno$COVAR1 <- as.numeric(pheno$COVAR1)

      # SPEED OPTIMIZATION: Pre-allocate df.res
      df.res <- data.frame(
        Quantile = rep(qntl, each = nrow(pheno)),
        IID = rep(pheno$IID, times = length(qntl)),
        Residual = numeric(length(qntl) * nrow(pheno))
      )

      fml.null <- "PHENOTYPE ~ COVAR1"
      for (idx.qntl in seq_along(qntl)) {
        fit.null = rq(data = pheno %>% select(PHENOTYPE, COVAR1),
                      formula = as.formula(fml.null),
                      tau = qntl[idx.qntl],
                      method = "fn")

        start_idx <- (idx.qntl - 1) * nrow(pheno) + 1
        end_idx <- idx.qntl * nrow(pheno)
        df.res$Residual[start_idx:end_idx] <- resid(fit.null)
      }

      # Parallel loop for SNPs
      clusterExport(main_cluster, c("pheno", "qntl", "df.res"), envir = environment())

      df.qr.parallel <- foreach(
        idx.geno = 1:ncol(geno.mat),
        .combine = function(...) {
          data.table::rbindlist(list(...), fill = TRUE)
        },
        .packages = c("quantreg", "dplyr", "data.table"),
        .errorhandling = "remove"
      ) %dopar% {
        df.qreg <- data.frame(
          PHENOTYPE = pheno$PHENOTYPE,
          COVAR1 = pheno$COVAR1,
          GENOTYPE = geno.mat[, idx.geno],
          check.names = FALSE
        )

        # Initialize beta vector
        beta <- double()

        # Run quantile regressions
        fml.qreg <- "PHENOTYPE ~ COVAR1 + GENOTYPE"
        for (idx.qntl in seq_along(qntl)) {
          qreg <- rq(
            as.formula(fml.qreg),
            tau = qntl[idx.qntl],
            data = df.qreg,
            method = "fn"
          )
          beta <- c(beta, qreg$coefficients["GENOTYPE"])
        }

        # Perform association test
        p <- qrank.test(
          phe = df.qreg %>% select(PHENOTYPE),
          cov = df.qreg %>% select(COVAR1),
          snp = df.qreg %>% select(GENOTYPE),
          fit.residual = df.res,
          tau = qntl
        )

        # Create result row
        result_row <- data.frame(
          p$cauchy.pvalue,
          t(p$quantile.pvalue),
          t(beta),
          stringsAsFactors = FALSE
        )

        result_row
      }

      # Convert back to data.frame
      if ("data.table" %in% class(df.qr.parallel)) {
        df.qr <- as.data.frame(df.qr.parallel)
      } else {
        df.qr <- df.qr.parallel
      }

      # Set column names
      colnames(df.qr) = c("P_QR",
                          paste0("P_QR", round(qntl, 3)),
                          paste0("BETA_QR", round(qntl, 3)))
      df.qr = df.qr %>% mutate(ID = colnames(geno.mat), .before = 1)

      # Update confusion matrix using corrected P_QR
      model_name <- sub("y.", "", m)
      df.qr$P_QR <- p.adjust(df.qr$P_QR, method = "fdr")
      significant <- df.qr$P_QR < 0.05
      is_causal <- df.qr$ID %in% causal_ids  # Check if ID is in causal_ids

      replicate_confusion_results[[model_name]]$TP <- replicate_confusion_results[[model_name]]$TP + sum(significant & is_causal, na.rm = TRUE)
      replicate_confusion_results[[model_name]]$FP <- replicate_confusion_results[[model_name]]$FP + sum(significant & !is_causal, na.rm = TRUE)
      replicate_confusion_results[[model_name]]$FN <- replicate_confusion_results[[model_name]]$FN + sum(!significant & is_causal, na.rm = TRUE)
      replicate_confusion_results[[model_name]]$TN <- replicate_confusion_results[[model_name]]$TN + sum(!significant & !is_causal, na.rm = TRUE)

      # Save per-replicate results
      fname <- paste0("QR_replicate_", replicate_num, "_seed_", n.seed, "_", m, "_num_quantili_", nq, "_num_varianti_", Num_varianti_local, ".tsv")
      data.table::fwrite(df.qr, fname, sep = "\t", row.names = FALSE, quote = FALSE)

      rm(df.res)
      gc()
    }

    # Stop cluster
    stopCluster(main_cluster)
    registerDoSEQ()

    return(replicate_confusion_results)
  }, error = function(e) {
    message("ERROR in replicate ", replicate_num, " : ", conditionMessage(e))
    # Return zeroed confusion to avoid breaking aggregation
    replicate_confusion_results <- lapply(models, function(m) list(TP = 0, FP = 0, FN = 0, TN = 0))
    names(replicate_confusion_results) <- models
    return(replicate_confusion_results)
  })
}

# ---------------------------
# PARALLEL EXECUTION of iterations
# ---------------------------

# choose number of cores carefully to be polite on shared server
n_cores_available <- detectCores()
n_cores <- as.numeric(Sys.getenv("N_CORES", min(8, n_cores_available - 4)))

print("caricando altre librerie")

cl <- makeCluster(n_cores)

# Export required objects and functions to workers
clusterExport(cl, varlist = c("geno", "models", "Num_varianti", "run_replication", "n_causal_snps"), envir = environment())
clusterEvalQ(cl, {
  library(data.table)
  library(QRank)
  library(quantreg)
  library(dplyr)
  library(tidyr)
  library(doParallel)
  library(foreach)
})

# Run in parallel, each worker executes run_replication(i)
results <- list()
for (i in 1:n_replicates) {
  results[[i]] <- run_replication(i, Num_varianti_local = Num_varianti)
}
print("ternimato il ciclo for")

# ---------------------------
# AGGREGATE results and compute performance metrics
# ---------------------------
overall_confusion_results <- lapply(models, function(m) list(TP = 0, FP = 0, FN = 0, TN = 0))
names(overall_confusion_results) <- models

# Create a list to store performance for each iteration
iteration_performance <- list()
print("stoccando le performance")

for (i in seq_along(results)) {
  replicate_results <- results[[i]]

  # Compute performance metrics for this iteration
  perf_table_iteration <- data.frame()

  for (model in models) {
    cm <- replicate_results[[model]]
    power <- cm$TP / (cm$TP + cm$FN + 1e-8)
    Type_I_Error <- cm$FP / (cm$FP + cm$TN + 1e-8)
    precision <- cm$TP / (cm$TP + cm$FP + 1e-8)
    f1_score <- 2 * (precision * power) / (precision + power + 1e-8)

    perf_table_iteration <- rbind(perf_table_iteration, data.frame(
      Iteration = i,
      Model = model,
      TP = cm$TP,
      FP = cm$FP,
      FN = cm$FN,
      TN = cm$TN,
      Power = round(power, 4),
      Type_I_Error = round(Type_I_Error, 4),
      Precision = round(precision, 4),
      F1_Score = round(f1_score, 4),
      N_quantiles = 9  # Fixed to 9 bins
    ))

    # Also aggregate for overall results
    overall_confusion_results[[model]]$TP <- overall_confusion_results[[model]]$TP + cm$TP
    overall_confusion_results[[model]]$FP <- overall_confusion_results[[model]]$FP + cm$FP
    overall_confusion_results[[model]]$FN <- overall_confusion_results[[model]]$FN + cm$FN
    overall_confusion_results[[model]]$TN <- overall_confusion_results[[model]]$TN + cm$TN
  }

  # Store iteration performance
  iteration_performance[[i]] <- perf_table_iteration
}

# Combine all iteration performance data
all_performance_data <- do.call(rbind, iteration_performance)

# ---------------------------
# Compute final performance table (mean across iterations)
# ---------------------------
perf_table <- data.frame(
  Model = character(),
  Mean_Power = numeric(),
  Mean_Type_I_Error = numeric(),
  Mean_F1_Score = numeric(),
  SD_Power = numeric(),
  SD_Type_I_Error = numeric(),
  SD_F1_Score = numeric(),
  stringsAsFactors = FALSE
)

for (model in models) {
  model_data <- all_performance_data %>% filter(Model == model)

  perf_table <- rbind(perf_table, data.frame(
    Model = model,
    Mean_Power = round(mean(model_data$Power, na.rm = TRUE), 4),
    Mean_Type_I_Error = round(mean(model_data$Type_I_Error, na.rm = TRUE), 4),
    Mean_F1_Score = round(mean(model_data$F1_Score, na.rm = TRUE), 4),
    SD_Power = round(sd(model_data$Power, na.rm = TRUE), 4),
    SD_Type_I_Error = round(sd(model_data$Type_I_Error, na.rm = TRUE), 4),
    SD_F1_Score = round(sd(model_data$F1_Score, na.rm = TRUE), 4)
  ))
}

print(perf_table)
write.table(perf_table, "MEAN_performance_metrics_QR.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# Save the raw performance data
save(results, file = "all_QR_iterations_confusion_results.RData")
write.table(all_performance_data, file = "ALL_QR_iteration_performance_data.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# ---------------------------
# BOXPLOTS
# ---------------------------
# Convert metrics to percentage for plotting
all_performance_data <- all_performance_data %>%
  mutate(
    Power_pct = Power * 100,
    F1_Score_pct = F1_Score * 100,
    Type_I_Error_pct = Type_I_Error * 100
  )

# Define consistent color palette for the five models (from SCRIPT 1)
model_colors <- c(
  "hom" = "#1f77b4",       # blue
  "cauchy" = "#ff7f0e",    # orange
  "locscale" = "#2ca02c",  # green
  "local" = "#d62728",     # red
  "dom" = "#9467bd"        # purple
)

print("iniziano i boxplot")
# ======== 1️⃣ POWER BOXPLOT ========
p1 <- ggplot(
  all_performance_data,
  aes(
    x = Model,
    y = Power_pct,
    fill = Model
  )
) +
  geom_boxplot(outlier.shape = 19, outlier.size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = model_colors) +
  labs(
    title = "Distribution of Power with 9 bins",
    x = "Model",
    y = "Power (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  "boxplot_power_QR.png",
  p1,
  width = 10,
  height = 6,
  dpi = 300
)

# ======== 2️⃣ F1 SCORE BOXPLOT ========
p2 <- ggplot(
  all_performance_data,
  aes(
    x = Model,
    y = F1_Score_pct,
    fill = Model
  )
) +
  geom_boxplot(outlier.shape = 19, outlier.size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = model_colors) +
  labs(
    title = "Distribution of F1 Score with 9 bins",
    x = "Model",
    y = "F1 Score (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  "boxplot_f1score_QR.png",
  p2,
  width = 10,
  height = 6,
  dpi = 300
)

# ======== 3️⃣ TYPE I ERROR BOXPLOT ========
p3 <- ggplot(
  all_performance_data,
  aes(
    x = Model,
    y = Type_I_Error_pct,
    fill = Model
  )
) +
  geom_boxplot(outlier.shape = 19, outlier.size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = model_colors) +
  labs(
    title = "Distribution of Type I Error with 9 bins",
    x = "Model",
    y = "Type I Error (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  "boxplot_typeIerror_QR.png",
  p3,
  width = 10,
  height = 6,
  dpi = 300
)

# ======== 4️⃣ COMBINED PLOT ========
# Convert to long format for combined plot
plot_df_long <- all_performance_data %>%
  select(Model, Power_pct, F1_Score_pct, Type_I_Error_pct) %>%
  pivot_longer(cols = c(Power_pct, F1_Score_pct, Type_I_Error_pct),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(Metric = factor(Metric,
                         levels = c("Power_pct", "F1_Score_pct", "Type_I_Error_pct"),
                         labels = c("Power", "F1 Score", "Type I Error")))

p_combined <- ggplot(
  plot_df_long,
  aes(
    x = Model,
    y = Value,
    fill = Model
  )
) +
  geom_boxplot(outlier.shape = 19, outlier.size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = model_colors) +
  labs(
    title = "Performance Metrics with 9 bins",
    x = "Model",
    y = "Value (%)"
  ) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  "boxplot_combined_QR.png",
  p_combined,
  width = 8,
  height = 12,
  dpi = 300
)

##
## Cluster-based regression with 10 significant SNPs
##
                                    
library(doParallel)
library(foreach)
library(data.table)
library(dplyr)
library(QRank)
library(quantreg)
library(dplyr)
library(ggplot2)
options(scipen = 999)
setwd(getwd())
load("simulation.RData")
Num_varianti <- 250 # in origine 250
ult_iterazione = 100 # in origine 100
n_causal_snps <- 10  # Number of causal SNPs to select

# Initialize confusion matrix storage
models <- c("hom", "cauchy", "locscale", "local", "dom")

# Create a list to store all performance data across iterations
all_performance_data <- list()

for (j in 1:ult_iterazione) {
  # Reset confusion results for each j iteration
  confusion_results <- lapply(models, function(m) {
    list(TP = 0, FP = 0, FN = 0, TN = 0)
  })
  names(confusion_results) <- models

  # Randomly select columns
  random_cols <- sample(ncol(geno), Num_varianti)
  print(paste("Iteration:", j))
  geno_subset <- geno[1:10000, random_cols] # in origine 10000

  # phenotype simulation function (modified for multiple causal SNPs)
  simu <- function(gt_matrix, covar, model = "hom", beta = rep(0.1, ncol(gt_matrix)), beta.dom = rep(0.2, ncol(gt_matrix)), gamma.cov = 0.5, local.tau = 0.7) {
    n.sample = nrow(gt_matrix)
    base_effect <- gamma.cov * covar

    # Calculate total genetic effect as sum of effects from all causal SNPs
    genetic_effect <- rep(0, n.sample)
    for (s in 1:ncol(gt_matrix)) {
      genetic_effect <- genetic_effect + beta[s] * gt_matrix[, s]
    }

    if (model == "hom") {
      y = base_effect + genetic_effect + rnorm(n.sample)
    } else if (model == "cauchy") {
      y = base_effect + genetic_effect + rcauchy(n.sample)
    } else if (model == "locscale") {
      u = seq(0, 1, length = 50)
      betax = function(u, u0) pmax(0, 5 * (u - u0) / (1 - u0))
      uu = runif(n.sample)
      slopex = betax(uu, local.tau)
      # Apply slope modifier to each causal SNP's effect
      modified_genetic_effect <- rep(0, n.sample)
      for (s in 1:ncol(gt_matrix)) {
        modified_genetic_effect <- modified_genetic_effect + beta[s] * slopex * gt_matrix[, s]
      }
      y = base_effect + modified_genetic_effect + (1 + 0.5 * modified_genetic_effect) * rnorm(n.sample)
    } else if (model == "local") {
      # Step 1: baseline phenotype without genetic effect
      y0 <- base_effect + rnorm(n.sample)

      # Step 2: empirical quantiles of y0
      q_low  <- quantile(y0, probs = 0.6)
      q_high <- quantile(y0, probs = 0.7)

      # Step 3: indicator of local region in the *outcome distribution*
      slopex <- ifelse(y0 > q_low & y0 < q_high, 1, 0)

      # Step 4: apply genetic effect ONLY in that region
      modified_genetic_effect <- rep(0, n.sample)
      for (s in 1:ncol(gt_matrix)) {
        modified_genetic_effect <- modified_genetic_effect +
          beta[s] * slopex * gt_matrix[, s]
      }

      # Step 5: final phenotype
      y <- y0 + modified_genetic_effect

    } else if (model == "dom") {
      # Add dominance effects for each causal SNP
      dominance_effect <- rep(0, n.sample)
      for (s in 1:ncol(gt_matrix)) {
        gt.dom = ifelse(gt_matrix[, s] == 1, 1, 0)
        dominance_effect <- dominance_effect + beta.dom[s] * gt.dom
      }
      y = base_effect + genetic_effect + dominance_effect + rnorm(n.sample)
    }
    return(y)
  }

  # Generate covariate
  COVAR1 <- sample(0:1, nrow(geno_subset), replace = TRUE)

  cat(">> Replicate number:", j, "\n")

  # beta erano 0.2
  beta_values <- rep(0.15, n_causal_snps)
  beta.dom_values <- rep(0.15, n_causal_snps)

  # prima era local.tau = 0.7
  local.tau = 0.85

  n.seed = j
  set.seed(seed = n.seed)

  # Select 5 causal variants (without replacement)
  idx.casual = sample(1:ncol(geno_subset), size = n_causal_snps, replace = FALSE)
  gt_matrix = geno_subset[, idx.casual]
  causal_ids <- colnames(geno_subset)[idx.casual]
  cat("Causal SNPs:", paste(causal_ids, collapse = ", "), "\n")

  # SPEED OPTIMIZATION: Generate all phenotypes at once
  set.seed(n.seed)
  phenotypes <- list(
    y.hom = simu(gt_matrix = gt_matrix, covar = COVAR1, model = "hom", beta = beta_values, beta.dom = beta.dom_values),
    y.cauchy = simu(gt_matrix = gt_matrix, covar = COVAR1, model = "cauchy", beta = beta_values, beta.dom = beta.dom_values),
    y.locscale = simu(gt_matrix = gt_matrix, covar = COVAR1, model = "locscale", beta = beta_values, beta.dom = beta.dom_values, local.tau = local.tau),
    y.local = simu(gt_matrix = gt_matrix, covar = COVAR1, model = "local", beta = beta_values, beta.dom = beta.dom_values, local.tau = local.tau),
    y.dom = simu(gt_matrix = gt_matrix, covar = COVAR1, model = "dom", beta = beta_values, beta.dom = beta.dom_values)
  )

  # Create IID
  protoIID <- seq(1:nrow(geno_subset))
  IID <- paste("HG", protoIID, sep = "")
  geno_current <- cbind(IID, geno_subset)

  # SPEED OPTIMIZATION: Create all phenotype data frames efficiently
  create_pheno_df <- function(y_values, model_name) {
    df <- data.frame(IID = IID, PHENOTYPE = y_values, COVAR1 = COVAR1)
    colnames(df)[2] <- "PHENOTYPE"
    return(df)
  }

  pheno_list <- list(
    y.hom = create_pheno_df(phenotypes$y.hom, "y.hom"),
    y.cauchy = create_pheno_df(phenotypes$y.cauchy, "y.cauchy"),
    y.locscale = create_pheno_df(phenotypes$y.locscale, "y.locscale"),
    y.local = create_pheno_df(phenotypes$y.local, "y.local"),
    y.dom = create_pheno_df(phenotypes$y.dom, "y.dom")
  )

  modelli <- c("y.hom", "y.cauchy", "y.locscale", "y.local", "y.dom")

  ### Analysis section
  # ACAT function with bounds to avoid precision issues
  acat <- function(pvals) {
  pvals = pvals[!is.na(pvals)]
  if (length(pvals) == 0) return(NA)
  pvals <- pmax(pmin(pvals, 1 - 1e-15), 1e-15)
  cauchy = tanpi(0.5 - pvals)
  stats = mean(cauchy)
  p = pcauchy(q = stats, lower.tail = FALSE)
  return(p)
  }

  ## CORRECTED QRank test function
  #qrank.test <- function(phe, cov, snp, fit.residual, tau) {
  #ltau = length(tau)
  #x = as.matrix(snp)
  #y = as.matrix(phe)
  #n <- nrow(y)
  #zz = cbind(rep(1, nrow(y)), as.matrix(cov))

  ## MODIFIED QRank test function - using mean(pvalue) instead of ACAT
  qrank.test <- function(phe, cov, snp, fit.residual, tau) {
    ltau = length(tau)
    x = as.matrix(snp)
    y = as.matrix(phe)
    n <- nrow(y)
    zz = cbind(rep(1, nrow(y)), as.matrix(cov))


    # SPEED OPTIMIZATION: Precompute VN matrix more efficiently
    VN = outer(tau, tau, function(t1, t2) pmin(t1, t2) - t1 * t2)

    # Get residuals for x
    xstar = lm(x ~ zz - 1)$residual

    SN = NULL
    for(i in 1:ltau) {
      start_idx <- (i-1)*n + 1
      end_idx <- i*n
      resid_vector <- fit.residual$Residual[start_idx:end_idx]

      dual = 1 * (resid_vector > 0)
      ranks = dual - (1 - tau[i])
      Sn = as.matrix(t(xstar) %*% ranks)
      SN = c(SN, Sn)
    }

    # Calculate VN2
    VN2 = VN * as.numeric(t(xstar) %*% xstar)

    # Calculate p-values
    pvalue = pchisq(SN^2/diag(VN2), 1, lower.tail = FALSE)
    pvalue.cauchy = acat(pvalue)

    # MODIFIED: Use mean(pvalue) instead of ACAT
    pvalue.mean = mean(pvalue, na.rm = TRUE)
    names(pvalue) = tau

    result = list(tau = tau, quantile.pvalue = pvalue, cauchy.pvalue = pvalue.cauchy)
    #result = list(tau = tau, quantile.pvalue = pvalue, cauchy.pvalue = pvalue.mean)
    return(result)
  }
  # in origine era seq(3, 100, by = 1)
  cluster_sizes <- seq(3, 100, by = 1)
  geno.mat <- geno_current[,-1]

  # SPEED OPTIMIZATION: Set up parallel cluster once per iteration
  ncores <- detectCores()
  main_cluster <- makeCluster(ncores)
  registerDoParallel(main_cluster)

  # SPEED OPTIMIZATION: Pre-export common objects to cluster
  clusterExport(main_cluster, c("qrank.test", "acat"), envir = environment())
  #clusterExport(main_cluster, c("qrank.test"), envir = environment())  # Removed "acat"
  clusterEvalQ(main_cluster, {
  })

  for (k in cluster_sizes) {
    # Reset confusion results for each cluster size
    # We need to track performance separately for each cluster size
    confusion_results <- lapply(models, function(m) {
      list(TP = 0, FP = 0, FN = 0, TN = 0)
    })
    names(confusion_results) <- models

    for (m in modelli) {
      print(paste("Cluster:", k, "Model:", m))
      pheno <- pheno_list[[m]]

      # K-means clustering to determine tau levels
      y <- pheno$PHENOTYPE
      kmeans_result <- kmeans(
        y,
        centers = k,
        iter.max = 500,
        nstart = 10,
        algorithm = "MacQueen"
      )

      # Create empirical CDF function
      ecdf_func <- ecdf(y)

      centroids <- sapply(sort(unique(kmeans_result$cluster)), function(cl) {
        mean(y[kmeans_result$cluster == cl])
      })

      centroids_sorted <- sort(centroids)
      tau_levels <- ecdf_func(centroids_sorted)
      tau_levels <- pmax(pmin(tau_levels, 0.99), 0.01)
      qntl <- tau_levels

      # ADD THE DIAGNOSTIC HERE - BEFORE THE FALLBACK CHECK
      cat("K =", k,
          "| Model =", m,
          "| Unique clusters found:", length(unique(kmeans_result$cluster)),
          "| Unique tau levels:", length(unique(tau_levels)),
          "| Fallback triggered:", length(unique(tau_levels)) < 5, "\n")


      # temporanemanete commentato
      # if (length(unique(tau_levels)) < 5) {
      #  tau_levels <- seq(0.1, 0.9, length.out = k)
      # }

      # SPEED OPTIMIZATION: Pre-allocate df.res
      df.res <- data.frame(
        Quantile = rep(tau_levels, each = nrow(pheno)),
        IID = rep(pheno$IID, times = length(tau_levels)),
        Residual = numeric(length(tau_levels) * nrow(pheno))
      )

      # Fit null models using corrected tau levels
      fml.null = "PHENOTYPE ~ COVAR1"
      for (idx.qntl in 1:length(tau_levels)) {
        fit.null = rq(data = pheno %>% select(PHENOTYPE, COVAR1),
                      formula = as.formula(fml.null),
                      tau = tau_levels[idx.qntl],
                      method = "fn")

        start_idx <- (idx.qntl - 1) * nrow(pheno) + 1
        end_idx <- idx.qntl * nrow(pheno)
        df.res$Residual[start_idx:end_idx] <- resid(fit.null)
      }

      # QRank analysis
      is.effect.estimated = TRUE

      # SPEED OPTIMIZATION: Use data.table for faster operations in parallel
      clusterExport(main_cluster, c("pheno", "tau_levels", "df.res"), envir = environment())

      # Parallel loop
      df.qr.parallel <- foreach(
        idx.geno = 1:ncol(geno.mat),
        .combine = function(...) {
          # SPEED OPTIMIZATION: Use rbindlist for faster combining
          data.table::rbindlist(list(...), fill = TRUE)
        },
        .packages = c("quantreg", "dplyr", "data.table"),
        .errorhandling = "remove"
      ) %dopar% {
        # Create data frame for current SNP
        df.qreg <- data.frame(
          PHENOTYPE = pheno$PHENOTYPE,
          COVAR1 = pheno$COVAR1,
          GENOTYPE = geno.mat[, idx.geno],
          check.names = FALSE
        )

        # Initialize beta vector
        beta <- double()

        # Run quantile regressions if needed
        if (is.effect.estimated) {
          fml.qreg <- "PHENOTYPE ~ COVAR1 + GENOTYPE"
          for (idx.qntl in seq_along(tau_levels)) {
            qreg <- rq(
              as.formula(fml.qreg),
              tau = tau_levels[idx.qntl],
              data = df.qreg,
              method = "fn"
            )
            beta <- c(beta, qreg$coefficients["GENOTYPE"])
          }
        }

        # Perform association test
        p <- qrank.test(
          phe = df.qreg %>% select(PHENOTYPE),
          cov = df.qreg %>% select(COVAR1),
          snp = df.qreg %>% select(GENOTYPE),
          fit.residual = df.res,
          tau = tau_levels
        )

        # Create result row
        if (is.effect.estimated) {
          result_row <- data.frame(
            p$cauchy.pvalue,
            t(p$quantile.pvalue),
            t(beta),
            stringsAsFactors = FALSE
          )
        } else {
          result_row <- data.frame(
            p$cauchy.pvalue,
            t(p$quantile.pvalue),
            stringsAsFactors = FALSE
          )
        }

        # Return the row
        result_row
      }

      # Convert back to data.frame if needed
      if ("data.table" %in% class(df.qr.parallel)) {
        df.qr <- as.data.frame(df.qr.parallel)
      } else {
        df.qr <- df.qr.parallel
      }

      # Use sorted centroids for labeling
      qntl_labels_rounded <- round(centroids_sorted, 3)
      colnames(df.qr) = c("P_QR",
                          paste0("P_QR", qntl_labels_rounded),
                          paste0("BETA_QR", qntl_labels_rounded))

      df.qr = df.qr %>% mutate(ID = colnames(geno.mat), .before = 1)

      # CORRECTED: Apply Bonferroni correction ONLY to the Cauchy combined p-value (P_QR)
      # sto provando fdr invece di bonferroni
      df.qr$P_QR <- p.adjust(df.qr$P_QR, method = "fdr")

      # Update confusion matrix using ONLY the corrected P_QR
      model_name <- sub("y.", "", m)

      significant <- df.qr$P_QR < 0.05

      is_causal <- df.qr$ID %in% causal_ids  # Changed to check if ID is in causal_ids

      confusion_results[[model_name]]$TP <- confusion_results[[model_name]]$TP + sum(significant & is_causal, na.rm = TRUE)
      confusion_results[[model_name]]$FP <- confusion_results[[model_name]]$FP + sum(significant & !is_causal, na.rm = TRUE)
      confusion_results[[model_name]]$FN <- confusion_results[[model_name]]$FN + sum(!significant & is_causal, na.rm = TRUE)
      confusion_results[[model_name]]$TN <- confusion_results[[model_name]]$TN + sum(!significant & !is_causal, na.rm = TRUE)

      # SPEED OPTIMIZATION: Use fwrite for faster file writing
      data.table::fwrite(df.qr, paste0(m, "_num_quantili_", length(qntl), "_CR_bonferroni_seed_", n.seed,"_num_varianti_", Num_varianti, ".tsv"),
                         sep = "\t", row.names = FALSE, quote = FALSE)

      # Clean up
      rm(df.res)
      gc()

    }  # End model loop

    # Performance table for this cluster size
    perf_table <- data.frame(
      Model = character(),
      TP = integer(),
      FP = integer(),
      FN = integer(),
      TN = integer(),
      Power = numeric(),
      Type_I_Error = numeric(),
      Type_II_Error = numeric(),
      Precision = numeric(),
      F1_Score = numeric(),
      stringsAsFactors = FALSE
    )

    for (model in models) {
      cm <- confusion_results[[model]]

      # Compute core metrics
      power <- cm$TP / (cm$TP + cm$FN + 1e-8)
      Type_I_Error <- cm$FP / (cm$FP + cm$TN + 1e-8)
      Type_II_Error <- cm$FN / (cm$TP + cm$FN + 1e-8)
      precision <- cm$TP / (cm$TP + cm$FP + 1e-8)
      f1_score <- 2 * (precision * power) / (precision + power + 1e-8)

      # Append row to the performance table
      perf_table <- rbind(perf_table, data.frame(
        Model = model,
        TP = cm$TP,
        FP = cm$FP,
        FN = cm$FN,
        TN = cm$TN,
        Power = round(power, 4),
        Type_I_Error = round(Type_I_Error, 4),
        Type_II_Error = round(Type_II_Error, 4),
        Precision = round(precision, 4),
        F1_Score = round(f1_score, 4),
        N_clusters = k,
        Iteration = j
      ))
    }

    # Store performance data for later aggregation
    all_performance_data[[length(all_performance_data) + 1]] <- perf_table

    # View table
    print(perf_table)

    # Save to CSV or TXT
    filename <- paste0("performance_metrics_seed", j, "_k", k, ".txt")
    write.table(perf_table, filename, row.names = FALSE, quote = FALSE, sep = "\t")

  }  # End cluster loop

  # SPEED OPTIMIZATION: Stop cluster only once per iteration
  stopCluster(main_cluster)
  registerDoSEQ()
}  # End iteration loop

# After all iterations are complete, create boxplots
if (length(all_performance_data) > 0) {
  # Combine all performance data
  all_perf_combined <- do.call(rbind, all_performance_data)

  # Convert metrics to percentage for plotting
  all_perf_combined <- all_perf_combined %>%
    mutate(
      Power_pct = Power * 100,
      F1_Score_pct = F1_Score * 100,
      Type_I_Error_pct = Type_I_Error * 100
    )

  # Convert N_clusters to factor for better boxplot grouping
  all_perf_combined$N_clusters_factor <- factor(all_perf_combined$N_clusters)

  # Define consistent color palette for the five models
  model_colors <- c(
    "hom" = "#1f77b4",       # blue
    "cauchy" = "#ff7f0e",    # orange
    "locscale" = "#2ca02c",  # green
    "local" = "#d62728",     # red
    "dom" = "#9467bd"        # purple
  )


  # ======== 1️⃣ POWER BOXPLOT ========
  p1 <- ggplot(
    all_perf_combined,
    aes(
      x = N_clusters_factor,
      y = Power_pct,
      fill = Model
    )
  ) +
    geom_boxplot(outlier.shape = 19, outlier.size = 1.5, alpha = 0.7) +
    geom_jitter(
      position = position_jitterdodge(jitter.width = 0.2),
      size = 0.8,
      alpha = 0.5
    ) +
    scale_fill_manual(values = model_colors) +
    scale_x_discrete(drop = FALSE) +
    labs(
      title = "Distribution of Power across clusters",
      x = "Number of clusters (k)",
      y = "Power (%)"
    ) +
    facet_wrap(~ Model, ncol = 1, scales = "free_y") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray50"),
      strip.text = element_text(face = "bold")
    )

  ggsave(
    "boxplot_power_vs_k.png",
    p1,
    width = 12,
    height = 18,
    dpi = 300,
    limitsize = FALSE
  )



  # ======== 2️⃣ F1 SCORE BOXPLOT ========
  p2 <- ggplot(
    all_perf_combined,
    aes(
      x = N_clusters_factor,
      y = F1_Score_pct,
      fill = Model
    )
  ) +
    geom_boxplot(outlier.shape = 19, outlier.size = 1.5, alpha = 0.7) +
    geom_jitter(
      position = position_jitterdodge(jitter.width = 0.2),
      size = 0.8,
      alpha = 0.5
    ) +
    scale_fill_manual(values = model_colors) +
    scale_x_discrete(drop = FALSE) +
    labs(
      title = "Distribution of F1 Score across clusters",
      x = "Number of clusters (k)",
      y = "F1 Score (%)"
    ) +
    facet_wrap(~ Model, ncol = 1, scales = "free_y") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.text = element_text(face = "bold")
    )

  ggsave(
    "boxplot_f1score_vs_k.png",
    p2,
    width = 12,
    height = 18,
    dpi = 300,
    limitsize = FALSE
  )

  # ======== 3️⃣ TYPE I ERROR BOXPLOT ========
  p3 <- ggplot(
    all_perf_combined,
    aes(
      x = N_clusters_factor,
      y = Type_I_Error_pct,
      fill = Model
    )
  ) +
    geom_boxplot(outlier.shape = 19, outlier.size = 1.5, alpha = 0.7) +
    geom_jitter(
      position = position_jitterdodge(jitter.width = 0.2),
      size = 0.8,
      alpha = 0.5
    ) +
    scale_fill_manual(values = model_colors) +
    scale_x_discrete(drop = FALSE) +
    labs(
      title = "Distribution of Type I Error across clusters",
      x = "Number of clusters (k)",
      y = "Type I Error (%)"
    ) +
    facet_wrap(~ Model, ncol = 1, scales = "free_y") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.text = element_text(face = "bold")
    )

  ggsave(
    "boxplot_typeIerror_vs_k.png",
    p3,
    width = 12,
    height = 18,
    dpi = 300,
    limitsize = FALSE
  )

  # Save the combined data for reference
  write.table(all_perf_combined, "all_performance_data_combined.txt",
              sep = "\t", row.names = FALSE, quote = FALSE)

  # Print summary statistics
  cat("\n📊 Summary statistics across all iterations:\n")
  summary_stats <- all_perf_combined %>%
    group_by(Model, N_clusters) %>%
    summarise(
      Mean_Power = mean(Power_pct, na.rm = TRUE),
      SD_Power = sd(Power_pct, na.rm = TRUE),
      Mean_F1 = mean(F1_Score_pct, na.rm = TRUE),
      SD_F1 = sd(F1_Score_pct, na.rm = TRUE),
      Mean_TypeI = mean(Type_I_Error_pct, na.rm = TRUE),
      SD_TypeI = sd(Type_I_Error_pct, na.rm = TRUE),
      .groups = "drop"
    )

  print(summary_stats)

  # Save summary statistics
  write.table(summary_stats, "performance_summary_statistics.txt",
              sep = "\t", row.names = FALSE, quote = FALSE)
}
               
# END OF SIMULATED DATA SECTION
# START OF A. thaliana DATA SECTION

               
##
## Linear regression with 9 reference SNPs
##

library(doParallel)
library(foreach)
library(data.table)
library(dplyr)
library(QRank)
library(quantreg)
library(ggplot2)
library(tidyr)


options(scipen = 999)
setwd(getwd())
# ============================================================
# Load real A. thaliana FT10 data
# ============================================================

input <- fread("A_thaliana_TF10_NO_NA.txt")

sample_ids <- input[[1]]
phenotype  <- input[[2]]

geno <- as.matrix(input[, 3:ncol(input)])
colnames(geno) <- colnames(input)[3:ncol(input)]

input <- NULL
gc()

cat("Samples:", nrow(geno), "\n")
cat("SNPs:", ncol(geno), "\n")


# ============================================================
# Reference SNPs to be identified
# ============================================================

reference_snps <- c(
  "Chr5_3188327",
  "Chr4_288749",
  "Chr1_24345355",
  "Chr5_18613760",
  "Chr2_11408771",
  "Chr1_29817195",
  "Chr2_9588628",
  "Chr4_1114351",
  "Chr4_226062"
)

missing <- setdiff(reference_snps, colnames(geno))
if (length(missing) > 0) {
  stop("Missing reference SNPs in genotype matrix: ",
       paste(missing, collapse = ", "))
}

ref_geno <- geno[, reference_snps, drop = FALSE]
is_mono  <- apply(ref_geno, 2, function(x) length(unique(x)) < 2)

if (any(is_mono)) {
  cat("\nMonomorphic reference SNPs (removed):\n")
  for (s in reference_snps[is_mono]) {
    cat("  -", s,
        " (unique value:", unique(ref_geno[, s]), ")\n")
  }
  reference_snps <- reference_snps[!is_mono]
}

n_causal_snps <- length(reference_snps)

cat("\nReference SNPs retained:", n_causal_snps, "\n")
cat(paste0("  ", reference_snps, collapse = "\n"), "\n")


# ============================================================
# Parameters
# ============================================================

Num_varianti   <- 250
N_RANDOM       <- Num_varianti - length(reference_snps)
ult_iterazione <- 100
FDR_THRESHOLD  <- 0.05

# Candidate PC counts to compare
N_PCS_GRID     <- c(4, 5)

USE_INT <- TRUE


# ============================================================
# Rank-based Inverse Normal Transformation
# ============================================================

rank_INT <- function(x) {
  n <- length(x)
  r <- rank(x, ties.method = "average")
  q <- (r - 0.375) / (n + 0.25)
  qnorm(q)
}


# ============================================================
# Apply INT if requested
# ============================================================

phenotype_original <- phenotype

if (USE_INT) {
  cat("Applying rank-based inverse normal transformation (INT) to FT10 phenotype\n")
  if (any(is.infinite(phenotype))) {
    phenotype[is.infinite(phenotype)] <-
      sign(phenotype[is.infinite(phenotype)]) *
      max(abs(phenotype[is.finite(phenotype)]), na.rm = TRUE) * 10
  }
  phenotype <- rank_INT(phenotype)
}


# ============================================================
# PCA (once, genome-wide, EXCLUDING reference SNPs)
# ============================================================

cat("\nRunning PCA on genotype matrix (all SNPs except reference)...\n")

snps_for_pca <- setdiff(colnames(geno), reference_snps)
X_pca        <- geno[, snps_for_pca, drop = FALSE]

cat("  SNPs used for PCA:", ncol(X_pca), "\n")

pca <- prcomp(X_pca, center = TRUE, scale. = FALSE)

## Keep all PCs; subset later per N_PCS value
PCs_all <- pca$x

# ============================================================
# Save covariates (PCs) for downstream analyses (e.g. QR)
# ============================================================

## Full PC matrix (all components)
fwrite(
  data.table(IID = sample_ids, as.data.table(PCs_all)),
  "covariates_PCs_all.tsv",
  sep = "\t"
)

## Convenience copies: first 4 and first 5 PCs (the ones actually used)
fwrite(
  data.table(IID = sample_ids,
             as.data.table(PCs_all[, 1:4, drop = FALSE])),
  "covariates_PCs4.tsv",
  sep = "\t"
)

fwrite(
  data.table(IID = sample_ids,
             as.data.table(PCs_all[, 1:5, drop = FALSE])),
  "covariates_PCs5.tsv",
  sep = "\t"
)

## Also save the phenotype used in the regression (post-INT),
## so QR sees the identical response variable
fwrite(
  data.table(IID = sample_ids, PHENOTYPE = phenotype),
  paste0("phenotype_used", ifelse(USE_INT, "_INT", ""), ".tsv"),
  sep = "\t"
)

cat("Saved covariates: covariates_PCs_all.tsv, covariates_PCs4.tsv, covariates_PCs5.tsv\n")
cat("Saved phenotype:", paste0("phenotype_used", ifelse(USE_INT, "_INT", ""), ".tsv"), "\n")

## Diagnostic table (informational only)
eig      <- pca$sdev^2
prop_var <- eig / sum(eig)

n_show <- min(10, length(eig))
cat("\nVariance explained by first", n_show, "PCs:\n")
print(data.frame(
  PC                = seq_len(n_show),
  Var_explained_pct = round(prop_var[seq_len(n_show)] * 100, 4),
  Cumulative_pct    = round(cumsum(prop_var[seq_len(n_show)]) * 100, 4)
))

rm(X_pca, pca); gc()


# ============================================================
# Non-reference SNP pool
# ============================================================

non_ref_snps <- setdiff(colnames(geno), reference_snps)

## ============================================================
## Save feature sets for reproducibility (independent of seed)
## ============================================================

## 1. Fixed reference SNPs
write.table(
  data.table(SNP = reference_snps),
  "reference_snps_used.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)

## 2. Pre-generate all random draws for all iterations
##    (must be done once, before the loop, so both N_PCS runs
##     use the identical background)
set.seed(123)
random_snps_matrix <- matrix(
  NA_character_,
  nrow = ult_iterazione,
  ncol = N_RANDOM,
  dimnames = list(paste0("Iter", seq_len(ult_iterazione)),
                  paste0("SNP_", seq_len(N_RANDOM)))
)
for (j in seq_len(ult_iterazione)) {
  set.seed(j)
  random_snps_matrix[j, ] <- sample(non_ref_snps, N_RANDOM, replace = FALSE)
}

## 3. Save in two formats:
##    - wide (one row per iteration, 241 columns)
##    - long (one row per iteration-SNP pair, with a type column)
fwrite(
  data.table(Iteration = seq_len(ult_iterazione),
             as.data.table(random_snps_matrix)),
  "random_snps_per_iteration.tsv",
  sep = "\t"
)

features_long <- rbindlist(lapply(seq_len(ult_iterazione), function(j) {
  rbind(
    data.table(Iteration = j, SNP = reference_snps,   Type = "reference"),
    data.table(Iteration = j, SNP = random_snps_matrix[j, ], Type = "random")
  )
}))
fwrite(features_long, "features_per_iteration.tsv", sep = "\t")

cat("Non-reference SNPs available for random sampling:",
    length(non_ref_snps), "\n")


# ============================================================
# Linear regression helper
# ============================================================

linear_regression_test <- function(phe, cov, snp) {
  
  df <- data.frame(
    PHENOTYPE = phe$PHENOTYPE,
    cov,
    GENOTYPE = snp$GENOTYPE,
    check.names = FALSE
  )
  
  pc_names     <- colnames(cov)
  formula_text <- paste(
    "PHENOTYPE ~",
    paste(c(pc_names, "GENOTYPE"), collapse = " + ")
  )
  
  lm_model   <- lm(as.formula(formula_text), data = df)
  summary_lm <- summary(lm_model)
  
  list(
    p_value = summary_lm$coefficients["GENOTYPE", "Pr(>|t|)"],
    beta    = summary_lm$coefficients["GENOTYPE", "Estimate"]
  )
}


# ============================================================
# Parallel cluster (setup once)
# ============================================================

ncores       <- detectCores()
main_cluster <- makeCluster(ncores)
registerDoParallel(main_cluster)


# ============================================================
# Storage for the comparison
# ============================================================

comparison_summary <- list()
comparison_perf    <- list()


# ============================================================
# OUTER LOOP: N_PCS values to compare
# ============================================================

for (N_PCS in N_PCS_GRID) {
  
  cat("\n\n############################################################\n")
  cat("RUNNING WITH N_PCS =", N_PCS, "\n")
  cat("############################################################\n")
  
  ## Build pheno with first N_PCS PCs
  PCs   <- PCs_all[, seq_len(N_PCS), drop = FALSE]
  pheno <- data.frame(
    IID       = sample_ids,
    PHENOTYPE = phenotype,
    PCs,
    check.names = FALSE
  )
  
  int_suffix <- ifelse(USE_INT, "_INT", "")
  pcs_suffix <- paste0("_PCs", N_PCS)
  run_tag    <- paste0(int_suffix, pcs_suffix)
  
  ## ---- Inner loop: 100 iterations ----
  all_performance_data <- list()
  all_results_data     <- list()
  
  for (j in seq_len(ult_iterazione)) {
    
    cat("\n----- [N_PCS =", N_PCS, "] Iteration", j, "/", ult_iterazione, "-----\n")
    set.seed(j)
    
    ## 1. Draw 241 random SNPs
    random_snps <- random_snps_matrix[j, ]
    
    ## 2. Build 250-SNP panel (9 reference + 241 random)
    geno.mat <- cbind(
      geno[, reference_snps, drop = FALSE],
      geno[, random_snps,    drop = FALSE]
    )
    colnames(geno.mat) <- c(reference_snps, random_snps)
    
    ## 3. Export to workers
    clusterExport(
      main_cluster,
      c("pheno", "linear_regression_test", "geno.mat", "N_PCS"),
      envir = environment()
    )
    
    ## 4. Per-SNP regression
    df.lr.parallel <- foreach(
      idx.geno = 1:ncol(geno.mat),
      .combine = function(...) data.table::rbindlist(list(...), fill = TRUE),
      .packages = c("dplyr", "data.table"),
      .errorhandling = "remove"
    ) %dopar% {
      
      snp <- data.frame(GENOTYPE = geno.mat[, idx.geno])
      cov <- pheno[, paste0("PC", 1:N_PCS), drop = FALSE]
      
      result <- linear_regression_test(
        phe = pheno[, "PHENOTYPE", drop = FALSE],
        cov = cov,
        snp = snp
      )
      
      data.frame(
        P_LR    = result$p_value,
        BETA_LR = result$beta,
        stringsAsFactors = FALSE
      )
    }
    
    ## 5. Convert and label
    if ("data.table" %in% class(df.lr.parallel)) {
      df.lr <- as.data.frame(df.lr.parallel)
    } else {
      df.lr <- df.lr.parallel
    }
    
    colnames(df.lr) <- c("P_LR", "BETA_LR")
    df.lr <- df.lr %>% mutate(ID = colnames(geno.mat), .before = 1)
    
    ## 6. FDR
    df.lr$P_LR_corrected           <- p.adjust(df.lr$P_LR, method = "fdr")
    df.lr$P_LR_nominal_significant <- df.lr$P_LR < 0.05
    df.lr$FDR_significant          <- df.lr$P_LR_corrected < FDR_THRESHOLD
    df.lr$Iteration                <- j
    df.lr$N_PCS                    <- N_PCS
    
    ## 7. Confusion matrix
    significant  <- df.lr$FDR_significant
    is_reference <- df.lr$ID %in% reference_snps
    
    cm <- list(
      TP = sum( significant &  is_reference, na.rm = TRUE),
      FP = sum( significant & !is_reference, na.rm = TRUE),
      FN = sum(!significant &  is_reference, na.rm = TRUE),
      TN = sum(!significant & !is_reference, na.rm = TRUE)
    )
    
    ## 8. Metrics
    power     <- cm$TP / (cm$TP + cm$FN + 1e-8)
    typeI     <- cm$FP / (cm$FP + cm$TN + 1e-8)
    typeII    <- cm$FN / (cm$TP + cm$FN + 1e-8)
    precision <- cm$TP / (cm$TP + cm$FP + 1e-8)
    f1        <- 2 * precision * power / (precision + power + 1e-8)
    
    n_ref_testable <- sum(reference_snps %in% df.lr$ID)
    power_nominal  <- sum(df.lr$P_LR_nominal_significant &
                            df.lr$ID %in% reference_snps) /
      (n_ref_testable + 1e-8)
    
    perf_table <- data.frame(
      Model         = "FT10",
      N_PCS         = N_PCS,
      TP = cm$TP, FP = cm$FP, FN = cm$FN, TN = cm$TN,
      Power         = round(power,         4),
      Power_nominal = round(power_nominal, 4),
      Type_I_Error  = round(typeI,         4),
      Type_II_Error = round(typeII,        4),
      Precision     = round(precision,     4),
      F1_Score      = round(f1,            4),
      Iteration     = j,
      INT           = USE_INT
    )
    
    all_performance_data[[j]] <- perf_table
    all_results_data[[j]]     <- df.lr
    
    print(perf_table)
    
    ## 9. Save per-iteration outputs
    data.table::fwrite(
      df.lr,
      paste0("FT10_linear_regression", run_tag,
             "_seed_", j, "_num_varianti_", Num_varianti, ".tsv"),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
    
    reference_results <- df.lr %>%
      filter(ID %in% reference_snps) %>%
      arrange(P_LR)
    
    data.table::fwrite(
      reference_results,
      paste0("FT10_reference_SNP_results", run_tag,
             "_seed_", j, ".tsv"),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
  }
  
  ## ---- Aggregate for this N_PCS ----
  perf_combined    <- do.call(rbind, all_performance_data)
  results_combined <- do.call(rbind, all_results_data)
  
  data.table::fwrite(
    perf_combined,
    paste0("all_performance_data_linear_regression", run_tag, "_combined.txt"),
    sep = "\t"
  )
  data.table::fwrite(
    results_combined,
    paste0("all_per_snp_results_linear_regression", run_tag, "_combined.txt"),
    sep = "\t"
  )
  
  ## ---- Summary for this N_PCS ----
  ## NOTE: N_PCS is NOT referenced inside summarise() because the column
  ## `N_PCS` already exists in perf_combined (100 rows) and would clash
  ## with the outer scalar. We compute the metrics first, then attach the
  ## scalar columns.
  summary_stats <- perf_combined %>%
    summarise(
      Mean_Power     = mean(Power,          na.rm = TRUE),
      SD_Power       = sd(Power,            na.rm = TRUE),
      Mean_Power_nom = mean(Power_nominal,  na.rm = TRUE),
      SD_Power_nom   = sd(Power_nominal,    na.rm = TRUE),
      Mean_F1        = mean(F1_Score,       na.rm = TRUE),
      SD_F1          = sd(F1_Score,         na.rm = TRUE),
      Mean_TypeI     = mean(Type_I_Error,   na.rm = TRUE),
      SD_TypeI       = sd(Type_I_Error,     na.rm = TRUE),
      Mean_TypeII    = mean(Type_II_Error,  na.rm = TRUE),
      Mean_Precision = mean(Precision,      na.rm = TRUE)
    )
  summary_stats$N_PCS       <- N_PCS
  summary_stats$N_reference <- n_causal_snps
  
  ## Reorder so N_PCS comes first
  summary_stats <- summary_stats %>%
    select(N_PCS, Mean_Power, SD_Power, Mean_Power_nom, SD_Power_nom,
           Mean_F1, SD_F1, Mean_TypeI, SD_TypeI, Mean_TypeII,
           Mean_Precision, N_reference)
  
  cat("\n===== Summary for N_PCS =", N_PCS, "=====\n")
  print(summary_stats)
  
  write.table(
    summary_stats,
    paste0("performance_summary_statistics_linear_regression", run_tag, ".txt"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  
  comparison_summary[[as.character(N_PCS)]] <- summary_stats
  comparison_perf[[as.character(N_PCS)]]    <- perf_combined
}


# ============================================================
# Stop cluster
# ============================================================

stopCluster(main_cluster)
registerDoSEQ()


# ============================================================
# Comparison across N_PCS values
# ============================================================

cat("\n\n############################################################\n")
cat("COMPARISON OF N_PCS VALUES\n")
cat("############################################################\n")

comparison_df <- do.call(rbind, comparison_summary)
print(comparison_df, row.names = FALSE)

write.table(
  comparison_df,
  paste0("N_PCS_comparison_summary", int_suffix, ".txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

## Decide the winner by a combined criterion
best_row <- comparison_df %>%
  mutate(score = Mean_Power + Mean_Power_nom - Mean_TypeI) %>%
  arrange(desc(score)) %>%
  slice(1)

cat("\n---- Best N_PCS by combined criterion (Power + Power_nom - TypeI) ----\n")
cat("N_PCS =", best_row$N_PCS, "\n")
cat("  Mean_Power       :", round(best_row$Mean_Power,     4), "\n")
cat("  Mean_Power_nom   :", round(best_row$Mean_Power_nom, 4), "\n")
cat("  Mean_F1          :", round(best_row$Mean_F1,        4), "\n")
cat("  Mean_TypeI_Error :", round(best_row$Mean_TypeI,     4), "\n")


# ============================================================
# Combined boxplots across both N_PCS values
# ============================================================

all_combined <- do.call(rbind, comparison_perf)

all_combined <- all_combined %>%
  mutate(
    Power_pct         = Power * 100,
    Power_nominal_pct = Power_nominal * 100,
    F1_Score_pct      = F1_Score * 100,
    Type_I_Error_pct  = Type_I_Error * 100,
    N_PCS_f           = factor(paste0("N_PCS = ", N_PCS),
                               levels = paste0("N_PCS = ", N_PCS_GRID))
  )

int_label    <- ifelse(USE_INT, "with INT", "without INT")
model_colors <- c("N_PCS = 4" = "#1f77b4",
                  "N_PCS = 5" = "#ff7f0e")


# ------------------------------------------------------------
# Power boxplot
# ------------------------------------------------------------
p1 <- ggplot(all_combined,
             aes(x = N_PCS_f, y = Power_pct, fill = N_PCS_f)) +
  geom_boxplot(outlier.shape = 19, outlier.size = 1.5, alpha = 0.7) +
  geom_jitter(position = position_jitter(width = 0.2),
              size = 0.8, alpha = 0.5) +
  scale_fill_manual(values = model_colors) +
  labs(title = paste("Power (FDR) — Linear Regression", int_label),
       x = NULL, y = "Power (%)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(paste0("boxplot_power_comparison", int_suffix, ".png"),
       p1, width = 8, height = 6, dpi = 300)


# ------------------------------------------------------------
# Nominal power boxplot
# ------------------------------------------------------------
p1b <- ggplot(all_combined,
              aes(x = N_PCS_f, y = Power_nominal_pct, fill = N_PCS_f)) +
  geom_boxplot(outlier.shape = 19, outlier.size = 1.5, alpha = 0.7) +
  geom_jitter(position = position_jitter(width = 0.2),
              size = 0.8, alpha = 0.5) +
  scale_fill_manual(values = model_colors) +
  labs(title = paste("Nominal Power (P<0.05) — Linear Regression", int_label),
       x = NULL, y = "Nominal Power (%)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(paste0("boxplot_power_nominal_comparison", int_suffix, ".png"),
       p1b, width = 8, height = 6, dpi = 300)


# ------------------------------------------------------------
# F1 score boxplot
# ------------------------------------------------------------
p2 <- ggplot(all_combined,
             aes(x = N_PCS_f, y = F1_Score_pct, fill = N_PCS_f)) +
  geom_boxplot(outlier.shape = 19, outlier.size = 1.5, alpha = 0.7) +
  geom_jitter(position = position_jitter(width = 0.2),
              size = 0.8, alpha = 0.5) +
  scale_fill_manual(values = model_colors) +
  labs(title = paste("F1 Score — Linear Regression", int_label),
       x = NULL, y = "F1 Score (%)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(paste0("boxplot_f1score_comparison", int_suffix, ".png"),
       p2, width = 8, height = 6, dpi = 300)


# ------------------------------------------------------------
# Type I error boxplot
# ------------------------------------------------------------
p3 <- ggplot(all_combined,
             aes(x = N_PCS_f, y = Type_I_Error_pct, fill = N_PCS_f)) +
  geom_boxplot(outlier.shape = 19, outlier.size = 1.5, alpha = 0.7) +
  geom_jitter(position = position_jitter(width = 0.2),
              size = 0.8, alpha = 0.5) +
  scale_fill_manual(values = model_colors) +
  labs(title = paste("Empirical FP proportion — Linear Regression", int_label),
       x = NULL, y = "Type I Error (%)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(paste0("boxplot_typeIerror_comparison", int_suffix, ".png"),
       p3, width = 8, height = 6, dpi = 300)



# ============================================================
# Combined boxplots: Power, F1, Type I Error in a single image
# ============================================================

plot_three <- all_combined %>%
  select(N_PCS_f, Power_pct, F1_Score_pct, Type_I_Error_pct) %>%
  pivot_longer(
    cols = c(Power_pct, F1_Score_pct, Type_I_Error_pct),
    names_to = "Metric",
    values_to = "Value_pct"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("Power_pct", "F1_Score_pct", "Type_I_Error_pct"),
      labels = c("Power (FDR)", "F1 Score", "Type I Error")
    )
  )

p_combined <- ggplot(plot_three,
                     aes(x = N_PCS_f, y = Value_pct, fill = N_PCS_f)) +
  geom_boxplot(outlier.shape = 19, outlier.size = 1.2,
               alpha = 0.7, width = 0.6) +
  geom_jitter(position = position_jitter(width = 0.15),
              size = 0.6, alpha = 0.45) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = model_colors) +
  labs(
    title = paste("Linear Regression performance —", int_label),
    x    = NULL,
    y    = "%"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    strip.text      = element_text(face = "bold", size = 12)
  )

ggsave(
  paste0("boxplot_power_F1_TypeI_combined", int_suffix, ".png"),
  p_combined,
  width  = 12,
  height = 5,
  dpi    = 300
)

# ============================================================
# Boxplots of Power, F1 and Type I Error — N_PCS = 4 only
# ============================================================

## Keep only the N_PCS = 4 rows
all_pcs4 <- all_combined %>%
  filter(N_PCS == 4)

## Reshape to long format for faceting
plot_pcs4 <- all_pcs4 %>%
  select(Power_pct, F1_Score_pct, Type_I_Error_pct) %>%
  pivot_longer(
    cols      = c(Power_pct, F1_Score_pct, Type_I_Error_pct),
    names_to  = "Metric",
    values_to = "Value_pct"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("Power_pct", "F1_Score_pct", "Type_I_Error_pct"),
      labels = c("Power (FDR)", "F1 Score", "Type I Error")
    )
  )

## Build the combined plot
p_pcs4 <- ggplot(plot_pcs4,
                 aes(x = Metric, y = Value_pct, fill = Metric)) +
  geom_boxplot(outlier.shape = 19, outlier.size = 1.2,
               alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.12, size = 0.6, alpha = 0.45) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(
    "Power (FDR)"  = "#1f77b4",
    "F1 Score"     = "#2ca02c",
    "Type I Error" = "#d62728"
  )) +
  labs(
    title = paste("Linear Regression"),
    x     = NULL,
    y     = "%"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    axis.text.x     = element_blank(),
    strip.text      = element_text(face = "bold", size = 12)
  )

## Save
ggsave(
  paste0("boxplot_power_F1_TypeI_PCs4", int_suffix, ".png"),
  p_pcs4,
  width  = 11,
  height = 5,
  dpi    = 300
)

cat("\nSaved: boxplot_power_F1_TypeI_PCs4", int_suffix, ".png\n", sep = "")

# ============================================================
# Save combined data
# ============================================================

write.table(
  all_combined,
  paste0("all_performance_data_linear_regression_combined", int_suffix, ".txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)


# ============================================================
# Final report
# ============================================================

cat("\n============================================================\n")
cat("Analysis completed.\n")
cat("Reference SNPs analyzed:", n_causal_snps, "\n")
cat("Total SNPs per iteration:", Num_varianti, "\n")
cat("Random SNPs per iteration:", N_RANDOM, "\n")
cat("Samples:", nrow(geno), "\n")
cat("Iterations:", ult_iterazione, "\n")
cat("N_PCS values compared:", paste(N_PCS_GRID, collapse = ", "), "\n")
cat("FDR threshold:", FDR_THRESHOLD, "\n")
cat("INT:", USE_INT, "\n")
cat("============================================================\n")

res <- fread("all_per_snp_results_linear_regression_INT_PCs4_combined.txt")
res[Iteration == 1 & ID %in% reference_snps,
    .(ID, P_LR, P_LR_corrected, BETA_LR, FDR_significant)][order(P_LR)]



##
## Quantile regression with 9 reference SNPs
##


library(data.table)
library(QRank)
library(quantreg)
library(reshape2)
library(viridis)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)
library(doParallel)
library(foreach)
options(scipen = 999)

setwd(getwd())
# ============================================================
# 1. Load data
# ============================================================
cat("Loading genotype data...\n")
input      <- fread("A_thaliana_TF10_NO_NA.txt")
sample_ids <- input[[1]]
geno       <- as.matrix(input[, 3:ncol(input)])
colnames(geno) <- colnames(input)[3:ncol(input)]
input <- NULL; gc()

cat("Samples:", nrow(geno), "  SNPs:", ncol(geno), "\n")

## Reference SNPs (9, fixed)
reference_snps <- read.delim("reference_snps_used.tsv")$SNP
n_ref <- length(reference_snps)
cat("Reference SNPs:", n_ref, "\n")

## Random SNPs matrix (100 x 241), pre-generated
random_snps_wide <- fread("random_snps_per_iteration.tsv")

## Covariates: IID + PC1..PC4
covariates <- fread("covariates_PCs4.tsv")

## Phenotype — toggle between INT and raw
USE_INT_FOR_QR <- FALSE   # TRUE = INT phenotype, FALSE = original untransformed

if (USE_INT_FOR_QR) {
  pheno_used  <- fread("phenotype_used_INT.tsv")
  pheno_label <- "INT"
} else {
  raw_input   <- fread("A_thaliana_TF10_NO_NA.txt", select = c(1, 2))
  pheno_used  <- data.frame(IID       = raw_input[[1]],
                            PHENOTYPE = raw_input[[2]])
  pheno_label <- "raw"
  rm(raw_input); gc()
}
cat("Phenotype used:", pheno_label, "\n")

match_cov <- match(sample_ids, covariates$IID)
match_phe <- match(sample_ids, pheno_used$IID)
if (any(is.na(match_cov)) || any(is.na(match_phe))) {
  stop("Sample IDs do not match between files.")
}
covariates <- covariates[match_cov, ]
pheno_used <- pheno_used[match_phe, ]

df_analysis <- data.frame(
  IID       = sample_ids,
  PHENOTYPE = pheno_used$PHENOTYPE,
  covariates[, paste0("PC", 1:4), with = FALSE]
)

cat("Covariate + phenotype alignment done.\n")

# ============================================================
# 2. Parameters
# ============================================================
Num_varianti  <- 250
N_RANDOM      <- ncol(random_snps_wide) - 1
n_replicates  <- nrow(random_snps_wide)
nq            <- 9
qntl          <- seq(1, nq) / (nq + 1)
FDR_THRESHOLD <- 0.05
N_PCS         <- 4

# ============================================================
# 3. ACAT
# ============================================================
acat <- function(pvals) {
  pvals = pvals[!is.na(pvals)]
  if (length(pvals) == 0) return(NA)
  pvals <- pmax(pmin(pvals, 1 - 1e-15), 1e-15)
  cauchy = tanpi(0.5 - pvals)
  stats  = mean(cauchy)
  pcauchy(q = stats, lower.tail = FALSE)
}

# ============================================================
# 4. QRank test
# ============================================================
qrank.test <- function(phe, cov, snp, fit.residual, tau) {
  ltau = length(tau)
  x = as.matrix(snp); y = as.matrix(phe); n <- nrow(y)
  zz = cbind(rep(1, n), as.matrix(cov))
  VN = outer(tau, tau, function(t1, t2) pmin(t1, t2) - t1 * t2)
  xstar = lm(x ~ zz - 1)$residual
  SN = NULL
  for (i in 1:ltau) {
    idx <- ((i - 1) * n + 1):(i * n)
    ranks = 1 * (fit.residual$Residual[idx] > 0) - (1 - tau[i])
    SN = c(SN, as.matrix(t(xstar) %*% ranks))
  }
  VN2 = VN * as.numeric(t(xstar) %*% xstar)
  pvalue = pchisq(SN^2 / diag(VN2), 1, lower.tail = FALSE)
  names(pvalue) = tau
  list(tau = tau, quantile.pvalue = pvalue,
       cauchy.pvalue = acat(pvalue))
}

# ============================================================
# 5. Null QR residuals
# ============================================================
cat("Computing null QR residuals...\n")
df.res <- data.frame(
  Quantile = rep(qntl, each = nrow(df_analysis)),
  IID      = rep(df_analysis$IID, times = length(qntl)),
  Residual = numeric(length(qntl) * nrow(df_analysis))
)
fml.null <- paste("PHENOTYPE ~",
                  paste(paste0("PC", 1:N_PCS), collapse = " + "))
for (idx.qntl in seq_along(qntl)) {
  fit.null = rq(as.formula(fml.null), tau = qntl[idx.qntl],
                data = df_analysis, method = "fn")
  idx <- ((idx.qntl - 1) * nrow(df_analysis) + 1):(idx.qntl * nrow(df_analysis))
  df.res$Residual[idx] <- resid(fit.null)
}
cat("Null residuals computed.\n")

# ============================================================
# 6. Parallel cluster
# ============================================================
n_cores_available <- detectCores()
n_cores <- as.numeric(Sys.getenv("N_CORES", min(8, n_cores_available - 4)))
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterExport(cl,
              c("df_analysis", "geno", "qntl", "df.res",
                "qrank.test", "acat", "reference_snps", "N_PCS"),
              envir = environment())
clusterEvalQ(cl, {
  library(quantreg); library(dplyr); library(data.table)
})

# ============================================================
# 7. Main loop  <-- per-iteration writes now INSIDE the loop
# ============================================================
all_performance_data <- list()
all_results_data     <- list()

for (j in seq_len(n_replicates)) {
  cat("\n===== Iteration", j, "/", n_replicates, "=====\n")
  
  random_snps_j <- as.character(random_snps_wide[j, -1, with = FALSE])
  
  geno.mat <- cbind(
    geno[, reference_snps, drop = FALSE],
    geno[, random_snps_j,  drop = FALSE]
  )
  colnames(geno.mat) <- c(reference_snps, random_snps_j)
  
  clusterExport(cl, "geno.mat", envir = environment())
  
  df.qr.parallel <- foreach(
    idx.geno = 1:ncol(geno.mat),
    .combine = function(...) data.table::rbindlist(list(...), fill = TRUE),
    .packages = c("quantreg", "dplyr", "data.table"),
    .errorhandling = "remove"
  ) %dopar% {
    df.qreg <- cbind(
      data.frame(PHENOTYPE = df_analysis$PHENOTYPE),
      df_analysis[, paste0("PC", 1:N_PCS)],
      data.frame(GENOTYPE = geno.mat[, idx.geno])
    )
    fml.qreg <- paste("PHENOTYPE ~",
                      paste(c(paste0("PC", 1:N_PCS), "GENOTYPE"),
                            collapse = " + "))
    beta <- double()
    for (idx.qntl in seq_along(qntl)) {
      qreg <- rq(as.formula(fml.qreg), tau = qntl[idx.qntl],
                 data = df.qreg, method = "fn")
      beta <- c(beta, qreg$coefficients["GENOTYPE"])
    }
    p <- qrank.test(
      phe = df.qreg %>% select(PHENOTYPE),
      cov = df.qreg %>% select(all_of(paste0("PC", 1:N_PCS))),
      snp = df.qreg %>% select(GENOTYPE),
      fit.residual = df.res, tau = qntl
    )
    data.frame(p$cauchy.pvalue, t(p$quantile.pvalue), t(beta),
               stringsAsFactors = FALSE)
  }
  
  if ("data.table" %in% class(df.qr.parallel)) {
    df.qr <- as.data.frame(df.qr.parallel)
  } else {
    df.qr <- df.qr.parallel
  }
  
  colnames(df.qr) <- c("P_QR",
                       paste0("P_QR",    round(qntl, 3)),
                       paste0("BETA_QR", round(qntl, 3)))
  df.qr <- df.qr %>% mutate(ID = colnames(geno.mat), .before = 1)
  df.qr$Iteration       <- j
  df.qr$P_QR_corrected  <- p.adjust(df.qr$P_QR, method = "fdr")
  df.qr$FDR_significant <- df.qr$P_QR_corrected < FDR_THRESHOLD
  df.qr$is_reference    <- df.qr$ID %in% reference_snps
  
  significant  <- df.qr$FDR_significant
  is_reference <- df.qr$is_reference
  TP <- sum( significant &  is_reference, na.rm = TRUE)
  FP <- sum( significant & !is_reference, na.rm = TRUE)
  FN <- sum(!significant &  is_reference, na.rm = TRUE)
  TN <- sum(!significant & !is_reference, na.rm = TRUE)
  power     <- TP / (TP + FN + 1e-8)
  typeI     <- FP / (FP + TN + 1e-8)
  typeII    <- FN / (TP + FN + 1e-8)
  precision <- TP / (TP + FP + 1e-8)
  f1        <- 2 * precision * power / (precision + power + 1e-8)
  n_ref_testable <- sum(reference_snps %in% df.qr$ID)
  power_nominal  <- sum(df.qr$P_QR < 0.05 & df.qr$is_reference) /
    (n_ref_testable + 1e-8)
  
  perf_table <- data.frame(
    Model = "FT10", Iteration = j,
    TP = TP, FP = FP, FN = FN, TN = TN,
    Power         = round(power,         4),
    Power_nominal = round(power_nominal, 4),
    Type_I_Error  = round(typeI,         4),
    Type_II_Error = round(typeII,        4),
    Precision     = round(precision,     4),
    F1_Score      = round(f1,            4),
    N_quantiles   = nq
  )
  all_performance_data[[j]] <- perf_table
  all_results_data[[j]]     <- df.qr
  print(perf_table)
  
  ## ---- Per-iteration writes (NOW INSIDE THE LOOP) ----
  fwrite(df.qr,
         paste0("QR_replicate_", j,
                "_pheno_", pheno_label,
                "_num_quantili_", nq,
                "_num_varianti_", Num_varianti, ".tsv"),
         sep = "\t")
  
  reference_results <- df.qr %>%
    filter(ID %in% reference_snps) %>%
    arrange(P_QR)
  
  fwrite(reference_results,
         paste0("QR_reference_SNP_results_pheno_", pheno_label,
                "_seed_", j, ".tsv"),
         sep = "\t")
}

stopCluster(cl)
registerDoSEQ()

# ============================================================
# 8. Aggregate
# ============================================================
all_performance_data <- do.call(rbind, all_performance_data)
all_results_combined <- do.call(rbind, all_results_data)

summary_stats <- all_performance_data %>%
  summarise(
    Mean_Power     = mean(Power,          na.rm = TRUE),
    SD_Power       = sd(Power,            na.rm = TRUE),
    Mean_Power_nom = mean(Power_nominal,  na.rm = TRUE),
    SD_Power_nom   = sd(Power_nominal,    na.rm = TRUE),
    Mean_F1        = mean(F1_Score,       na.rm = TRUE),
    SD_F1          = sd(F1_Score,         na.rm = TRUE),
    Mean_TypeI     = mean(Type_I_Error,   na.rm = TRUE),
    SD_TypeI       = sd(Type_I_Error,     na.rm = TRUE),
    Mean_TypeII    = mean(Type_II_Error,  na.rm = TRUE),
    Mean_Precision = mean(Precision,      na.rm = TRUE),
    N_reference    = n_ref
  )
print(summary_stats)

fwrite(all_performance_data,
       paste0("ALL_QR_iteration_performance_data_pheno_", pheno_label, ".tsv"),
       sep = "\t")
fwrite(all_results_combined,
       paste0("all_per_snp_QR_results_combined_pheno_", pheno_label, ".tsv"),
       sep = "\t")
write.table(summary_stats,
            paste0("QR_performance_summary_statistics_pheno_", pheno_label, ".txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

# ============================================================
# 9. Boxplots (with paste0 filenames)
# ============================================================
all_perf_pct <- all_performance_data %>%
  mutate(
    Power_pct         = Power * 100,
    Power_nominal_pct = Power_nominal * 100,
    F1_Score_pct      = F1_Score * 100,
    Type_I_Error_pct  = Type_I_Error * 100
  )

make_box <- function(df, col, title, fill) {
  ggplot(df, aes(x = "", y = .data[[col]])) +
    geom_boxplot(fill = fill, alpha = 0.7,
                 outlier.shape = 19, outlier.size = 1.5) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
    labs(title = title, x = NULL, y = "%") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

ggsave(paste0("boxplot_power_QR_pheno_", pheno_label, ".png"),
       make_box(all_perf_pct, "Power_pct",
                "Distribution of Power — QR", "#1f77b4"),
       width = 6, height = 6, dpi = 300)

ggsave(paste0("boxplot_power_nominal_QR_pheno_", pheno_label, ".png"),
       make_box(all_perf_pct, "Power_nominal_pct",
                "Distribution of Nominal Power — QR", "#9467bd"),
       width = 6, height = 6, dpi = 300)

ggsave(paste0("boxplot_f1score_QR_pheno_", pheno_label, ".png"),
       make_box(all_perf_pct, "F1_Score_pct",
                "Distribution of F1 Score — QR", "#2ca02c"),
       width = 6, height = 6, dpi = 300)

ggsave(paste0("boxplot_typeIerror_QR_pheno_", pheno_label, ".png"),
       make_box(all_perf_pct, "Type_I_Error_pct",
                "Distribution of Empirical FP proportion — QR", "#d62728"),
       width = 6, height = 6, dpi = 300)

plot_long <- all_perf_pct %>%
  select(Power_pct, F1_Score_pct, Type_I_Error_pct) %>%
  pivot_longer(everything(), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(
    Metric,
    levels = c("Power_pct", "F1_Score_pct", "Type_I_Error_pct"),
    labels = c("Power (FDR)", "F1 Score", "Type I Error")
  ))

p_combined <- ggplot(plot_long, aes(x = Metric, y = Value, fill = Metric)) +
  geom_boxplot(outlier.shape = 19, outlier.size = 1.2,
               alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.12, size = 0.6, alpha = 0.45) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(
    "Power (FDR)"  = "#1f77b4",
    "F1 Score"     = "#2ca02c",
    "Type I Error" = "#d62728"
  )) +
  labs(title = paste0("QR performance — N_PCS = 4, ",
                      ifelse(USE_INT_FOR_QR, "INT", "raw"), " phenotype"),
       x = NULL, y = "%") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    axis.text.x     = element_blank(),
    strip.text      = element_text(face = "bold", size = 12)
  )

ggsave(paste0("boxplot_power_F1_TypeI_combined_QR_pheno_", pheno_label, ".png"),
       p_combined, width = 11, height = 5, dpi = 300)

cat("\nAll outputs saved (phenotype =", pheno_label, ").\n")

# ============================================================
# 10. Diagnostic — merge QR vs LR (uses correct filename)
# ============================================================
qr  <- fread(paste0("all_per_snp_QR_results_combined_pheno_",
                    pheno_label, ".tsv"))
lr  <- fread("all_per_snp_results_linear_regression_INT_PCs4_combined.txt")
ref <- read.delim("reference_snps_used.tsv")$SNP

qr_ref <- qr[Iteration == 1 & ID %in% ref,
             .(ID, P_QR, FDR_significant)]
lr_ref <- lr[Iteration == 1 & ID %in% ref,
             .(ID, P_LR, FDR_significant)]

merged <- merge(qr_ref, lr_ref, by = "ID",
                suffixes = c("_QR", "_LR"))
print(merged)

##
## Cluster-based regression with 9 significant SNPs
##

library(doParallel)
library(foreach)
library(data.table)
library(dplyr)
library(QRank)
library(quantreg)
library(ggplot2)
options(scipen = 999)

setwd(getwd())

# ============================================================
# 1. Load data
# ============================================================
cat("Loading genotype data...\n")
input      <- fread("A_thaliana_TF10_NO_NA.txt")
sample_ids <- input[[1]]
geno       <- as.matrix(input[, 3:ncol(input)])
colnames(geno) <- colnames(input)[3:ncol(input)]
raw_pheno  <- input[[2]]                    # RAW (non-INT) phenotype
input <- NULL; gc()

cat("Samples:", nrow(geno), "  SNPs:", ncol(geno), "\n")

## Reference SNPs (9, fixed)
reference_snps <- read.delim("reference_snps_used.tsv")$SNP
n_ref <- length(reference_snps)
cat("Reference SNPs:", n_ref, "\n")

## Random SNP panels (100 x 241)
random_snps_wide <- fread("random_snps_per_iteration.tsv")

## Covariates (IID + PC1..PC4)
covariates <- fread("covariates_PCs4.tsv")
match_cov  <- match(sample_ids, covariates$IID)
if (any(is.na(match_cov))) stop("IID mismatch in covariates.")
covariates <- covariates[match_cov, ]

## Analysis dataframe
df_analysis <- data.frame(
  IID       = sample_ids,
  PHENOTYPE = raw_pheno,
  covariates[, paste0("PC", 1:4), with = FALSE]
)

cat("Data aligned. Phenotype is RAW (non-INT).\n")

# ============================================================
# 2. Parameters
# ============================================================
Num_varianti  <- 250
N_RANDOM      <- ncol(random_snps_wide) - 1        # 241
ult_iterazione <- nrow(random_snps_wide)           # 100
FDR_THRESHOLD <- 0.05
N_PCS         <- 4

## Cap the k-means grid at the number of distinct phenotype values
n_unique_pheno <- length(unique(df_analysis$PHENOTYPE))
cat("Distinct phenotype values:", n_unique_pheno, "\n")

cluster_sizes <- seq(3, min(100, n_unique_pheno), by = 1)
cat("k-means grid:", min(cluster_sizes), "to", max(cluster_sizes),
    "(", length(cluster_sizes), "values)\n")

# ============================================================
# 3. ACAT + QRank test (unchanged)
# ============================================================
acat <- function(pvals) {
  pvals = pvals[!is.na(pvals)]
  if (length(pvals) == 0) return(NA)
  pvals <- pmax(pmin(pvals, 1 - 1e-15), 1e-15)
  cauchy = tanpi(0.5 - pvals)
  stats  = mean(cauchy)
  pcauchy(q = stats, lower.tail = FALSE)
}

qrank.test <- function(phe, cov, snp, fit.residual, tau) {
  ltau = length(tau)
  x = as.matrix(snp); y = as.matrix(phe); n <- nrow(y)
  zz = cbind(rep(1, n), as.matrix(cov))
  VN = outer(tau, tau, function(t1, t2) pmin(t1, t2) - t1 * t2)
  xstar = lm(x ~ zz - 1)$residual
  SN = NULL
  for (i in 1:ltau) {
    idx <- ((i - 1) * n + 1):(i * n)
    ranks = 1 * (fit.residual$Residual[idx] > 0) - (1 - tau[i])
    SN = c(SN, as.matrix(t(xstar) %*% ranks))
  }
  VN2 = VN * as.numeric(t(xstar) %*% xstar)
  pvalue = pchisq(SN^2 / diag(VN2), 1, lower.tail = FALSE)
  names(pvalue) = tau
  list(tau = tau, quantile.pvalue = pvalue,
       cauchy.pvalue = acat(pvalue))
}

# ============================================================
# 4. Storage
# ============================================================
all_performance_data <- list()

# ============================================================
# 5. Main iteration loop
# ============================================================
for (j in seq_len(ult_iterazione)) {

  cat("\n============================================================\n")
  cat(">> Iteration", j, "/", ult_iterazione, "\n")
  cat("============================================================\n")

  ## Random SNPs for this iteration
  random_snps_j <- as.character(random_snps_wide[j, -1, with = FALSE])

  ## Build 250-SNP matrix
  geno.mat <- cbind(
    geno[, reference_snps, drop = FALSE],
    geno[, random_snps_j,  drop = FALSE]
  )
  colnames(geno.mat) <- c(reference_snps, random_snps_j)

  ## Parallel cluster (per iteration, as in original)
  ncores       <- detectCores()
  main_cluster <- makeCluster(ncores)
  registerDoParallel(main_cluster)

  clusterExport(main_cluster, c("qrank.test", "acat", "df_analysis",
                                "N_PCS"),
                envir = environment())

  clusterEvalQ(main_cluster, {
    library(quantreg); library(dplyr); library(data.table)
  })

  # ============================================================
  # Inner loop: k-means cluster sizes
  # ============================================================
  for (k in cluster_sizes) {

    ## Reset confusion results for this (j, k)
    confusion_results <- list(
      FT10 = list(TP = 0, FP = 0, FN = 0, TN = 0)
    )

    cat("\n  Cluster k =", k, "\n")

    y <- df_analysis$PHENOTYPE

    ## K-means on the phenotype
    kmeans_result <- kmeans(
      y,
      centers   = k,
      iter.max  = 500,
      nstart    = 10,
      algorithm = "MacQueen"
    )

    ## Empirical CDF of phenotype
    ecdf_func <- ecdf(y)

    ## Sort centroids and convert to tau levels
    centroids <- sapply(sort(unique(kmeans_result$cluster)), function(cl) {
      mean(y[kmeans_result$cluster == cl])
    })
    centroids_sorted <- sort(centroids)
    tau_levels       <- ecdf_func(centroids_sorted)
    tau_levels       <- pmax(pmin(tau_levels, 0.99), 0.01)

    ## Diagnostic (as in original)
    cat("    Unique clusters found:", length(unique(kmeans_result$cluster)),
        "| Unique tau levels:", length(unique(tau_levels)), "\n")

    ## Pre-allocate residual data frame
    df.res <- data.frame(
      Quantile = rep(tau_levels, each = nrow(df_analysis)),
      IID      = rep(df_analysis$IID, times = length(tau_levels)),
      Residual = numeric(length(tau_levels) * nrow(df_analysis))
    )

    ## Null model: PHENOTYPE ~ PC1 + PC2 + PC3 + PC4
    fml.null <- paste("PHENOTYPE ~",
                      paste(paste0("PC", 1:N_PCS), collapse = " + "))
    for (idx.qntl in seq_along(tau_levels)) {
      fit.null = rq(as.formula(fml.null),
                    tau = tau_levels[idx.qntl],
                    data = df_analysis,
                    method = "fn")
      idx <- ((idx.qntl - 1) * nrow(df_analysis) + 1):
             (idx.qntl * nrow(df_analysis))
      df.res$Residual[idx] <- resid(fit.null)
    }

    ## Export per-k objects
    clusterExport(main_cluster,
                  c("geno.mat", "tau_levels", "df.res"),
                  envir = environment())

    ## Parallel QRank over all 250 SNPs
    df.qr.parallel <- foreach(
      idx.geno = 1:ncol(geno.mat),
      .combine = function(...) data.table::rbindlist(list(...), fill = TRUE),
      .packages = c("quantreg", "dplyr", "data.table"),
      .errorhandling = "remove"
    ) %dopar% {

      df.qreg <- cbind(
        data.frame(PHENOTYPE = df_analysis$PHENOTYPE),
        df_analysis[, paste0("PC", 1:N_PCS)],
        data.frame(GENOTYPE = geno.mat[, idx.geno])
      )

      ## Per-quantile QR coefficients
      beta <- double()
      fml.qreg <- paste("PHENOTYPE ~",
                        paste(c(paste0("PC", 1:N_PCS), "GENOTYPE"),
                              collapse = " + "))
      for (idx.qntl in seq_along(tau_levels)) {
        qreg <- rq(as.formula(fml.qreg),
                   tau = tau_levels[idx.qntl],
                   data = df.qreg,
                   method = "fn")
        beta <- c(beta, qreg$coefficients["GENOTYPE"])
      }

      ## QRank test
      p <- qrank.test(
        phe          = df.qreg %>% select(PHENOTYPE),
        cov          = df.qreg %>% select(all_of(paste0("PC", 1:N_PCS))),
        snp          = df.qreg %>% select(GENOTYPE),
        fit.residual = df.res,
        tau          = tau_levels
      )

      data.frame(
        p$cauchy.pvalue,
        t(p$quantile.pvalue),
        t(beta),
        stringsAsFactors = FALSE
      )
    }

    if ("data.table" %in% class(df.qr.parallel)) {
      df.qr <- as.data.frame(df.qr.parallel)
    } else {
      df.qr <- df.qr.parallel
    }

    ## Column names: use sorted centroid values as labels (raw-phenotype units)
    qntl_labels <- round(centroids_sorted, 3)
    colnames(df.qr) <- c("P_QR",
                         paste0("P_QR",   qntl_labels),
                         paste0("BETA_QR", qntl_labels))
    df.qr <- df.qr %>% mutate(ID = colnames(geno.mat), .before = 1)

    ## FDR on the combined p-value
    df.qr$P_QR_corrected <- p.adjust(df.qr$P_QR, method = "fdr")
    df.qr$FDR_significant <- df.qr$P_QR_corrected < FDR_THRESHOLD
    df.qr$is_reference    <- df.qr$ID %in% reference_snps

    ## Confusion matrix
    significant  <- df.qr$FDR_significant
    is_reference <- df.qr$is_reference

    TP <- sum( significant &  is_reference, na.rm = TRUE)
    FP <- sum( significant & !is_reference, na.rm = TRUE)
    FN <- sum(!significant &  is_reference, na.rm = TRUE)
    TN <- sum(!significant & !is_reference, na.rm = TRUE)

    power     <- TP / (TP + FN + 1e-8)
    typeI     <- FP / (FP + TN + 1e-8)
    typeII    <- FN / (TP + FN + 1e-8)
    precision <- TP / (TP + FP + 1e-8)
    f1        <- 2 * precision * power / (precision + power + 1e-8)

    perf_table <- data.frame(
      Model         = "FT10",
      Iteration     = j,
      N_clusters    = k,
      TP = TP, FP = FP, FN = FN, TN = TN,
      Power         = round(power,     4),
      Type_I_Error  = round(typeI,     4),
      Type_II_Error = round(typeII,    4),
      Precision     = round(precision, 4),
      F1_Score      = round(f1,        4)
    )
    print(perf_table)

    all_performance_data[[length(all_performance_data) + 1]] <- perf_table

    ## Save per-(j,k) results
    fwrite(df.qr,
           paste0("QR_kmeans_seed_", j, "_k", k,
                  "_num_varianti_", Num_varianti,
                  "_pheno_raw.tsv"),
           sep = "\t")

    rm(df.res); gc()
  }  # end k-means cluster-size loop

  stopCluster(main_cluster)
  registerDoSEQ()
}  # end iteration loop

# ============================================================
# 6. Aggregate + save
# ============================================================
all_perf_combined <- do.call(rbind, all_performance_data)

fwrite(all_perf_combined,
       "ALL_QR_kmeans_performance_pheno_raw.tsv",
       sep = "\t")

summary_stats <- all_perf_combined %>%
  group_by(N_clusters) %>%
  summarise(
    Mean_Power    = mean(Power,        na.rm = TRUE),
    SD_Power      = sd(Power,          na.rm = TRUE),
    Mean_F1       = mean(F1_Score,     na.rm = TRUE),
    SD_F1         = sd(F1_Score,       na.rm = TRUE),
    Mean_TypeI    = mean(Type_I_Error, na.rm = TRUE),
    SD_TypeI      = sd(Type_I_Error,   na.rm = TRUE),
    .groups = "drop"
  )
print(summary_stats)
write.table(summary_stats,
            "QR_kmeans_summary_pheno_raw.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")

# ============================================================
# 7. Plots (performance vs k)
# ============================================================
all_perf_combined <- all_perf_combined %>%
  mutate(
    Power_pct        = Power        * 100,
    F1_Score_pct     = F1_Score     * 100,
    Type_I_Error_pct = Type_I_Error * 100,
    N_clusters_f     = factor(N_clusters)
  )

## Power vs k
p1 <- ggplot(all_perf_combined,
             aes(x = N_clusters_f, y = Power_pct)) +
  geom_boxplot(outlier.shape = 19, outlier.size = 1.2,
               fill = "#1f77b4", alpha = 0.7) +
  scale_x_discrete(drop = FALSE) +
  labs(title = "Power vs k (raw phenotype)",
       x = "Number of clusters (k)", y = "Power (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        plot.title  = element_text(hjust = 0.5, face = "bold"))
ggsave("boxplot_power_vs_k_pheno_raw.png",
       p1, width = 14, height = 5, dpi = 300)

## F1 vs k
p2 <- ggplot(all_perf_combined,
             aes(x = N_clusters_f, y = F1_Score_pct)) +
  geom_boxplot(outlier.shape = 19, outlier.size = 1.2,
               fill = "#2ca02c", alpha = 0.7) +
  scale_x_discrete(drop = FALSE) +
  labs(title = "F1 Score vs k (raw phenotype)",
       x = "Number of clusters (k)", y = "F1 Score (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        plot.title  = element_text(hjust = 0.5, face = "bold"))
ggsave("boxplot_f1score_vs_k_pheno_raw.png",
       p2, width = 14, height = 5, dpi = 300)

## Type I vs k
p3 <- ggplot(all_perf_combined,
             aes(x = N_clusters_f, y = Type_I_Error_pct)) +
  geom_boxplot(outlier.shape = 19, outlier.size = 1.2,
               fill = "#d62728", alpha = 0.7) +
  scale_x_discrete(drop = FALSE) +
  labs(title = "Empirical FP proportion vs k (raw phenotype)",
       x = "Number of clusters (k)", y = "Type I Error (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        plot.title  = element_text(hjust = 0.5, face = "bold"))
ggsave("boxplot_typeIerror_vs_k_pheno_raw.png",
       p3, width = 14, height = 5, dpi = 300)

cat("\nAll outputs saved.\n")
