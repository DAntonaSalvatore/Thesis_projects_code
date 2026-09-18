##
## Quantile regression with 10 significant SNPs
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
#setwd("C:\\Users\\Salvatore\\Desktop\\Dottorato\\studio_DL\\Quantile_Regression\\Nuovi_dati_tutorial_Laza\\prova_CR_funzionera")
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
