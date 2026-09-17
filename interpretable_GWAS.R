# INIZIO PARTE FENOTIPO CONTINUO

###
### PARTE RELATIVA ALLA FEATURE SELECTION E ALLA LINEAR REGRESSION (OLS)
###

library(data.table)
library(nestedcv)
library(glmnet)

# Load data
final_df <- fread("A_thaliana_TF10_NO_NA.txt")

# Train / test split
set.seed(123)
n <- nrow(final_df)
train_idx <- sample(seq_len(n), size = floor(0.8 * n))

train_df <- final_df[train_idx, ]
test_df  <- final_df[-train_idx, ]

# Separate X and y
y_train <- train_df[[2]]
X_train <- as.matrix(train_df[, 3:ncol(train_df)])

y_test <- test_df[[2]]
X_test <- as.matrix(test_df[, 3:ncol(test_df)])

# Nested CV LASSO. FEATURE SELECTION
fit <- nestcv.glmnet(
  x = X_train,
  y = y_train,
  family = "gaussian",
  alpha = 1,
  outer_k = 10,
  inner_k = 10,
  min_1se = 1
)

# Extract selected SNPs
coef_mat <- fit$final_coef
features_selected <- rownames(coef_mat)[coef_mat[, 1] != 0]
features_selected <- setdiff(features_selected, "(Intercept)")

# Reduced datasets
X_train_fs <- as.matrix(train_df[, ..features_selected])
X_test_fs  <- as.matrix(test_df[,  ..features_selected])

# turning in data.frame as required by lm()
df_train_fs <- as.data.frame(X_train_fs)
df_train_fs$y <- y_train

# LINEAR REGRESSION
lm_fit <- lm(y ~ ., data = df_train_fs)

summary(lm_fit)


# turn in data.frame
df_test_fs <- as.data.frame(X_test_fs)

# columns need to be in the same order of the training dataset
df_test_fs <- df_test_fs[, colnames(X_train_fs)]

# Predict sul test set
y_pred <- predict(lm_fit, newdata = df_test_fs)

# RMSE (Root Mean Squared Error)
rmse <- sqrt(mean((y_test - y_pred)^2))

# MAE (Mean Absolute Error)
mae <- mean(abs(y_test - y_pred))

# R² (proportion of explained variance in test set)
r2 <- 1 - sum((y_test - y_pred)^2) / sum((y_test - mean(y_test))^2)

cat("RMSE:", rmse, "\nMAE:", mae, "\nR²:", r2, "\n")

# for the interpretability. Mathematic formula of the model.
coefs <- coef(lm_fit)
coefs <- coefs[!is.na(coefs)]

formula_string <- paste0(
  "y = ",
  round(coefs[1], 3),
  paste(
    sprintf(" %+0.3f * %s", coefs[-1], names(coefs)[-1]),
    collapse = ""
  )
)

cat(formula_string)



split_df <- data.table(
  
  sample_id = final_df[[1]],
  
  set = ifelse(seq_len(n) %in% train_idx, "train", "test")
  
)



fwrite(split_df, "train_test_split.txt", sep = "\t")

fwrite(
  
  data.table(SNP = features_selected),
  
  "selected_snps.txt",
  
  sep = "\t"
  
)

coef_df <- data.table(
  
  SNP = names(coefs),
  
  beta = as.numeric(coefs)
  
)



fwrite(coef_df, "linear_model_coefficients.txt", sep = "\t")

coef_summary <- summary(lm_fit)$coefficients

coef_summary <- coef_summary[!is.na(coef_summary[,1]), ]



stats_df <- data.table(
  
  SNP = rownames(coef_summary),
  
  beta = coef_summary[,1],
  
  std_error = coef_summary[,2],
  
  t_value = coef_summary[,3],
  
  p_value = coef_summary[,4]
  
)



fwrite(stats_df, "linear_model_statistics.txt", sep = "\t")

perf_df <- data.table(
  
  RMSE = rmse,
  
  MAE = mae,
  
  R2 = r2
  
)



fwrite(perf_df, "test_performance.txt", sep = "\t")

pred_df <- data.table(
  
  sample_id = test_df[[1]],
  
  y_true = y_test,
  
  y_pred = y_pred
  
)



fwrite(pred_df, "test_predictions.txt", sep = "\t")

writeLines(formula_string, "linear_model_formula.txt")


###
### PARTE RELATIVA AL DECISION TREE
###


set.seed(123)
library(gramEvol)
library(data.table)
library(foreach)
library(doParallel)
library(parallel)
library(rpart)
library(rpart.plot)
library(caret)
library(Metrics)
library(ggplot2)

# ------------------------------------------------------------
# Set up parallel processing
# ------------------------------------------------------------
# Detect number of cores
num_cores <- detectCores()

# Register parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

setwd("C:\\Users\\Salvatore\\Desktop\\Dottorato\\Creta\\Progetto_Ricerca\\AtPolyDB_(call_method_75__Horton_et_al)\\Nuovo_tentativo_results\\linear_regression")
ID_train_test <- read.delim("train_test_split.txt")
input <- fread("A_thaliana_TF10_NO_NA.txt")
feature_selected <- read.delim("selected_snps.txt")
input <- as.data.frame(input)

# ------------------------------------------------------------
# Data Preparation
# ------------------------------------------------------------
# rebuild the original training and testing set
keep_cols <- c("IID", "FT10", feature_selected$SNP)
input_selected <- input[, colnames(input) %in% keep_cols, drop = FALSE]

training_list <- ID_train_test[ID_train_test$set == "train",]
testing_list <- ID_train_test[ID_train_test$set == "test",]
training_list <- training_list[,-2]
testing_list <- testing_list[,-2]

train_df <- input_selected[input_selected$IID %in% training_list, ]
test_df  <- input_selected[input_selected$IID %in% testing_list, ]

y_train <- train_df$FT10
X_train <- train_df[, !(names(train_df) %in% c("IID", "FT10"))]

y_test <- test_df$FT10
X_test <- test_df[, !(names(test_df) %in% c("IID", "FT10"))]

lm_data  <- cbind(phenotype = y_train, X_train)
test_data <- cbind(phenotype = y_test, X_test)

# Extract SNP names for grammar
snp_names <- colnames(X_train)

# Convert all SNPs to factors
for (snp in snp_names) {
  lm_data[[snp]]  <- as.factor(lm_data[[snp]])
  test_data[[snp]] <- as.factor(test_data[[snp]])
}

# ------------------------------------------------------------
# 1. Fix the error and analyze what's happening
# ------------------------------------------------------------
# Fix the error by handling constant predictions
calculate_safe_r2 <- function(actual, predicted) {
  if (sd(predicted, na.rm = TRUE) == 0) {
    return(0)  # If predictions are constant, R² is 0
  } else {
    return(cor(actual, predicted, use = "complete.obs")^2)
  }
}

# Function to calculate MAE
calculate_mae <- function(actual, predicted) {
  mean(abs(actual - predicted), na.rm = TRUE)
}

# ------------------------------------------------------------
# 3. Modified Decision Tree with Safe R2 Calculation
# ------------------------------------------------------------
cat("\n================== DECISION TREE ANALYSIS ==================\n")

# Create a safe prediction function
safe_predict <- function(model, data) {
  pred <- tryCatch({
    predict(model, data)
  }, error = function(e) {
    warning("Prediction failed, returning mean of training data")
    rep(mean(lm_data$phenotype, na.rm = TRUE), nrow(data))
  })
  return(pred)
}

# Revised parameter search with safe functions
best_cv_error <- Inf          # cross-validated error (lower is better)
best_params <- NULL
best_model <- NULL
results_grid <- data.frame()

# Try a range of parameters
for(cp_val in seq(0.00001, 0.1, length.out = 10)) {
  for(maxd in c(2, 3, 4, 5)) {
    for(mins in c(20, 50, 100)) {
      for(minb in c(round(mins/3), round(mins/2))) {
        
        tryCatch({
          dt_model <- rpart(
            phenotype ~ .,
            data = lm_data,
            method = "anova",
            control = rpart.control(
              minsplit = mins,
              minbucket = minb,
              cp = cp_val,
              maxdepth = maxd,
              xval = 5
            )
          )
          
          # >>> CHANGE 1: add this block right here <<<
          cv_error <- tryCatch({
            cpt <- dt_model$cptable
            if (is.null(cpt) || nrow(cpt) == 0) {
              Inf
            } else {
              min(cpt[, "xerror"], na.rm = TRUE)
            }
          }, error = function(e) Inf)
          
          # Predict with safety
          train_pred <- safe_predict(dt_model, lm_data)
          test_pred <- safe_predict(dt_model, test_data)
          
          # Calculate safe R2
          train_r2 <- calculate_safe_r2(lm_data$phenotype, train_pred)
          test_r2 <- calculate_safe_r2(test_data$phenotype, test_pred)
          
          # Calculate MAE
          train_mae <- calculate_mae(lm_data$phenotype, train_pred)
          test_mae <- calculate_mae(test_data$phenotype, test_pred)
          
          # Calculate RMSE
          train_rmse <- sqrt(mean((lm_data$phenotype - train_pred)^2, na.rm = TRUE))
          test_rmse <- sqrt(mean((test_data$phenotype - test_pred)^2, na.rm = TRUE))
          
          # Store results
          results_grid <- rbind(results_grid, data.frame(
            cp = cp_val,
            maxdepth = maxd,
            minsplit = mins,
            minbucket = minb,
            cv_error = cv_error,       # >>> CHANGE 2: new column <<<
            train_r2 = train_r2,
            test_r2 = test_r2,
            train_mae = train_mae,
            test_mae = test_mae,
            train_rmse = train_rmse,
            test_rmse = test_rmse,
            tree_size = ifelse(is.null(dt_model$frame), 1, nrow(dt_model$frame))
          ))
          
          # Update best model (still based on R², but could be changed to MAE)
          if(!is.na(cv_error) && cv_error < best_cv_error) {
            best_cv_error <- cv_error
            best_params <- list(cp = cp_val, maxdepth = maxd, 
                                minsplit = mins, minbucket = minb)
            best_model <- dt_model
            
            cat(sprintf("New best: cp=%.6f, depth=%d, minsplit=%d, minbucket=%d\n",
                        cp_val, maxd, mins, minb))
            cat(sprintf("  CV error=%.4f  (lower is better)\n", cv_error))
            cat(sprintf("  Train R2=%.4f, Test R2=%.4f\n",
                        train_r2, test_r2))
            cat(sprintf("  Train MAE=%.4f, Test MAE=%.4f\n",
                        train_mae, test_mae))
            cat(sprintf("  Train RMSE=%.4f, Test RMSE=%.4f, Size=%d\n\n",
                        train_rmse, test_rmse, ifelse(is.null(dt_model$frame), 1, nrow(dt_model$frame))))
          }
        }, error = function(e) {
          # Skip this parameter combination if it fails
          cat(sprintf("Failed for cp=%.6f, depth=%d: %s\n", 
                      cp_val, maxd, e$message))
        })
      }
    }
  }
}

# ------------------------------------------------------------
# 3.5 Distribution analysis of all grid-search runs
# ------------------------------------------------------------
cat("\n================== DISTRIBUTION ANALYSIS ==================\n")

if (nrow(results_grid) > 0) {
  
  # ---- 1. Build long-format data (base R, no extra deps) ----
  metrics_to_plot <- c("train_r2", "test_r2",
                       "train_rmse", "test_rmse",
                       "train_mae",  "test_mae")
  
  results_long <- data.frame(
    Metric = factor(rep(metrics_to_plot, each = nrow(results_grid)),
                    levels = metrics_to_plot,
                    labels = c("Train R²", "Test R²",
                               "Train RMSE", "Test RMSE",
                               "Train MAE",  "Test MAE")),
    Value  = unlist(results_grid[, metrics_to_plot], use.names = FALSE)
  )
  
  # ---- 2. Build the distribution plot ----
  p <- ggplot(results_long, aes(x = Metric, y = Value, fill = Metric)) +
    geom_violin(alpha = 0.55, trim = FALSE, colour = NA) +
    geom_boxplot(width = 0.18, alpha = 0.85, outlier.size = 0.8,
                 colour = "grey20") +
    stat_summary(fun = mean, geom = "point", shape = 23,
                 size = 2.6, fill = "white", colour = "black") +
    facet_wrap(~ Metric, scales = "free", ncol = 3) +
    theme_bw() +
    labs(
      title    = "Distribution of Metrics Across All Grid-Search Runs",
      subtitle = "White diamond = mean | Boxplot = median & quartiles | Violin = density",
      x = NULL, y = "Value"
    ) +
    theme(
      legend.position = "none",
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle   = element_text(hjust = 0.5),
      strip.text      = element_text(face = "bold")
    )
  
  ggsave("distribution_tree_graph.png", plot = p,
         width = 12, height = 8, dpi = 300)
  cat("Distribution graph saved to 'distribution_tree_graph.png'\n")
  
  # ---- 3. Summary statistics (mean / median / sd / min / max) ----
  summary_stats <- data.frame(
    Metric = metrics_to_plot,
    Mean   = sapply(results_grid[, metrics_to_plot], mean,   na.rm = TRUE),
    Median = sapply(results_grid[, metrics_to_plot], median, na.rm = TRUE),
    SD     = sapply(results_grid[, metrics_to_plot], sd,     na.rm = TRUE),
    Min    = sapply(results_grid[, metrics_to_plot], min,    na.rm = TRUE),
    Max    = sapply(results_grid[, metrics_to_plot], max,    na.rm = TRUE)
  )
  rownames(summary_stats) <- NULL
  
  cat("\nSummary statistics across all runs:\n")
  print(summary_stats)
  
  write.table(summary_stats, file = "distribution_summary_stats.txt",
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("\nSummary statistics saved to 'distribution_summary_stats.txt'\n")
  
  # ---- 4. A table in the same shape as 'decision_tree_metrics.txt'
  #         but with mean / median of the distribution instead of the best run ----
  distribution_metrics <- data.frame(
    Metric          = c("R²", "RMSE", "MAE"),
    Training_Mean   = c(mean(results_grid$train_r2,   na.rm = TRUE),
                        mean(results_grid$train_rmse, na.rm = TRUE),
                        mean(results_grid$train_mae,  na.rm = TRUE)),
    Training_Median = c(median(results_grid$train_r2,   na.rm = TRUE),
                        median(results_grid$train_rmse, na.rm = TRUE),
                        median(results_grid$train_mae,  na.rm = TRUE)),
    Testing_Mean    = c(mean(results_grid$test_r2,   na.rm = TRUE),
                        mean(results_grid$test_rmse, na.rm = TRUE),
                        mean(results_grid$test_mae,  na.rm = TRUE)),
    Testing_Median  = c(median(results_grid$test_r2,   na.rm = TRUE),
                        median(results_grid$test_rmse, na.rm = TRUE),
                        median(results_grid$test_mae,  na.rm = TRUE))
  )
  
  cat("\nMean / median metrics across the distribution:\n")
  print(distribution_metrics)
  
  write.table(distribution_metrics,
              file = "distribution_metrics.txt",
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Mean/median metrics saved to 'distribution_metrics.txt'\n")
  
} else {
  cat("results_grid is empty — nothing to plot.\n")
}

# ------------------------------------------------------------
# 4. Analyze and Save Results
# ------------------------------------------------------------
cat("\n================== FINAL RESULTS ==================\n")

if(!is.null(best_model)) {
  # Make final predictions
  final_train_pred <- safe_predict(best_model, lm_data)
  final_test_pred <- safe_predict(best_model, test_data)
  
  # Calculate metrics
  final_train_r2 <- calculate_safe_r2(lm_data$phenotype, final_train_pred)
  final_test_r2 <- calculate_safe_r2(test_data$phenotype, final_test_pred)
  
  final_train_rmse <- sqrt(mean((lm_data$phenotype - final_train_pred)^2, na.rm = TRUE))
  final_test_rmse <- sqrt(mean((test_data$phenotype - final_test_pred)^2, na.rm = TRUE))
  
  # Calculate MAE
  final_train_mae <- calculate_mae(lm_data$phenotype, final_train_pred)
  final_test_mae <- calculate_mae(test_data$phenotype, final_test_pred)
  
  cat(sprintf("Best Model Parameters:\n"))
  cat(sprintf("  cp: %.6f\n", best_params$cp))
  cat(sprintf("  maxdepth: %d\n", best_params$maxdepth))
  cat(sprintf("  minsplit: %d\n", best_params$minsplit))
  cat(sprintf("  minbucket: %d\n", best_params$minbucket))
  
  cat(sprintf("\nFinal Model Performance:\n"))
  cat(sprintf("  Training R²:    %.6f\n", final_train_r2))
  cat(sprintf("  Testing R²:     %.6f\n", final_test_r2))
  cat(sprintf("  Training RMSE:  %.6f\n", final_train_rmse))
  cat(sprintf("  Testing RMSE:   %.6f\n", final_test_rmse))
  cat(sprintf("  Training MAE:   %.6f\n", final_train_mae))
  cat(sprintf("  Testing MAE:    %.6f\n", final_test_mae))
  cat(sprintf("  Tree size:      %d nodes\n", 
              ifelse(is.null(best_model$frame), 1, nrow(best_model$frame))))
  
  # ------------------------------------------------------------
  # Save Results to Files
  # ------------------------------------------------------------
  
  # 1. Save performance metrics
  metrics_df <- data.frame(
    Metric = c("R²", "RMSE", "MAE"),
    Training = c(final_train_r2, final_train_rmse, final_train_mae),
    Testing = c(final_test_r2, final_test_rmse, final_test_mae)
  )
  
  write.table(metrics_df, file = "decision_tree_metrics.txt", 
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("\nMetrics saved to 'decision_tree_metrics.txt'\n")
  
  # 2. Save detailed results including parameters
  detailed_results <- data.frame(
    Parameter = c("cp", "maxdepth", "minsplit", "minbucket", 
                  "train_r2", "test_r2", "train_rmse", "test_rmse", 
                  "train_mae", "test_mae", "tree_size"),
    Value = c(best_params$cp, best_params$maxdepth, best_params$minsplit, best_params$minbucket,
              final_train_r2, final_test_r2, final_train_rmse, final_test_rmse,
              final_train_mae, final_test_mae, ifelse(is.null(best_model$frame), 1, nrow(best_model$frame)))
  )
  
  write.table(detailed_results, file = "decision_tree_detailed_results.txt", 
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Detailed results saved to 'decision_tree_detailed_results.txt'\n")
  
  # 3. Save the tree rules
  # Capture the tree summary
  tree_summary <- capture.output(print(best_model))
  tree_rules <- capture.output(rpart.rules(best_model))
  
  # Save tree structure
  sink("decision_tree_rules.txt")
  cat("================== DECISION TREE RULES ==================\n\n")
  
  cat("Tree Summary:\n")
  cat("-------------\n")
  for (line in tree_summary) {
    cat(line, "\n")
  }
  
  cat("\n\nTree Rules:\n")
  cat("-----------\n")
  for (line in tree_rules) {
    cat(line, "\n")
  }
  
  cat("\n\nVariable Importance:\n")
  cat("-------------------\n")
  if (!is.null(best_model$variable.importance)) {
    importance_df <- data.frame(
      Variable = names(best_model$variable.importance),
      Importance = as.numeric(best_model$variable.importance)
    )
    importance_df <- importance_df[order(-importance_df$Importance), ]
    for (i in 1:nrow(importance_df)) {
      cat(sprintf("%-20s: %.4f\n", importance_df$Variable[i], importance_df$Importance[i]))
    }
  } else {
    cat("No variable importance available\n")
  }
  
  # Also save the actual predictions for further analysis
  cat("\n\nPredictions (first 20 rows):\n")
  cat("------------------------------\n")
  pred_df <- data.frame(
    Actual = test_data$phenotype[1:min(20, length(test_data$phenotype))],
    Predicted = final_test_pred[1:min(20, length(final_test_pred))],
    Residual = test_data$phenotype[1:min(20, length(test_data$phenotype))] - 
      final_test_pred[1:min(20, length(final_test_pred))]
  )
  print(pred_df)
  
  sink()
  cat("Tree rules saved to 'decision_tree_rules.txt'\n")
  
  # 4. Save the full results grid for analysis
  write.table(results_grid, file = "decision_tree_grid_search_results.txt", 
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Grid search results saved to 'decision_tree_grid_search_results.txt'\n")
  
  # 5. Save the model object for future use
  saveRDS(best_model, file = "best_decision_tree_model.rds")
  cat("Model object saved to 'best_decision_tree_model.rds'\n")
  
  # 6. Save predictions for further analysis
  predictions_df <- data.frame(
    IID = c(train_df$IID, test_df$IID),
    Actual = c(lm_data$phenotype, test_data$phenotype),
    Predicted = c(final_train_pred, final_test_pred),
    Set = rep(c("train", "test"), c(nrow(train_df), nrow(test_df)))
  )
  
  write.table(predictions_df, file = "decision_tree_predictions.txt", 
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Predictions saved to 'decision_tree_predictions.txt'\n")
  
} else {
  cat("No valid model found. All parameter combinations failed.\n")
  
  # Try a simple baseline
  cat("\nTrying simple baseline model (mean prediction):\n")
  mean_pred <- mean(lm_data$phenotype, na.rm = TRUE)
  baseline_train_r2 <- calculate_safe_r2(lm_data$phenotype, 
                                         rep(mean_pred, nrow(lm_data)))
  baseline_test_r2 <- calculate_safe_r2(test_data$phenotype, 
                                        rep(mean_pred, nrow(test_data)))
  baseline_train_mae <- calculate_mae(lm_data$phenotype, 
                                      rep(mean_pred, nrow(lm_data)))
  baseline_test_mae <- calculate_mae(test_data$phenotype, 
                                     rep(mean_pred, nrow(test_data)))
  baseline_train_rmse <- sqrt(mean((lm_data$phenotype - rep(mean_pred, nrow(lm_data)))^2, na.rm = TRUE))
  baseline_test_rmse <- sqrt(mean((test_data$phenotype - rep(mean_pred, nrow(test_data)))^2, na.rm = TRUE))
  
  cat(sprintf("  Baseline (mean) Training R²:  %.6f\n", baseline_train_r2))
  cat(sprintf("  Baseline (mean) Testing R²:   %.6f\n", baseline_test_r2))
  cat(sprintf("  Baseline (mean) Training RMSE: %.6f\n", baseline_train_rmse))
  cat(sprintf("  Baseline (mean) Testing RMSE:  %.6f\n", baseline_test_rmse))
  cat(sprintf("  Baseline (mean) Training MAE: %.6f\n", baseline_train_mae))
  cat(sprintf("  Baseline (mean) Testing MAE:  %.6f\n", baseline_test_mae))
  
  # Save baseline results
  baseline_metrics <- data.frame(
    Metric = c("R²", "RMSE", "MAE"),
    Training = c(baseline_train_r2, baseline_train_rmse, baseline_train_mae),
    Testing = c(baseline_test_r2, baseline_test_rmse, baseline_test_mae)
  )
  
  write.table(baseline_metrics, file = "baseline_metrics.txt", 
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("\nBaseline metrics saved to 'baseline_metrics.txt'\n")
}

# ------------------------------------------------------------
# 5. Additional Diagnostic: Check if SNPs are informative
# ------------------------------------------------------------
cat("\n================== ADDITIONAL DIAGNOSTICS ==================\n")

# Check for constant predictions in the best model
if(!is.null(best_model)) {
  test_pred <- safe_predict(best_model, test_data)
  if(sd(test_pred, na.rm = TRUE) == 0) {
    cat("WARNING: Best model produces constant predictions on test set.\n")
    cat("This suggests the model is not learning meaningful patterns.\n")
    
    # Save warning to file
    sink("model_warnings.txt")
    cat("WARNING: Best model produces constant predictions on test set.\n")
    cat("This suggests the model is not learning meaningful patterns.\n")
    sink()
  }
}

# ------------------------------------------------------------
# 6. Summary of saved files
# ------------------------------------------------------------
cat("\n================== FILES SAVED ==================\n")
cat("1. decision_tree_metrics.txt - Performance metrics\n")
cat("2. decision_tree_detailed_results.txt - Detailed results with parameters\n")
cat("3. decision_tree_rules.txt - Tree structure and rules\n")
cat("4. decision_tree_grid_search_results.txt - All parameter combinations tried\n")
cat("5. best_decision_tree_model.rds - Saved model object (RDS format)\n")
cat("6. decision_tree_predictions.txt - Actual vs predicted values\n")
if (exists("model_warnings")) {
  cat("7. model_warnings.txt - Any warnings about the model\n")
}

# Stop parallel cluster
stopCluster(cl)


##
## PARTE RELATIVA A SVM LINEAR REGRESSION
##



set.seed(123)

# Libraries
library(data.table)
library(foreach)
library(doParallel)
library(parallel)
library(caret)
library(kernlab)
library(ggplot2)
library(gridExtra)
library(Metrics)


############################################################
# 2. PARALLEL COMPUTING
############################################################

num_cores <- detectCores() - 1

# Safety check
if (is.na(num_cores) || num_cores < 1) {
  num_cores <- 1
}

cl <- makeCluster(num_cores)
registerDoParallel(cl)


############################################################
# 3. WORKING DIRECTORY
############################################################

setwd(
  "C:/Users/Salvatore/Desktop/Dottorato/Creta/Progetto_Ricerca/AtPolyDB_(call_method_75__Horton_et_al)/Nuovo_tentativo_results/linear_regression"
)


############################################################
# 4. INPUT FILES
############################################################

# Train/test assignment
ID_train_test <- read.delim(
  "train_test_split.txt",
  stringsAsFactors = FALSE
)

# Phenotype + SNP dataset
input <- fread(
  "A_thaliana_TF10_NO_NA.txt",
  data.table = FALSE
)

# Selected SNP list
feature_selected <- read.delim(
  "selected_snps.txt",
  stringsAsFactors = FALSE
)


############################################################
# 5. STANDARDIZE SAMPLE-ID COLUMN NAME
############################################################

# Identify sample ID column in train/test file
id_col_tt <- intersect(
  c("IID", "id", "ID", "sample", "sample_id"),
  colnames(ID_train_test)
)

if (length(id_col_tt) != 1) {
  stop(
    "Could not uniquely identify the sample ID column in train_test_split.txt"
  )
}

colnames(ID_train_test)[
  colnames(ID_train_test) == id_col_tt
] <- "IID"


# Identify sample ID column in input file
id_col_input <- intersect(
  c("IID", "id", "ID", "sample", "sample_id"),
  colnames(input)
)

if (length(id_col_input) != 1) {
  stop(
    "Could not uniquely identify the sample ID column in A_thaliana_TF10_NO_NA.txt"
  )
}

colnames(input)[
  colnames(input) == id_col_input
] <- "IID"


# Make IDs character
ID_train_test$IID <- as.character(ID_train_test$IID)
input$IID <- as.character(input$IID)


############################################################
# 6. CHECK SELECTED SNP COLUMN
############################################################

if (!"SNP" %in% colnames(feature_selected)) {
  stop(
    "The file selected_snps.txt must contain a column called 'SNP'."
  )
}

selected_snps <- as.character(
  feature_selected$SNP
)

selected_snps <- selected_snps[
  !is.na(selected_snps) &
    selected_snps != ""
]


cat("\nNumber of selected SNPs:", length(selected_snps), "\n")


############################################################
# 7. CHECK THAT SNPs EXIST IN INPUT DATA
############################################################

missing_snps <- setdiff(
  selected_snps,
  colnames(input)
)

if (length(missing_snps) > 0) {
  
  cat(
    "\nERROR: The following selected SNPs are missing from the input file:\n"
  )
  
  print(missing_snps)
  
  stop(
    "Some selected SNPs are not present in A_thaliana_TF10_NO_NA.txt"
  )
}


############################################################
# 8. SELECT PHENOTYPE + SNPs
############################################################

keep_cols <- c(
  "IID",
  "FT10",
  selected_snps
)

input_selected <- input[
  ,
  keep_cols,
  drop = FALSE
]


############################################################
# 9. IDENTIFY TRAINING AND TESTING SAMPLES
############################################################

if (!"set" %in% colnames(ID_train_test)) {
  stop(
    "train_test_split.txt must contain a column called 'set'."
  )
}


training_ids <- ID_train_test$IID[
  ID_train_test$set == "train"
]

testing_ids <- ID_train_test$IID[
  ID_train_test$set == "test"
]


############################################################
# 10. CREATE TRAINING AND TESTING DATASETS
############################################################

train_df <- input_selected[
  input_selected$IID %in% training_ids,
  ,
  drop = FALSE
]

test_df <- input_selected[
  input_selected$IID %in% testing_ids,
  ,
  drop = FALSE
]


############################################################
# 11. CHECK SAMPLE NUMBERS
############################################################

cat("\nTraining samples:", nrow(train_df), "\n")
cat("Testing samples :", nrow(test_df), "\n")

if (nrow(train_df) == 0) {
  stop("No training samples were found.")
}

if (nrow(test_df) == 0) {
  stop("No testing samples were found.")
}


############################################################
# 12. EXTRACT PHENOTYPE
############################################################

y_train <- train_df$FT10
y_test <- test_df$FT10


############################################################
# 13. EXTRACT SNP MATRICES
############################################################

X_train <- train_df[
  ,
  !(colnames(train_df) %in% c("IID", "FT10")),
  drop = FALSE
]

X_test <- test_df[
  ,
  !(colnames(test_df) %in% c("IID", "FT10")),
  drop = FALSE
]


############################################################
# 14. CONVERT SNPs TO NUMERIC
############################################################

convert_snp <- function(x) {
  
  if (is.factor(x)) {
    return(as.numeric(as.character(x)))
  }
  
  return(as.numeric(x))
}


X_train[] <- lapply(
  X_train,
  convert_snp
)

X_test[] <- lapply(
  X_test,
  convert_snp
)


############################################################
# 15. CHECK FOR NA VALUES
############################################################

if (anyNA(X_train)) {
  
  cat("\nNA values detected in X_train.\n")
  
  print(
    colSums(is.na(X_train))
  )
  
  stop(
    "Training SNP matrix contains NA values."
  )
}


if (anyNA(X_test)) {
  
  cat("\nNA values detected in X_test.\n")
  
  print(
    colSums(is.na(X_test))
  )
  
  stop(
    "Testing SNP matrix contains NA values."
  )
}


############################################################
# 16. CHECK PHENOTYPE
############################################################

if (anyNA(y_train)) {
  stop("NA values detected in training FT10.")
}

if (anyNA(y_test)) {
  stop("NA values detected in testing FT10.")
}


############################################################
# 17. CREATE CARET DATA FRAMES
############################################################

train_data <- data.frame(
  FT10 = y_train,
  X_train,
  check.names = FALSE
)

test_data <- data.frame(
  FT10 = y_test,
  X_test,
  check.names = FALSE
)


############################################################
# 18. CROSS-VALIDATION SETTINGS
############################################################

ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 5,
  allowParallel = TRUE
)


############################################################
# 19. SVM HYPERPARAMETER GRID
############################################################

svm_grid <- expand.grid(
  cost = c(
    0.001,
    0.01,
    0.1,
    1
  )
)


############################################################
# 20. TRAIN LINEAR SVM
############################################################

cat("\n========================================\n")
cat("TRAINING LINEAR SVM\n")
cat("========================================\n\n")


svm_model <- train(
  FT10 ~ .,
  data = train_data,
  method = "svmLinear2",
  trControl = ctrl,
  tuneGrid = svm_grid,
  
  # IMPORTANT:
  # caret calculates these using the TRAINING DATA ONLY
  preProcess = c(
    "center",
    "scale"
  )
)


############################################################
# 21. DISPLAY BEST MODEL
############################################################

cat("\n========================================\n")
cat("BEST MODEL\n")
cat("========================================\n\n")

print(
  svm_model$bestTune
)


############################################################
# 22. EXTRACT FINAL KERNLAB MODEL
############################################################

final_model <- svm_model$finalModel


cat("\n========================================\n")
cat("FINAL MODEL CLASS\n")
cat("========================================\n\n")

print(
  class(final_model)
)


############################################################
# 23. DISPLAY MODEL STRUCTURE
############################################################

cat("\n========================================\n")
cat("FINAL MODEL STRUCTURE\n")
cat("========================================\n\n")

print(
  slotNames(final_model)
)


############################################################
# 24. IMPORTANT:
#     DO NOT USE:
#
#     final_model$coef
#
#     It returns NULL because coef is not stored
#     as a normal list/data-frame element.
############################################################


############################################################
# 25. EXTRACT PREPROCESSING PARAMETERS
############################################################

preprocess_mean <- svm_model$preProcess$mean
preprocess_sd <- svm_model$preProcess$std


############################################################
# 26. CHECK PREPROCESSING PARAMETERS
############################################################

cat("\n========================================\n")
cat("CENTERING VALUES (MEANS)\n")
cat("========================================\n\n")

print(
  preprocess_mean
)


cat("\n========================================\n")
cat("SCALING VALUES (SDs)\n")
cat("========================================\n\n")

print(
  preprocess_sd
)


############################################################
# 27. CHECK THAT THE ORDER IS IDENTICAL
############################################################

if (!identical(
  names(preprocess_mean),
  names(preprocess_sd)
)) {
  
  stop(
    "The names/order of the preprocessing means and SDs do not match."
  )
}


############################################################
# 28. CHECK THAT ALL SNPs HAVE PREPROCESSING PARAMETERS
############################################################

if (!all(
  colnames(X_train) %in%
  names(preprocess_mean)
)) {
  
  stop(
    "Some SNPs do not have preprocessing means."
  )
}

if (!all(
  colnames(X_train) %in%
  names(preprocess_sd)
)) {
  
  stop(
    "Some SNPs do not have preprocessing SDs."
  )
}


############################################################
# 29. REORDER PARAMETERS TO MATCH X MATRIX
############################################################

preprocess_mean <- preprocess_mean[
  colnames(X_train)
]

preprocess_sd <- preprocess_sd[
  colnames(X_train)
]


############################################################
# 30. CHECK FOR ZERO SD
############################################################

if (any(preprocess_sd == 0)) {
  
  cat(
    "\nVariables with zero standard deviation:\n"
  )
  
  print(
    names(preprocess_sd)[preprocess_sd == 0]
  )
  
  stop(
    "At least one SNP has zero variance and cannot be scaled."
  )
}


############################################################
# 31. PREDICTIONS FROM CARET SVM
############################################################

pred_train <- predict(
  svm_model,
  newdata = train_data
)

pred_test <- predict(
  svm_model,
  newdata = test_data
)


############################################################
# 32. PERFORMANCE METRICS
############################################################

train_R2 <- R2(
  pred_train,
  y_train
)

test_R2 <- R2(
  pred_test,
  y_test
)

train_RMSE <- rmse(
  y_train,
  pred_train
)

test_RMSE <- rmse(
  y_test,
  pred_test
)

train_MAE <- mae(
  y_train,
  pred_train
)

test_MAE <- mae(
  y_test,
  pred_test
)


############################################################
# 33. PRINT PERFORMANCE
############################################################

cat("\n========================================\n")
cat("MODEL PERFORMANCE\n")
cat("========================================\n\n")

cat(
  "Train R2   :", train_R2, "\n"
)

cat(
  "Test R2    :", test_R2, "\n"
)

cat(
  "Train RMSE :", train_RMSE, "\n"
)

cat(
  "Test RMSE  :", test_RMSE, "\n"
)

cat(
  "Train MAE  :", train_MAE, "\n"
)

cat(
  "Test MAE   :", test_MAE, "\n"
)


############################################################
# 34. EXTRACT EFFECTIVE LINEAR COEFFICIENTS


reference_X <- as.data.frame(
  as.list(preprocess_mean),
  check.names = FALSE
)

reference_X <- reference_X[
  ,
  colnames(X_train),
  drop = FALSE
]


############################################################
# 36. REFERENCE PREDICTION


reference_prediction <- as.numeric(
  predict(
    svm_model,
    newdata = reference_X
  )
)


############################################################
# 37. INITIALIZE COEFFICIENT VECTOR
############################################################

beta_standardized <- numeric(
  length = ncol(X_train)
)

names(beta_standardized) <- colnames(X_train)


############################################################
# 38. CALCULATE EACH STANDARDIZED SNP COEFFICIENT
############################################################

for (j in seq_along(beta_standardized)) {
  
  snp_name <- colnames(X_train)[j]
  
  # Start with all SNPs at training means
  test_point <- reference_X
  
  # Move SNP j upward by exactly one SD
  test_point[[snp_name]] <-
    preprocess_mean[j] +
    preprocess_sd[j]
  
  # Prediction at Z_j = 1
  prediction_plus_one_sd <- as.numeric(
    predict(
      svm_model,
      newdata = test_point
    )
  )
  
  # Difference = beta_j
  beta_standardized[j] <-
    prediction_plus_one_sd -
    reference_prediction
}


############################################################
# 39. DISPLAY STANDARDIZED COEFFICIENTS
############################################################

cat("\n========================================\n")
cat("STANDARDIZED SVM COEFFICIENTS\n")
cat("========================================\n\n")

print(
  beta_standardized
)


############################################################
# 40. ORIGINAL-SCALE COEFFICIENTS


beta_original <- beta_standardized /
  preprocess_sd

names(beta_original) <- colnames(X_train)


############################################################
# 41. ORIGINAL-SCALE INTERCEPT


intercept_original <-
  reference_prediction -
  sum(
    beta_original *
      preprocess_mean
  )


############################################################
# 42. CREATE COEFFICIENT TABLE
############################################################

coefficient_table <- data.frame(
  
  SNP = colnames(X_train),
  
  Mean = as.numeric(
    preprocess_mean
  ),
  
  SD = as.numeric(
    preprocess_sd
  ),
  
  Coefficient_standardized =
    as.numeric(
      beta_standardized
    ),
  
  Coefficient_original =
    as.numeric(
      beta_original
    ),
  
  stringsAsFactors = FALSE,
  
  check.names = FALSE
)


############################################################
# 43. ADD ABSOLUTE COEFFICIENT
############################################################

coefficient_table$Absolute_coefficient_standardized <-
  abs(
    coefficient_table$Coefficient_standardized
  )


coefficient_table$Absolute_coefficient_original <-
  abs(
    coefficient_table$Coefficient_original
  )


############################################################
# 44. SORT BY ABSOLUTE EFFECT
############################################################

coefficient_table_sorted <-
  coefficient_table[
    order(
      -coefficient_table$Absolute_coefficient_standardized
    ),
    ,
    drop = FALSE
  ]


############################################################
# 45. DISPLAY SNP COEFFICIENTS
############################################################

cat("\n========================================\n")
cat("SNP COEFFICIENTS\n")
cat("========================================\n\n")

print(
  coefficient_table_sorted,
  row.names = FALSE
)


############################################################
# 46. DISPLAY INTERCEPTS
############################################################

cat("\n========================================\n")
cat("INTERCEPTS\n")
cat("========================================\n\n")

cat(
  "Standardized-space intercept:",
  reference_prediction,
  "\n"
)

cat(
  "Original-scale intercept:",
  intercept_original,
  "\n"
)


############################################################
# 47. BUILD STANDARDIZED MATHEMATICAL FORMULA
############################################################

standardized_formula <- paste0(
  "FT10_hat = ",
  format(
    reference_prediction,
    digits = 12
  )
)


for (j in seq_along(beta_standardized)) {
  
  coefficient <- beta_standardized[j]
  
  sign_symbol <-
    ifelse(
      coefficient >= 0,
      " + ",
      " - "
    )
  
  standardized_formula <- paste0(
    standardized_formula,
    sign_symbol,
    format(
      abs(coefficient),
      digits = 12
    ),
    " * (",
    colnames(X_train)[j],
    " - ",
    format(
      preprocess_mean[j],
      digits = 12
    ),
    ") / ",
    format(
      preprocess_sd[j],
      digits = 12
    )
  )
}


############################################################
# 48. BUILD ORIGINAL-SCALE MATHEMATICAL FORMULA
############################################################

original_formula <- paste0(
  "FT10_hat = ",
  format(
    intercept_original,
    digits = 12
  )
)


for (j in seq_along(beta_original)) {
  
  coefficient <- beta_original[j]
  
  sign_symbol <-
    ifelse(
      coefficient >= 0,
      " + ",
      " - "
    )
  
  original_formula <- paste0(
    original_formula,
    sign_symbol,
    format(
      abs(coefficient),
      digits = 12
    ),
    " * ",
    colnames(X_train)[j]
  )
}


############################################################
# 49. PRINT THE FORMULAS
############################################################

cat("\n========================================\n")
cat("MATHEMATICAL FORMULA\n")
cat("========================================\n\n")


cat(
  "STANDARDIZED FORM:\n\n"
)

cat(
  standardized_formula,
  "\n\n"
)


cat(
  "ORIGINAL SNP SCALE:\n\n"
)

cat(
  original_formula,
  "\n\n"
)


############################################################
# 50. RECONSTRUCT TRAINING DATA IN STANDARDIZED SPACE
############################################################

# IMPORTANT:
# Convert explicitly to a numeric matrix.
# This is required for matrix multiplication (%*%).

X_train_matrix <- as.matrix(
  data.frame(
    lapply(
      X_train,
      function(x) as.numeric(x)
    ),
    check.names = FALSE
  )
)

storage.mode(X_train_matrix) <- "numeric"


############################################################
# 51. STANDARDIZE TRAINING MATRIX
############################################################

X_train_standardized <- sweep(
  X_train_matrix,
  2,
  as.numeric(preprocess_mean),
  FUN = "-"
)

X_train_standardized <- sweep(
  X_train_standardized,
  2,
  as.numeric(preprocess_sd),
  FUN = "/"
)


############################################################
# 52. FORCE COEFFICIENTS TO NUMERIC VECTOR
############################################################

beta_standardized <- as.numeric(
  beta_standardized
)

names(beta_standardized) <- colnames(
  X_train
)


############################################################
# 53. CHECK DIMENSIONS AND DATA TYPES
############################################################

cat("\n========================================\n")
cat("FORMULA MATRIX CHECK\n")
cat("========================================\n\n")

cat(
  "X_train_standardized class:",
  class(X_train_standardized),
  "\n"
)

cat(
  "X_train_standardized mode:",
  mode(X_train_standardized),
  "\n"
)

cat(
  "Number of rows:",
  nrow(X_train_standardized),
  "\n"
)

cat(
  "Number of columns:",
  ncol(X_train_standardized),
  "\n"
)

cat(
  "Number of coefficients:",
  length(beta_standardized),
  "\n"
)

cat(
  "Coefficient mode:",
  mode(beta_standardized),
  "\n\n"
)


if (!is.numeric(X_train_standardized)) {
  stop(
    "X_train_standardized is still not numeric."
  )
}

if (!is.numeric(beta_standardized)) {
  stop(
    "beta_standardized is still not numeric."
  )
}

if (
  ncol(X_train_standardized) !=
  length(beta_standardized)
) {
  
  stop(
    paste0(
      "Dimension mismatch: ",
      ncol(X_train_standardized),
      " SNP columns but ",
      length(beta_standardized),
      " coefficients."
    )
  )
}


############################################################
# 54. FORMULA-BASED TRAIN PREDICTION
############################################################

pred_train_formula <-
  as.numeric(
    reference_prediction +
      X_train_standardized %*%
      beta_standardized
  )


############################################################
# 55. CHECK TRAIN PREDICTIONS
############################################################

cat("\nFirst five SVM predictions:\n")

print(
  head(pred_train, 5)
)

cat("\nFirst five formula predictions:\n")

print(
  head(pred_train_formula, 5)
)

cat("\nDifference between SVM and formula:\n")

print(
  head(
    pred_train - pred_train_formula,
    5
  )
)


############################################################
# 56. RECONSTRUCT TEST DATA IN STANDARDIZED SPACE
############################################################

X_test_matrix <- as.matrix(
  data.frame(
    lapply(
      X_test,
      function(x) as.numeric(x)
    ),
    check.names = FALSE
  )
)

storage.mode(X_test_matrix) <- "numeric"


############################################################
# 57. STANDARDIZE TEST MATRIX
############################################################

X_test_standardized <- sweep(
  X_test_matrix,
  2,
  as.numeric(preprocess_mean),
  FUN = "-"
)

X_test_standardized <- sweep(
  X_test_standardized,
  2,
  as.numeric(preprocess_sd),
  FUN = "/"
)


############################################################
# 58. FORMULA-BASED TEST PREDICTION
############################################################

pred_test_formula <-
  as.numeric(
    reference_prediction +
      X_test_standardized %*%
      beta_standardized
  )


############################################################
# 59. CHECK TEST PREDICTIONS
############################################################

cat("\nFirst five SVM test predictions:\n")

print(
  head(pred_test, 5)
)

cat("\nFirst five formula test predictions:\n")

print(
  head(pred_test_formula, 5)
)

cat("\nDifference between SVM and formula:\n")

print(
  head(
    pred_test - pred_test_formula,
    5
  )
)


############################################################
# 54. COMPARE SVM VS RECONSTRUCTED FORMULA
############################################################

train_difference <-
  pred_train -
  pred_train_formula

test_difference <-
  pred_test -
  pred_test_formula


############################################################
# 55. FORMULA RECONSTRUCTION ACCURACY
############################################################

max_train_difference <-
  max(
    abs(train_difference)
  )

max_test_difference <-
  max(
    abs(test_difference)
  )

mean_train_difference <-
  mean(
    abs(train_difference)
  )

mean_test_difference <-
  mean(
    abs(test_difference)
  )


############################################################
# 56. PRINT RECONSTRUCTION CHECK
############################################################

cat("\n========================================\n")
cat("FORMULA RECONSTRUCTION CHECK\n")
cat("========================================\n\n")


cat(
  "Maximum absolute difference - TRAIN:",
  max_train_difference,
  "\n"
)

cat(
  "Mean absolute difference - TRAIN:",
  mean_train_difference,
  "\n\n"
)

cat(
  "Maximum absolute difference - TEST:",
  max_test_difference,
  "\n"
)

cat(
  "Mean absolute difference - TEST:",
  mean_test_difference,
  "\n\n"
)


############################################################
# 57. CHECK IF RECONSTRUCTION IS NUMERICALLY IDENTICAL
############################################################

tolerance <- 1e-8

formula_matches_train <-
  max_train_difference < tolerance

formula_matches_test <-
  max_test_difference < tolerance


cat(
  "Formula reproduces TRAIN predictions:",
  formula_matches_train,
  "\n"
)

cat(
  "Formula reproduces TEST predictions:",
  formula_matches_test,
  "\n"
)


############################################################
# 60. CALCULATE TRAINING SNP CONTRIBUTIONS
############################################################

train_contributions <- sweep(
  X_train_standardized,
  2,
  beta_standardized,
  FUN = "*"
)

train_contributions <- as.data.frame(
  train_contributions,
  check.names = FALSE
)

colnames(train_contributions) <- colnames(
  X_train
)


############################################################
# 61. CALCULATE TEST SNP CONTRIBUTIONS
############################################################

test_contributions <- sweep(
  X_test_standardized,
  2,
  beta_standardized,
  FUN = "*"
)

test_contributions <- as.data.frame(
  test_contributions,
  check.names = FALSE
)

colnames(test_contributions) <- colnames(
  X_test
)


############################################################
# 62. CHECK CONTRIBUTION SUMS
############################################################

train_contribution_sum <-
  reference_prediction +
  rowSums(
    train_contributions
  )

test_contribution_sum <-
  reference_prediction +
  rowSums(
    test_contributions
  )


cat("\n========================================\n")
cat("CONTRIBUTION CHECK\n")
cat("========================================\n\n")

cat(
  "Maximum train difference:",
  max(
    abs(
      pred_train -
        train_contribution_sum
    )
  ),
  "\n"
)

cat(
  "Maximum test difference:",
  max(
    abs(
      pred_test -
        test_contribution_sum
    )
  ),
  "\n"
)


############################################################
# 63. ADD SAMPLE IDs AND PREDICTIONS
############################################################

train_contribution_table <- data.frame(
  
  IID = train_df$IID,
  
  FT10_observed = y_train,
  
  Intercept = reference_prediction,
  
  train_contributions,
  
  Prediction_SVM = pred_train,
  
  Prediction_formula = pred_train_formula,
  
  stringsAsFactors = FALSE,
  
  check.names = FALSE
)


test_contribution_table <- data.frame(
  
  IID = test_df$IID,
  
  FT10_observed = y_test,
  
  Intercept = reference_prediction,
  
  test_contributions,
  
  Prediction_SVM = pred_test,
  
  Prediction_formula = pred_test_formula,
  
  stringsAsFactors = FALSE,
  
  check.names = FALSE
)

############################################################
# 61. ADD PREDICTION ERROR
############################################################

train_contribution_table$Error <-
  train_contribution_table$FT10_observed -
  train_contribution_table$Prediction_SVM


train_contribution_table$Absolute_error <-
  abs(
    train_contribution_table$Error
  )


test_contribution_table$Error <-
  test_contribution_table$FT10_observed -
  test_contribution_table$Prediction_SVM


test_contribution_table$Absolute_error <-
  abs(
    test_contribution_table$Error
  )


############################################################
# 62. CREATE STANDARDIZED GENOTYPE TABLES
############################################################

train_standardized_table <-
  data.frame(
    IID = train_df$IID,
    X_train_standardized,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


test_standardized_table <-
  data.frame(
    IID = test_df$IID,
    X_test_standardized,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


############################################################
# 63. CREATE ORIGINAL GENOTYPE TABLES
############################################################

train_genotype_table <-
  data.frame(
    IID = train_df$IID,
    X_train,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


test_genotype_table <-
  data.frame(
    IID = test_df$IID,
    X_test,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


############################################################
# 64. CREATE PREDICTION TABLE
############################################################

train_prediction_table <-
  data.frame(
    IID = train_df$IID,
    Observed_FT10 = y_train,
    Predicted_FT10 = pred_train,
    Formula_FT10 = pred_train_formula,
    Error = y_train - pred_train,
    Absolute_error = abs(y_train - pred_train),
    stringsAsFactors = FALSE
  )


test_prediction_table <-
  data.frame(
    IID = test_df$IID,
    Observed_FT10 = y_test,
    Predicted_FT10 = pred_test,
    Formula_FT10 = pred_test_formula,
    Error = y_test - pred_test,
    Absolute_error = abs(y_test - pred_test),
    stringsAsFactors = FALSE
  )


############################################################
# 65. CREATE PERFORMANCE TABLE
############################################################

performance_table <- data.frame(
  
  Dataset = c(
    "Training",
    "Testing"
  ),
  
  R2 = c(
    train_R2,
    test_R2
  ),
  
  RMSE = c(
    train_RMSE,
    test_RMSE
  ),
  
  MAE = c(
    train_MAE,
    test_MAE
  ),
  
  N = c(
    length(y_train),
    length(y_test)
  ),
  
  stringsAsFactors = FALSE
)


############################################################
# 66. CREATE MODEL INFORMATION TABLE
############################################################

model_information <- data.frame(
  
  Parameter = c(
    "Model",
    "Kernel",
    "Cost",
    "Number_of_SNPs",
    "Training_samples",
    "Testing_samples",
    "Standardized_intercept",
    "Original_scale_intercept",
    "Formula_matches_training",
    "Formula_matches_testing",
    "Max_training_difference",
    "Max_testing_difference"
  ),
  
  Value = c(
    
    "svmLinear2",
    "linear",
    
    as.character(
      svm_model$bestTune$cost
    ),
    
    as.character(
      ncol(X_train)
    ),
    
    as.character(
      nrow(X_train)
    ),
    
    as.character(
      nrow(X_test)
    ),
    
    as.character(
      reference_prediction
    ),
    
    as.character(
      intercept_original
    ),
    
    as.character(
      formula_matches_train
    ),
    
    as.character(
      formula_matches_test
    ),
    
    as.character(
      max_train_difference
    ),
    
    as.character(
      max_test_difference
    )
  ),
  
  stringsAsFactors = FALSE
)


############################################################
# 67. SAVE COEFFICIENT TABLE
############################################################

write.csv(
  coefficient_table,
  "SVM_SNP_coefficients.csv",
  row.names = FALSE
)


write.csv(
  coefficient_table_sorted,
  "SVM_SNP_coefficients_sorted_by_effect.csv",
  row.names = FALSE
)


############################################################
# 68. SAVE PREPROCESSING PARAMETERS
############################################################

preprocessing_table <- data.frame(
  
  SNP = colnames(X_train),
  
  Mean = as.numeric(
    preprocess_mean
  ),
  
  SD = as.numeric(
    preprocess_sd
  ),
  
  stringsAsFactors = FALSE,
  
  check.names = FALSE
)


write.csv(
  preprocessing_table,
  "SVM_preprocessing_parameters.csv",
  row.names = FALSE
)


############################################################
# 69. SAVE STANDARDIZED GENOTYPES
############################################################

write.csv(
  train_standardized_table,
  "SVM_train_standardized_SNPs.csv",
  row.names = FALSE
)


write.csv(
  test_standardized_table,
  "SVM_test_standardized_SNPs.csv",
  row.names = FALSE
)


############################################################
# 70. SAVE ORIGINAL GENOTYPES
############################################################

write.csv(
  train_genotype_table,
  "SVM_train_original_SNPs.csv",
  row.names = FALSE
)


write.csv(
  test_genotype_table,
  "SVM_test_original_SNPs.csv",
  row.names = FALSE
)


############################################################
# 71. SAVE INDIVIDUAL SNP CONTRIBUTIONS
############################################################

write.csv(
  train_contribution_table,
  "SVM_train_SNP_contributions.csv",
  row.names = FALSE
)


write.csv(
  test_contribution_table,
  "SVM_test_SNP_contributions.csv",
  row.names = FALSE
)


############################################################
# 72. SAVE PREDICTIONS
############################################################

write.csv(
  train_prediction_table,
  "SVM_train_predictions.csv",
  row.names = FALSE
)


write.csv(
  test_prediction_table,
  "SVM_test_predictions.csv",
  row.names = FALSE
)


############################################################
# 73. SAVE PERFORMANCE
############################################################

write.csv(
  performance_table,
  "SVM_performance.csv",
  row.names = FALSE
)


############################################################
# 74. SAVE MODEL INFORMATION
############################################################

write.csv(
  model_information,
  "SVM_model_information.csv",
  row.names = FALSE
)


############################################################
# 75. SAVE THE STANDARDIZED FORMULA
############################################################

formula_file <- file(
  "SVM_formula.txt",
  open = "wt"
)

writeLines(
  "============================================================",
  formula_file
)

writeLines(
  "LINEAR SVM REGRESSION FORMULA",
  formula_file
)

writeLines(
  "============================================================",
  formula_file
)

writeLines(
  "",
  formula_file
)

writeLines(
  "STANDARDIZED FORM:",
  formula_file
)

writeLines(
  standardized_formula,
  formula_file
)

writeLines(
  "",
  formula_file
)

writeLines(
  "ORIGINAL SNP SCALE:",
  formula_file
)

writeLines(
  original_formula,
  formula_file
)

writeLines(
  "",
  formula_file
)

writeLines(
  "============================================================",
  formula_file
)

writeLines(
  "INTERCEPTS",
  formula_file
)

writeLines(
  "============================================================",
  formula_file
)

writeLines(
  paste0(
    "Standardized intercept = ",
    reference_prediction
  ),
  formula_file
)

writeLines(
  paste0(
    "Original-scale intercept = ",
    intercept_original
  ),
  formula_file
)

writeLines(
  "",
  formula_file
)

writeLines(
  "============================================================",
  formula_file
)

writeLines(
  "FORMULA RECONSTRUCTION CHECK",
  formula_file
)

writeLines(
  "============================================================",
  formula_file
)

writeLines(
  paste0(
    "Maximum absolute difference, training = ",
    max_train_difference
  ),
  formula_file
)

writeLines(
  paste0(
    "Maximum absolute difference, testing = ",
    max_test_difference
  ),
  formula_file
)

writeLines(
  paste0(
    "Formula reproduces training predictions = ",
    formula_matches_train
  ),
  formula_file
)

writeLines(
  paste0(
    "Formula reproduces testing predictions = ",
    formula_matches_test
  ),
  formula_file
)

close(
  formula_file
)


############################################################
# 76. SAVE A DETAILED COEFFICIENT TEXT FILE
############################################################

coefficient_file <- file(
  "SVM_detailed_SNP_coefficients.txt",
  open = "wt"
)

writeLines(
  "============================================================",
  coefficient_file
)

writeLines(
  "DETAILED LINEAR SVM SNP COEFFICIENTS",
  coefficient_file
)

writeLines(
  "============================================================",
  coefficient_file
)

writeLines(
  "",
  coefficient_file
)

for (j in seq_along(beta_standardized)) {
  
  snp_name <- colnames(X_train)[j]
  
  writeLines(
    paste0(
      "SNP: ",
      snp_name
    ),
    coefficient_file
  )
  
  writeLines(
    paste0(
      "  Mean = ",
      preprocess_mean[j]
    ),
    coefficient_file
  )
  
  writeLines(
    paste0(
      "  SD = ",
      preprocess_sd[j]
    ),
    coefficient_file
  )
  
  writeLines(
    paste0(
      "  Standardized coefficient = ",
      beta_standardized[j]
    ),
    coefficient_file
  )
  
  writeLines(
    paste0(
      "  Original-scale coefficient = ",
      beta_original[j]
    ),
    coefficient_file
  )
  
  writeLines(
    "",
    coefficient_file
  )
}

close(
  coefficient_file
)


############################################################
# 77. SAVE TRAINED MODEL
############################################################

saveRDS(
  svm_model,
  "SVM_linear_model.rds"
)


############################################################
# 78. SAVE THE FINAL KERNLAB MODEL
############################################################

saveRDS(
  final_model,
  "SVM_linear_final_kernlab_model.rds"
)


############################################################
# 79. SAVE COMPLETE ANALYSIS OBJECT
############################################################

analysis_results <- list(
  
  model = svm_model,
  
  final_model = final_model,
  
  best_tune = svm_model$bestTune,
  
  selected_snps = colnames(X_train),
  
  preprocessing_mean = preprocess_mean,
  
  preprocessing_sd = preprocess_sd,
  
  standardized_coefficients =
    beta_standardized,
  
  original_coefficients =
    beta_original,
  
  standardized_intercept =
    reference_prediction,
  
  original_intercept =
    intercept_original,
  
  standardized_formula =
    standardized_formula,
  
  original_formula =
    original_formula,
  
  train_predictions =
    pred_train,
  
  test_predictions =
    pred_test,
  
  train_formula_predictions =
    pred_train_formula,
  
  test_formula_predictions =
    pred_test_formula,
  
  train_R2 =
    train_R2,
  
  test_R2 =
    test_R2,
  
  train_RMSE =
    train_RMSE,
  
  test_RMSE =
    test_RMSE,
  
  train_MAE =
    train_MAE,
  
  test_MAE =
    test_MAE,
  
  coefficient_table =
    coefficient_table,
  
  train_contributions =
    train_contribution_table,
  
  test_contributions =
    test_contribution_table
)


saveRDS(
  analysis_results,
  "SVM_complete_analysis_results.rds"
)


############################################################
# 80. OPTIONAL: PREDICTED VS OBSERVED PLOT
############################################################

train_plot_data <- data.frame(
  Observed = y_train,
  Predicted = pred_train,
  Dataset = "Training"
)

test_plot_data <- data.frame(
  Observed = y_test,
  Predicted = pred_test,
  Dataset = "Testing"
)


plot_data <- rbind(
  train_plot_data,
  test_plot_data
)


p <- ggplot(
  plot_data,
  aes(
    x = Observed,
    y = Predicted
  )
) +
  
  geom_point() +
  
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  
  facet_wrap(
    ~Dataset
  ) +
  
  labs(
    title = "Observed vs Predicted FT10",
    x = "Observed FT10",
    y = "Predicted FT10"
  ) +
  
  theme_bw()


ggsave(
  "SVM_observed_vs_predicted.png",
  plot = p,
  width = 10,
  height = 5,
  dpi = 300
)


############################################################
# 81. PRINT TOP SNPs BY ABSOLUTE EFFECT
############################################################

cat("\n========================================\n")
cat("TOP SNPs BY ABSOLUTE STANDARDIZED EFFECT\n")
cat("========================================\n\n")


print(
  head(
    coefficient_table_sorted,
    21
  ),
  row.names = FALSE
)


############################################################
# 82. FINAL SUMMARY
############################################################

cat("\n\n")
cat("============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Number of training samples:",
  nrow(X_train),
  "\n"
)

cat(
  "Number of testing samples:",
  nrow(X_test),
  "\n"
)

cat(
  "Number of SNPs:",
  ncol(X_train),
  "\n"
)

cat(
  "Best SVM cost:",
  svm_model$bestTune$cost,
  "\n\n"
)


cat(
  "Training R2:",
  train_R2,
  "\n"
)

cat(
  "Testing R2:",
  test_R2,
  "\n"
)

cat(
  "Training RMSE:",
  train_RMSE,
  "\n"
)

cat(
  "Testing RMSE:",
  test_RMSE,
  "\n"
)

cat(
  "Training MAE:",
  train_MAE,
  "\n"
)

cat(
  "Testing MAE:",
  test_MAE,
  "\n\n"
)


cat(
  "Standardized intercept:",
  reference_prediction,
  "\n"
)

cat(
  "Original-scale intercept:",
  intercept_original,
  "\n\n"
)


cat(
  "Formula reconstruction max error - TRAIN:",
  max_train_difference,
  "\n"
)

cat(
  "Formula reconstruction max error - TEST:",
  max_test_difference,
  "\n\n"
)


cat(
  "Formula reproduces TRAIN predictions:",
  formula_matches_train,
  "\n"
)

cat(
  "Formula reproduces TEST predictions:",
  formula_matches_test,
  "\n\n"
)


cat(
  "Formula saved to: SVM_formula.txt\n"
)

cat(
  "Coefficients saved to: SVM_SNP_coefficients.csv\n"
)

cat(
  "Individual contributions saved to: SVM_test_SNP_contributions.csv\n"
)

cat(
  "Complete model saved to: SVM_complete_analysis_results.rds\n"
)

cat(
  "Prediction table saved to: SVM_test_predictions.csv\n"
)

cat(
  "Plot saved to: SVM_observed_vs_predicted.png\n"
)


############################################################
# 83. STOP PARALLEL CLUSTER
############################################################

stopCluster(cl)

registerDoSEQ()

##
## PARTE RELATIVA ALLA SYMBOLIC REGRESSION
##
set.seed(123)
library(gramEvol)
library(data.table)
library(foreach)
library(doParallel)
library(parallel)

# ------------------------------------------------------------
# Set up parallel processing
# ------------------------------------------------------------
# Detect number of cores
num_cores <- detectCores()

# Register parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

setwd(getwd())
ID_train_test <- read.delim("train_test_split.txt")
input <- fread("A_thaliana_TF10_NO_NA.txt")
feature_selected <- read.delim("selected_snps.txt")
input <- as.data.frame(input)


# rebuild the original training and testing set
keep_cols <- c("IID", "FT10", feature_selected$SNP)
input_selected <- input[, colnames(input) %in% keep_cols, drop = FALSE]

training_list <- ID_train_test[ID_train_test$set == "train",]
testing_list <- ID_train_test[ID_train_test$set == "test",]
training_list <- training_list[,-2]
testing_list <- testing_list[,-2]

train_df <- input_selected[input_selected$IID %in% training_list, ]
test_df  <- input_selected[input_selected$IID %in% testing_list, ]

y_train <- train_df$FT10
X_train <- train_df[, !(names(train_df) %in% c("IID", "FT10"))]

y_test <- test_df$FT10
X_test <- test_df[, !(names(test_df) %in% c("IID", "FT10"))]
X_train_scaled <- scale(X_train)
X_test_scaled  <- scale(X_test,
                        center = attr(X_train_scaled, "scaled:center"),
                        scale  = attr(X_train_scaled, "scaled:scale"))

lm_data  <- cbind(phenotype = y_train, X_train_scaled)
test_data <- cbind(phenotype = y_test, X_test_scaled)

# Extract SNP names for grammar
snp_names <- colnames(X_train)

# ------------------------------------------------------------
# 1. Prepare data for symbolic regression
# ------------------------------------------------------------
marker_names <- snp_names

# Create training data with phenotype column (rename FT10 to phenotype)
lm_data <- data.frame(phenotype = y_train, X_train_scaled, check.names = FALSE)

# Create test data with phenotype column
test_data <- data.frame(phenotype = y_test, X_test_scaled, check.names = FALSE)

# ------------------------------------------------------------
# 2. Define safe mathematical functions
# ------------------------------------------------------------
safe_log <- function(x) ifelse(x <= 0, 0, log(x))
safe_sqrt <- function(x) ifelse(x < 0, 0, sqrt(x))
safe_div <- function(x, y) ifelse(y == 0, 0, x / y)

# ------------------------------------------------------------
# 3. Define grammar
# ------------------------------------------------------------
ruleDef <- list(
  expr = gsrule("<expr><op><expr>", "<func>(<expr>)", "<var>", "<const>"),
  op = gsrule("+", "-", "*"),
  func = grule(safe_log, exp, safe_sqrt),
  var = do.call(grule, as.list(marker_names)),  # dynamically include all markers
  const = grule(1, 2, 3)
)

grammarDef <- CreateGrammar(ruleDef)

# ------------------------------------------------------------
# 4. Define evaluation function with regularization and quote cleaning
# ------------------------------------------------------------
evalFunc <- function(expr) {
  tryCatch({
    # Clean expression from quotes
    expr_clean <- gsub('[\"\'`]', '', expr)
    
    # Evaluate expression
    e <- parse(text = expr_clean)
    y_pred <- eval(e, envir = lm_data)
    y <- lm_data$phenotype
    
    # Validate predictions
    if (length(y_pred) != length(y) ||
        !is.numeric(y_pred) ||
        any(is.na(y_pred)) ||
        any(is.infinite(y_pred))) {
      return(1e10)
    }
    
    # RMSE
    rmse <- sqrt(mean((y - y_pred)^2))
    
    # Complexity penalty
    n_terms <- length(strsplit(expr_clean, split = "\\+|-|\\*")[[1]])
    n_ops <- length(gregexpr("\\+|-|\\*", expr_clean)[[1]])
    if (n_ops < 0) n_ops <- 0
    
    lambda <- 0.01  # regularization strength
    cost <- rmse + lambda * (n_terms + n_ops)
    
    return(cost)
    
  }, error = function(e) {
    return(1e10)
  })
}

# ------------------------------------------------------------
# 5. Define metrics computation function
# ------------------------------------------------------------
compute_metrics <- function(expr, data) {
  tryCatch({
    y <- data$phenotype
    
    # Clean expression from quotes
    expr_clean <- gsub('[\"\'`]', '', expr)
    
    e <- parse(text = expr_clean)
    y_pred <- eval(e, envir = data)
    
    if (length(y_pred) != length(y) || 
        !is.numeric(y_pred) || 
        any(is.na(y_pred)) ||
        any(is.infinite(y_pred))) {
      return(list(
        RMSE = NA,
        MAE = NA,
        R_squared = NA,
        Pearson_r_squared = NA
      ))
    }
    
    rmse <- sqrt(mean((y - y_pred)^2))
    mae <- mean(abs(y - y_pred))
    ss_res <- sum((y - y_pred)^2)
    ss_tot <- sum((y - mean(y))^2)
    r_squared <- ifelse(ss_tot == 0, 0, 1 - (ss_res / ss_tot))
    pearson_r <- cor(y, y_pred, use = "complete.obs")
    pearson_r_squared <- ifelse(is.na(pearson_r), 0, pearson_r^2)
    
    list(
      RMSE = rmse,
      MAE = mae,
      R_squared = r_squared,
      Pearson_r_squared = pearson_r_squared
    )
  }, error = function(e) {
    return(list(
      RMSE = NA,
      MAE = NA,
      R_squared = NA,
      Pearson_r_squared = NA
    ))
  })
}

# ------------------------------------------------------------
# 6. Run parallel symbolic regression
# ------------------------------------------------------------
cat("Starting parallel symbolic regression with", num_cores, "cores...\n")

# Number of runs
n_runs <- 1000

# Run parallel symbolic regression
results_list <- foreach(i = 1:n_runs, 
                        .packages = c("gramEvol", "data.table"),
                        .combine = rbind,
                        .errorhandling = "remove") %dopar% {
                          cat("Starting run", i, "of", n_runs, "\n")
                          set.seed(i)
                          
                          # Run grammatical evolution
                          ge <- GrammaticalEvolution(
                            grammarDef,
                            evalFunc,
                            iterations = 1000,
                            popSize = 300,
                            newPerGen = 15,
                            elitism = 3,
                            mutationChance = 0.3,
                            optimizer = "ga",
                            verbose = FALSE
                          )
                          
                          # Extract best expression
                          best_expr <- as.character(ge$best$expressions)
                          best_expr <- gsub('[\"\'`]', '', best_expr)  # Clean quotes
                          
                          # Compute metrics on training set
                          train_metrics <- compute_metrics(best_expr, lm_data)
                          
                          # Compute metrics on test set
                          test_metrics <- compute_metrics(best_expr, test_data)
                          
                          # Return results as a data frame row
                          data.frame(
                            Seed = paste0("Seed N:", i),
                            Train_R_squared = train_metrics$R_squared,
                            Train_Pearson_r_squared = train_metrics$Pearson_r_squared,
                            Train_RMSE = train_metrics$RMSE,
                            Train_MAE = train_metrics$MAE,
                            Test_R_squared = test_metrics$R_squared,
                            Test_Pearson_r_squared = test_metrics$Pearson_r_squared,
                            Test_RMSE = test_metrics$RMSE,
                            Test_MAE = test_metrics$MAE,
                            Best_expression = best_expr,
                            stringsAsFactors = FALSE
                          )
                        }

# Stop the parallel cluster
stopCluster(cl)

# ------------------------------------------------------------
# 7. Process and summarize results
# ------------------------------------------------------------
# Convert to data frame if not already
if (is.list(results_list) && !is.data.frame(results_list)) {
  results_table <- do.call(rbind, results_list)
} else {
  results_table <- results_list
}

# Find best overall expression based on test R-squared
valid_results <- results_table[!is.na(results_table$Test_R_squared), ]
if (nrow(valid_results) > 0) {
  best_idx <- which.max(valid_results$Test_R_squared)
  best_overall_expr <- valid_results$Best_expression[best_idx]
  best_train_metrics <- list(
    R_squared = valid_results$Train_R_squared[best_idx],
    Pearson_r_squared = valid_results$Train_Pearson_r_squared[best_idx],
    RMSE = valid_results$Train_RMSE[best_idx],
    MAE = valid_results$Train_MAE[best_idx]
  )
  best_test_metrics <- list(
    R_squared = valid_results$Test_R_squared[best_idx],
    Pearson_r_squared = valid_results$Test_Pearson_r_squared[best_idx],
    RMSE = valid_results$Test_RMSE[best_idx],
    MAE = valid_results$Test_MAE[best_idx]
  )
} else {
  best_overall_expr <- "No valid expressions found"
  best_train_metrics <- list(R_squared = NA, Pearson_r_squared = NA, RMSE = NA, MAE = NA)
  best_test_metrics <- list(R_squared = NA, Pearson_r_squared = NA, RMSE = NA, MAE = NA)
}

# ------------------------------------------------------------
# 8. Print summary statistics
# ------------------------------------------------------------
cat("\n=== BEST OVERALL EXPRESSION ===\n")
cat("Expression:", best_overall_expr, "\n\n")
cat("=== TRAINING SET METRICS (Best Expression) ===\n")
cat("R-squared:", round(best_train_metrics$R_squared, 4), "\n")
cat("Pearson r-squared:", round(best_train_metrics$Pearson_r_squared, 4), "\n")
cat("MAE:", round(best_train_metrics$MAE, 4), "\n\n")



cat("=== TEST SET METRICS (Best Expression) ===\n")
cat("R-squared:", round(best_test_metrics$R_squared, 4), "\n")
cat("Pearson r-squared:", round(best_test_metrics$Pearson_r_squared, 4), "\n")
cat("RMSE:", round(best_test_metrics$RMSE, 4), "\n")
cat("MAE:", round(best_test_metrics$MAE, 4), "\n\n")



# ------------------------------------------------------------
# 9. Save full results to file
# ------------------------------------------------------------
# Convert expressions to character before saving
results_table$Best_expression <- as.character(results_table$Best_expression)

# Save detailed results
write.table(results_table, 
            file = "symbolic_regression_parallel_results.txt", 
            row.names = FALSE, 
            quote = FALSE, 
            sep = "\t")

# Save summary of best expression
best_summary <- data.frame(
  Metric = c("Best_Expression", 
             "Train_R_squared", "Train_Pearson_r_squared", "Train_RMSE", "Train_MAE",
             "Test_R_squared", "Test_Pearson_r_squared", "Test_RMSE", "Test_MAE"),
  Value = c(best_overall_expr,
            best_train_metrics$R_squared, best_train_metrics$Pearson_r_squared,
            best_train_metrics$RMSE, best_train_metrics$MAE,
            best_test_metrics$R_squared, best_test_metrics$Pearson_r_squared,
            best_test_metrics$RMSE, best_test_metrics$MAE)
)

write.table(best_summary, 
            file = "best_symbolic_model_summary.txt", 
            row.names = FALSE, 
            quote = FALSE, 
            sep = "\t")

cat("\nResults saved to:\n")

library(ggplot2)

# Read the results table
results_table <- read.delim(
  "symbolic_regression_parallel_results.txt",
  stringsAsFactors = FALSE
)

# Which metrics to plot
metrics_to_plot <- c(
  "Train_R_squared", "Test_R_squared",
  "Train_RMSE",      "Test_RMSE",
  "Train_MAE",       "Test_MAE"
)

# Long format for ggplot
results_long <- data.frame(
  Metric = factor(
    rep(metrics_to_plot, each = nrow(results_table)),
    levels = metrics_to_plot,
    labels = c("Train R²", "Test R²",
               "Train RMSE", "Test RMSE",
               "Train MAE",  "Test MAE")
  ),
  Value = unlist(results_table[, metrics_to_plot], use.names = FALSE)
)

# Violin + boxplot + mean diamond
p <- ggplot(results_long, aes(x = Metric, y = Value, fill = Metric)) +
  geom_violin(alpha = 0.55, trim = FALSE, colour = NA) +
  geom_boxplot(width = 0.18, alpha = 0.85, outlier.size = 0.8,
               colour = "grey20") +
  stat_summary(fun = mean, geom = "point", shape = 23,
               size = 2.6, fill = "white", colour = "black") +
  facet_wrap(~ Metric, scales = "free", ncol = 3) +
  theme_bw() +
  labs(
    title    = "Distribution of Metrics Across All Symbolic-Regression Runs",
    subtitle = "White diamond = mean | Boxplot = median & quartiles | Violin = density",
    x = NULL, y = "Value"
  ) +
  theme(
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5),
    strip.text      = element_text(face = "bold")
  )

ggsave("symbolic_regression_distribution.png", plot = p,
       width = 12, height = 8, dpi = 300)

# Summary stats
summary_stats <- data.frame(
  Metric = metrics_to_plot,
  Mean   = sapply(results_table[, metrics_to_plot], mean,   na.rm = TRUE),
  Median = sapply(results_table[, metrics_to_plot], median, na.rm = TRUE),
  SD     = sapply(results_table[, metrics_to_plot], sd,     na.rm = TRUE),
  Q1     = sapply(results_table[, metrics_to_plot],
                  quantile, probs = 0.25, na.rm = TRUE),
  Q3     = sapply(results_table[, metrics_to_plot],
                  quantile, probs = 0.75, na.rm = TRUE)
)
summary_stats$IQR <- summary_stats$Q3 - summary_stats$Q1

write.csv(summary_stats,
          "symbolic_regression_distribution_summary.csv",
          row.names = FALSE)

print(summary_stats, row.names = FALSE)


# FINE PARTE FENOTIPO CONTINUO
# INIZIO PARTE FENOTIPO BINARIO

###
##FEATURE SELECTION
###

setwd("C:\\Users\\Salvatore\\Desktop\\Dottorato\\richiesta_dati_reali\\EGAD00000000003_decriptati\\effettivamente_utili")

# 1. Load required libraries
library(data.table)
library(nestedcv)
library(glmnet)
library(pROC)

# ====================================================================
# 2. Load data
# ====================================================================
setwd("C:\\Users\\Salvatore\\Desktop\\Dottorato\\richiesta_dati_reali\\EGAD00000000003_decriptati\\effettivamente_utili")

final_df <- fread("FINAL_shuffled_CASI_E_CONTROLLI_BD_58C_chr22.txt")

# ====================================================================
# 3. Extract phenotype and convert to numeric 0/1
# ====================================================================
# The 2nd column contains the phenotype (already 0/1 as per your check)
y <- as.numeric(final_df[[2]])   # numeric 0 and 1

# Double‑check that y is indeed 0/1
if (!all(y %in% c(0, 1))) {
  stop("Phenotype must contain only 0 and 1.")
}
cat("Phenotype distribution:\n")
print(table(y))

# ====================================================================
# 4. Train / Test split (80/20)
# ====================================================================
set.seed(123)
n <- nrow(final_df)
train_idx <- sample(seq_len(n), size = floor(0.8 * n))
test_idx  <- setdiff(seq_len(n), train_idx)

# Training and test data frames (full rows)
train_df <- final_df[train_idx, ]
test_df  <- final_df[test_idx, ]

# Separate y and X (SNPs start at column 3)
y_train <- y[train_idx]
y_test  <- y[test_idx]

X_train <- as.matrix(train_df[, 3:ncol(train_df), with = FALSE])
X_test  <- as.matrix(test_df[,  3:ncol(test_df),  with = FALSE])

# Force numeric (just in case)
X_train <- apply(X_train, 2, as.numeric)
X_test  <- apply(X_test,  2, as.numeric)

# ====================================================================
# 5. Remove zero‑variance SNPs (and any constant columns)
# ====================================================================
train_var <- apply(X_train, 2, var, na.rm = TRUE)
zero_var   <- which(is.na(train_var) | train_var == 0)

if (length(zero_var) > 0) {
  cat("Removing", length(zero_var), "zero‑variance SNPs.\n")
  X_train <- X_train[, -zero_var, drop = FALSE]
  X_test  <- X_test[,  -zero_var, drop = FALSE]
}

cat("Number of SNPs retained:", ncol(X_train), "\n")

# ====================================================================
# 6. Check for missing/infinite values
# ====================================================================
if (any(!is.finite(X_train))) stop("X_train contains NA/NaN/Inf.")
if (any(!is.finite(X_test)))  stop("X_test contains NA/NaN/Inf.")

# ====================================================================
# 7. Nested CV LASSO (binomial) – feature selection
# ====================================================================
# IMPORTANT: for binomial family, y must be a factor
y_train_factor <- factor(y_train, levels = c(0, 1))
y_test_factor  <- factor(y_test,  levels = c(0, 1))

set.seed(123)
fit <- nestcv.glmnet(
  x = X_train,
  y = y_train_factor,               # factor with levels 0,1
  family = "binomial",
  alpha = 1,
  outer_k = 10,
  inner_k = 10,
  min_1se = 1,
  type.measure = "auc"              # use AUC for inner CV; can also use "deviance"
)

# ====================================================================
# 8. Extract selected SNPs (coefficients != 0)
# ====================================================================
coef_mat <- fit$final_coef
features_selected <- rownames(coef_mat)[coef_mat[, 1] != 0]
features_selected <- setdiff(features_selected, "(Intercept)")

cat("\nNumber of selected features:", length(features_selected), "\n")

# ====================================================================
# 9. Build reduced datasets using only selected SNPs
# ====================================================================
X_train_fs <- X_train[, features_selected, drop = FALSE]
X_test_fs  <- X_test[,  features_selected, drop = FALSE]

# Create data frames for logistic regression
df_train_fs <- as.data.frame(X_train_fs)
df_train_fs$y <- y_train

df_test_fs <- as.data.frame(X_test_fs)
df_test_fs$y <- y_test

# ====================================================================
# SAVE REDUCED DATASETS AND SELECTED SNPS
# ====================================================================

# 1. Save reduced data frames (selected SNPs + phenotype)
# ------------------------------------------
# CSV format (human‑readable, portable)
write.csv(df_train_fs, file = "train_selected_SNPs.csv", row.names = FALSE)
write.csv(df_test_fs,  file = "test_selected_SNPs.csv",  row.names = FALSE)

# Optional: also save as RData to preserve column classes
# save(df_train_fs, df_test_fs, file = "selected_datasets.RData")

# 2. Save the list of selected SNP names
# ------------------------------------------
# Plain text (one per line)
writeLines(features_selected, con = "selected_SNPs.txt")

# R binary format (for easy reloading in R)
saveRDS(features_selected, file = "selected_SNPs.rds")

# 3. (Optional) Save train/test indices for reproducibility
# ------------------------------------------
saveRDS(list(train_idx = train_idx, test_idx = test_idx), 
        file = "train_test_indices.rds")

cat("\nAll files saved successfully.\n")


##
## LOGISTIC REGRESSION
##
library(pROC)

features_selected <- read.delim("selected_SNPs.txt", header = F,
                                stringsAsFactors = FALSE)
df_train_fs <- read.csv("train_selected_SNPs.csv")
df_test_fs <- read.csv("test_selected_SNPs.csv")



# Check number of selected SNPs
cat("Number of selected SNPs:", length(features_selected), "\n")

if (length(features_selected) == 0) {
  stop("No SNPs were selected by LASSO.")
}


# --------------------------------------------------------------------
# Create logistic regression formula
# --------------------------------------------------------------------
features_selected <- features_selected$V1

logistic_formula <- reformulate(
  features_selected,
  response = "y"
)


cat("\nLogistic regression formula:\n")


# --------------------------------------------------------------------
# Fit logistic regression on TRAINING data only
# --------------------------------------------------------------------

logistic_model <- glm(
  formula = logistic_formula,
  data = df_train_fs,
  family = binomial(link = "logit")
)

cat("\nLogistic regression fitted successfully.\n")

# ====================================================================
# 10B. EXTRACT LOGISTIC REGRESSION MATHEMATICAL FORMULA
# ====================================================================

cat("\n============================================\n")


cat("LOGISTIC REGRESSION MATHEMATICAL FORMULA\n")
cat("============================================\n")


# Extract coefficients
coef_logistic <- coef(logistic_model)

# Intercept
intercept <- unname(coef_logistic["(Intercept)"])

# SNP coefficients
beta <- coef_logistic[names(coef_logistic) != "(Intercept)"]


# --------------------------------------------------------------------
# Create the linear predictor formula
# --------------------------------------------------------------------

formula_terms <- paste0(
  sprintf("%.10f", beta),
  " * ",
  names(beta)
)

linear_formula <- paste(
  c(sprintf("%.10f", intercept), formula_terms),
  collapse = " + "
)

# Replace "+ -" with "- " for negative coefficients
linear_formula <- gsub(
  "\\+ -",
  "- ",
  linear_formula
)

cat("\nLinear predictor (log-odds) formula:\n\n")


cat("eta = ", linear_formula, "\n", sep = "")


# --------------------------------------------------------------------
# Create probability formula
# --------------------------------------------------------------------

cat("\nProbability formula:\n\n")



cat(
  "p = 1 / (1 + exp(-eta))\n"
)

cat("\nwhere eta is:\n\n")


cat(
  "eta = ",
  linear_formula,
  "\n",
  sep = ""
)


# --------------------------------------------------------------------
# Save coefficients
# --------------------------------------------------------------------

logistic_coefficients <- data.frame(
  Variable = names(coef_logistic),
  Coefficient = as.numeric(coef_logistic),
  Odds_Ratio = exp(as.numeric(coef_logistic))
)

write.csv(
  logistic_coefficients,
  "logistic_regression_coefficients.csv",
  row.names = FALSE
)


# ====================================================================
# 11. PREDICT PROBABILITIES
# ====================================================================

# IMPORTANT:
# type = "response" gives probabilities between 0 and 1

p_train <- predict(
  logistic_model,
  newdata = df_train_fs,
  type = "response"
)

p_test <- predict(
  logistic_model,
  newdata = df_test_fs,
  type = "response"
)


# Check predictions
cat("\nTraining predicted probabilities:\n")

print(head(p_train))


cat("\nTesting predicted probabilities:\n")

print(head(p_test))


# ====================================================================
# 12. AUC-ROC
# ====================================================================

roc_train <- roc(
  response = df_train_fs$y,
  predictor = p_train,
  quiet = TRUE
)

roc_test <- roc(
  response = df_test_fs$y,
  predictor = p_test,
  quiet = TRUE
)

auc_train <- as.numeric(auc(roc_train))
auc_test  <- as.numeric(auc(roc_test))


# ====================================================================
# 13. BRIER SCORE
# ====================================================================

# Brier score:
#
#       1
# BS = --- sum((y - p)^2)
#       n

brier_score <- function(obs, pred) {
  
  mean((obs - pred)^2)
}

brier_train <- brier_score(
  df_train_fs$y,
  p_train
)

brier_test <- brier_score(
  df_test_fs$y,
  p_test
)


# ====================================================================
# 14. LOG LOSS
# ====================================================================

# Log loss:
#
#        1
# LL = - - sum[y log(p) + (1-y) log(1-p)]
#        n

log_loss <- function(obs, pred) {
  
  # Prevent log(0)
  eps <- 1e-15
  
  pred <- pmin(
    pmax(pred, eps),
    1 - eps
  )
  
  -mean(
    obs * log(pred) +
      (1 - obs) * log(1 - pred)
  )
}

logloss_train <- log_loss(
  df_train_fs$y,
  p_train
)

logloss_test <- log_loss(
  df_test_fs$y,
  p_test
)


# ====================================================================
# 15. COMBINE PERFORMANCE RESULTS
# ====================================================================

performance <- data.frame(
  
  Dataset = c(
    "Training",
    "Testing"
  ),
  
  AUC_ROC = c(
    auc_train,
    auc_test
  ),
  
  Brier_Score = c(
    brier_train,
    brier_test
  ),
  
  Log_Loss = c(
    logloss_train,
    logloss_test
  )
)


# Print results
cat("\n============================================\n")


cat("LOGISTIC REGRESSION PERFORMANCE\n")
cat("============================================\n")

print(performance)

##
## SYMBOLIC REGRESSION
##

library(gramEvol)
library(data.table)
library(parallel)
library(pbapply)
library(pROC)

# ------------------------------------------------------------
# 1. Imposta numero di thread e cluster parallelo
# ------------------------------------------------------------
total_threads <- parallel::detectCores(logical = TRUE)
num_cores <- total_threads - 1

# One data.table thread per worker to avoid oversubscription
setDTthreads(1)

print(paste0(
  "Available logical threads: ", total_threads,
  " | Parallel workers: ", num_cores
))

print("Creating cluster...")
cl <- makeCluster(num_cores)
print("Cluster created successfully.")

# ------------------------------------------------------------
# 2. Carica dataset
# ------------------------------------------------------------
print("Loading datasets...")

features_selected <- fread("selected_SNPs.txt", header = FALSE)
marker_names <- features_selected$V1

df_train_fs <- fread("train_selected_SNPs_header_fenotipo_corretto_in_phenotype_da_y.csv")
df_test_fs  <- fread("test_selected_SNPs_header_fenotipo_corretto_in_phenotype_da_y.csv")

# Remove any literal quote characters left in the column names
names(df_train_fs) <- gsub('"', '', names(df_train_fs))
names(df_test_fs)  <- gsub('"', '', names(df_test_fs))

print("Datasets loaded successfully.")

# Check phenotype
if (!"phenotype" %in% names(df_train_fs)) {
  stop("Column 'phenotype' not found in train_selected_SNPs.csv")
}

if (!"phenotype" %in% names(df_test_fs)) {
  stop("Column 'phenotype' not found in test_selected_SNPs.csv")
}

# Keep only selected SNPs + phenotype
missing_train <- setdiff(marker_names, names(df_train_fs))
missing_test  <- setdiff(marker_names, names(df_test_fs))

if (length(missing_train) > 0) {
  stop("The following selected SNPs are missing from the training dataset: ",
       paste(missing_train, collapse = ", "))
}

if (length(missing_test) > 0) {
  stop("The following selected SNPs are missing from the test dataset: ",
       paste(missing_test, collapse = ", "))
}

lm_data <- df_train_fs[, c(marker_names, "phenotype"), with = FALSE]
test_data <- df_test_fs[, c(marker_names, "phenotype"), with = FALSE]

# Make sure phenotype is binary numeric 0/1
lm_data[, phenotype := as.numeric(phenotype)]
test_data[, phenotype := as.numeric(phenotype)]

if (!all(na.omit(lm_data$phenotype) %in% c(0, 1))) {
  stop("Training phenotype must contain only 0/1 values.")
}

if (!all(na.omit(test_data$phenotype) %in% c(0, 1))) {
  stop("Test phenotype must contain only 0/1 values.")
}

print(paste0("Training samples: ", nrow(lm_data)))
print(paste0("Test samples: ", nrow(test_data)))
print(paste0("Selected SNPs: ", length(marker_names)))

# ------------------------------------------------------------
# 3. Funzione principale per ogni seed
# ------------------------------------------------------------
run_single_seed <- function(i, lm_data, test_data, marker_names) {
  
  library(gramEvol)
  library(data.table)
  
  print(paste0("Running seed ", i))
  set.seed(i)
  
  # --- Funzioni sicure ---
  safe_log <- function(x) ifelse(x <= 0, 0, log(x))
  safe_sqrt <- function(x) ifelse(x < 0, 0, sqrt(x))
  safe_div <- function(x, y) ifelse(y == 0, 0, x / y)
  
  # --- Crea grammatica ---
  ruleDef <- list(
    expr = gsrule(
      "<expr><op><expr>",
      "<func>(<expr>)",
      "<var>",
      "<const>"
    ),
    op = gsrule("+", "-", "*", "/"),
    func = grule(sin, safe_log, exp, safe_sqrt),
    var = do.call(grule, as.list(marker_names)),
    const = grule(1, 2, 3)
  )
  
  grammarDef <- CreateGrammar(ruleDef)
  
  # --- Funzione di valutazione ---
  evalFunc <- function(expr) {
    
    tryCatch({
      
      for (m in marker_names) {
        expr <- gsub(paste0('"', m, '"'), m, expr, fixed = TRUE)
        expr <- gsub(paste0("'", m, "'"), m, expr, fixed = TRUE)
      }
      
      e <- parse(text = expr)
      score <- eval(e, envir = lm_data)
      y <- lm_data$phenotype
      
      if (length(score) != length(y) ||
          !is.numeric(score) ||
          any(!is.finite(score))) {
        return(1e10)
      }
      
      # Convert symbolic-regression score to probability
      prob <- plogis(score)
      
      if (any(!is.finite(prob))) {
        return(1e10)
      }
      
      # Log Loss
      eps <- 1e-15
      prob <- pmin(pmax(prob, eps), 1 - eps)
      
      logloss <- -mean(
        y * log(prob) + (1 - y) * log(1 - prob)
      )
      
      # Penalità per complessità
      n_terms <- length(strsplit(expr, split = "\\+|-|\\*")[[1]])
      n_ops <- length(gregexpr("\\+|-|\\*", expr)[[1]])
      
      if (n_ops < 0) n_ops <- 0
      
      lambda <- 3
      cost <- logloss + lambda * (n_terms + n_ops)
      
      return(cost)
      
    }, error = function(e) {
      return(1e10)
    })
  }
  
  # --- Esecuzione grammaticale ---
  ge <- GrammaticalEvolution(
    grammarDef,
    evalFunc,
    iterations = 1000,
    popSize = 300,
    newPerGen = 15,
    elitism = 3,
    mutationChance = 0.3,
    optimizer = "ga"
  )
  
  best_expr <- as.character(ge$best$expressions)
  
  # ----------------------------------------------------------
  # Funzione per calcolare metriche
  # ----------------------------------------------------------
  compute_metrics <- function(expr, data) {
    
    tryCatch({
      
      y <- data$phenotype
      
      for (m in marker_names) {
        expr <- gsub(paste0('"', m, '"'), m, expr, fixed = TRUE)
        expr <- gsub(paste0("'", m, "'"), m, expr, fixed = TRUE)
      }
      
      e <- parse(text = expr)
      score <- eval(e, envir = data)
      
      if (length(score) != length(y) ||
          !is.numeric(score) ||
          any(!is.finite(score))) {
        
        return(list(
          AUC_ROC = NA,
          Brier_Score = NA,
          Log_Loss = NA
        ))
      }
      
      # Convert symbolic score to probability
      prob <- plogis(score)
      
      if (any(!is.finite(prob))) {
        return(list(
          AUC_ROC = NA,
          Brier_Score = NA,
          Log_Loss = NA
        ))
      }
      
      # Keep probabilities away from exactly 0/1
      eps <- 1e-15
      prob <- pmin(pmax(prob, eps), 1 - eps)
      
      # AUC ROC
      auc_value <- tryCatch(
        as.numeric(
          pROC::auc(
            response = y,
            predictor = prob,
            direction = "<"
          )
        ),
        error = function(e) NA_real_
      )
      
      # Brier Score
      brier_score <- mean((prob - y)^2)
      
      # Log Loss
      log_loss <- -mean(
        y * log(prob) + (1 - y) * log(1 - prob)
      )
      
      return(list(
        AUC_ROC = auc_value,
        Brier_Score = brier_score,
        Log_Loss = log_loss
      ))
      
    }, error = function(e) {
      
      return(list(
        AUC_ROC = NA,
        Brier_Score = NA,
        Log_Loss = NA
      ))
    })
  }
  
  # ----------------------------------------------------------
  # Calcolo metriche TRAIN e TEST
  # ----------------------------------------------------------
  train_metrics <- compute_metrics(best_expr, lm_data)
  test_metrics  <- compute_metrics(best_expr, test_data)
  
  return(list(
    Seed = paste0("Seed N:", i),
    
    Train_AUC_ROC = train_metrics$AUC_ROC,
    Train_Brier_Score = train_metrics$Brier_Score,
    Train_Log_Loss = train_metrics$Log_Loss,
    
    Test_AUC_ROC = test_metrics$AUC_ROC,
    Test_Brier_Score = test_metrics$Brier_Score,
    Test_Log_Loss = test_metrics$Log_Loss,
    
    Best_expression = best_expr
  ))
}

# ------------------------------------------------------------
# 4. Esporta funzione e variabili ai worker
# ------------------------------------------------------------
clusterExport(
  cl,
  varlist = c(
    "run_single_seed",
    "lm_data",
    "test_data",
    "marker_names"
  )
)

# ------------------------------------------------------------
# 5. Esegui in parallelo con barra di avanzamento
# ------------------------------------------------------------
seeds_to_run <- 1000:1099

print("Starting parallel runs with progress bar...")

results_list <- pblapply(
  seeds_to_run,
  function(i) {
    run_single_seed(
      i,
      lm_data,
      test_data,
      marker_names
    )
  },
  cl = cl
)

print("All parallel runs completed.")

# ------------------------------------------------------------
# 6. Salva risultati
# ------------------------------------------------------------
results_table <- rbindlist(results_list, fill = TRUE)

print("Saving results...")

fwrite(
  results_table,
  "symbolic_results_summary_binary_1_100.csv"
)

print("Results saved successfully.")

# ------------------------------------------------------------
# 7. Chiudi cluster
# ------------------------------------------------------------
stopCluster(cl)

print("Cluster stopped. All done.")

# ====================================================================
# SYMBOLIC REGRESSION (BINARY) - DISTRIBUTION VIOLIN PLOT
# Option A: drop outliers with implausible metric values
# ====================================================================

setwd("C:\\Users\\Salvatore\\Desktop\\Dottorato\\richiesta_dati_reali\\EGAD00000000003_decriptati\\effettivamente_utili")
library(ggplot2)

# ------------------------------------------------------------
# 1. Load the results file and fix the shifted columns
# ------------------------------------------------------------
raw <- read.delim(
  "symbolic_results_summary_binary_1_1000.txt",
  stringsAsFactors = FALSE
)

results_table <- data.frame(
  Seed              = raw[[2]],
  Train_AUC_ROC     = as.numeric(raw[[3]]),
  Train_Brier_Score = as.numeric(raw[[4]]),
  Train_Log_Loss    = as.numeric(raw[[5]]),
  Test_AUC_ROC      = as.numeric(raw[[6]]),
  Test_Brier_Score  = as.numeric(raw[[7]]),
  Test_Log_Loss     = as.numeric(raw[[8]]),
  Best_expression   = raw[[9]],
  stringsAsFactors  = FALSE
)

# Drop rows where the key numeric metrics are NA
results_table <- results_table[
  !is.na(results_table$Train_AUC_ROC) &
    !is.na(results_table$Test_AUC_ROC),
]

# Drop runs where grammatical evolution failed to produce a
# meaningful expression (empty or containing no SNP term).
valid_expr <- !is.na(results_table$Best_expression) &
  nchar(results_table$Best_expression) > 0 &
  grepl("SNP_", results_table$Best_expression)

cat("Rows with valid expressions:", sum(valid_expr), "of", nrow(results_table), "\n")
results_table <- results_table[valid_expr, ]

# ------------------------------------------------------------
# 2. Option A - keep only domain-plausible runs
# ------------------------------------------------------------
# A valid binary classifier must satisfy:
#   0.5 <= AUC <= 1
#   0.1 <= Brier <= 0.4
#   0.3 <= LogLoss <= 1.0
# Everything outside these ranges comes from degenerate symbolic
# expressions and is dropped before plotting.

metrics_to_plot <- c(
  "Train_AUC_ROC", "Test_AUC_ROC",
  "Train_Brier_Score", "Test_Brier_Score",
  "Train_Log_Loss", "Test_Log_Loss"
)

meaningful <- apply(results_table[, metrics_to_plot], 1, function(x) {
  all(!is.na(x)) &&
    x["Train_AUC_ROC"]     >= 0.5 && x["Train_AUC_ROC"]     <= 1.0 &&
    x["Test_AUC_ROC"]      >= 0.5 && x["Test_AUC_ROC"]      <= 1.0 &&
    x["Train_Brier_Score"] >= 0.1 && x["Train_Brier_Score"] <= 0.4 &&
    x["Test_Brier_Score"]  >= 0.1 && x["Test_Brier_Score"]  <= 0.4 &&
    x["Train_Log_Loss"]    >= 0.3 && x["Train_Log_Loss"]    <= 1.0 &&
    x["Test_Log_Loss"]     >= 0.3 && x["Test_Log_Loss"]     <= 1.0
})

cat("Rows before cleaning:", nrow(results_table), "\n")
results_table <- results_table[meaningful, ]
cat("Rows after cleaning: ", nrow(results_table), "\n")

# ------------------------------------------------------------
# 3. Sanity check
# ------------------------------------------------------------
cat("\nStructure of results_table:\n")
str(results_table)

cat("\nSummary of the six metrics (after cleaning):\n")
print(summary(results_table[, metrics_to_plot]))

# ------------------------------------------------------------
# 4. Reshape into (Family, Split) long format
# ------------------------------------------------------------
family_map <- c(
  Train_AUC_ROC      = "AUC",
  Test_AUC_ROC       = "AUC",
  Train_Brier_Score  = "Brier",
  Test_Brier_Score   = "Brier",
  Train_Log_Loss     = "LogLoss",
  Test_Log_Loss      = "LogLoss"
)

split_map <- c(
  Train_AUC_ROC      = "Train",
  Test_AUC_ROC       = "Test",
  Train_Brier_Score  = "Train",
  Test_Brier_Score   = "Test",
  Train_Log_Loss     = "Train",
  Test_Log_Loss      = "Test"
)

long <- do.call(rbind, lapply(metrics_to_plot, function(m) {
  data.frame(
    Family = family_map[[m]],
    Split  = split_map[[m]],
    Value  = results_table[[m]]
  )
}))

long$Family <- factor(long$Family, levels = c("AUC", "Brier", "LogLoss"))
long$Split  <- factor(long$Split,  levels = c("Train", "Test"))

# ------------------------------------------------------------
# 5. Build the violin plot
# ------------------------------------------------------------
p <- ggplot(long, aes(x = Split, y = Value, fill = Split)) +
  geom_violin(alpha = 0.55, trim = FALSE, colour = NA) +
  geom_boxplot(width = 0.15, alpha = 0.85, outlier.size = 0.8,
               colour = "grey20") +
  stat_summary(fun = mean, geom = "point", shape = 23,
               size = 2.6, fill = "white", colour = "black") +
  facet_wrap(~ Family, scales = "free_y", ncol = 3) +
  theme_bw() +
  labs(
    title    = "Distribution of Metrics - Train vs Test",
    subtitle = "White diamond = mean | Boxplot = median & quartiles | Violin = density",
    x = NULL, y = "Value"
  ) +
  theme(
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5),
    strip.text      = element_text(face = "bold")
  )

ggsave(
  "symbolic_regression_binary_distribution.png",
  plot = p,
  width = 11, height = 5, dpi = 300
)

cat("\nPlot saved to 'symbolic_regression_binary_distribution.png'\n")

# ------------------------------------------------------------
# 6. Summary statistics table
# ------------------------------------------------------------
summary_stats <- data.frame(
  Metric = metrics_to_plot,
  Mean   = sapply(results_table[, metrics_to_plot], mean,   na.rm = TRUE),
  Median = sapply(results_table[, metrics_to_plot], median, na.rm = TRUE),
  SD     = sapply(results_table[, metrics_to_plot], sd,     na.rm = TRUE),
  Q1     = sapply(results_table[, metrics_to_plot],
                  quantile, probs = 0.25, na.rm = TRUE),
  Q3     = sapply(results_table[, metrics_to_plot],
                  quantile, probs = 0.75, na.rm = TRUE)
)
summary_stats$IQR <- summary_stats$Q3 - summary_stats$Q1

cat("\nSummary statistics (after cleaning):\n")
print(summary_stats, row.names = FALSE)

write.csv(
  summary_stats,
  "symbolic_regression_binary_distribution_summary.csv",
  row.names = FALSE
)

cat("\nSummary saved to 'symbolic_regression_binary_distribution_summary.csv'\n")

##
## DECISION TREE - BINARY PHENOTYPE
##


setwd("C:\\Users\\Salvatore\\Desktop\\Dottorato\\richiesta_dati_reali\\EGAD00000000003_decriptati\\effettivamente_utili")


# ====================================================================
# 1. LOAD LIBRARIES
# ====================================================================

library(data.table)
library(rpart)
library(rpart.plot)
library(pROC)
library(ggplot2)


# ====================================================================
# 2. LOAD SELECTED SNPs AND DATA
# ====================================================================

features_selected <- read.delim(
  "selected_SNPs.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)

# Convert data frame -> character vector
features_selected <- features_selected$V1

df_train_fs <- read.csv(
  "train_selected_SNPs_no_double_quotes.csv",
  stringsAsFactors = FALSE
)

df_test_fs <- read.csv(
  "test_selected_SNPs_no_double_quotes.csv",
  stringsAsFactors = FALSE
)


# ====================================================================
# 3. CHECK DATA
# ====================================================================

cat("\n============================================\n")
cat("DATA CHECK\n")
cat("============================================\n")

cat("Number of selected SNPs:", length(features_selected), "\n")

if (length(features_selected) == 0) {
  stop("No SNPs were selected.")
}

# Check that all selected SNPs are present
missing_train <- setdiff(
  features_selected,
  names(df_train_fs)
)

missing_test <- setdiff(
  features_selected,
  names(df_test_fs)
)

if (length(missing_train) > 0) {
  stop(
    "The following SNPs are missing from training data:\n",
    paste(missing_train, collapse = ", ")
  )
}

if (length(missing_test) > 0) {
  stop(
    "The following SNPs are missing from testing data:\n",
    paste(missing_test, collapse = ", ")
  )
}


# Check phenotype
if (!"y" %in% names(df_train_fs)) {
  stop("Column 'y' not found in training data.")
}

if (!"y" %in% names(df_test_fs)) {
  stop("Column 'y' not found in testing data.")
}

if (!all(df_train_fs$y %in% c(0, 1))) {
  stop("Training phenotype y must contain only 0 and 1.")
}

if (!all(df_test_fs$y %in% c(0, 1))) {
  stop("Testing phenotype y must contain only 0 and 1.")
}


cat("\nTraining phenotype:\n")
print(table(df_train_fs$y))

cat("\nTesting phenotype:\n")
print(table(df_test_fs$y))


# ====================================================================
# 4. CONVERT PHENOTYPE TO FACTOR
# ====================================================================
#
# For rpart classification:
#
# 0 = Control
# 1 = Case
#
# We explicitly define the levels so that the probability
# extracted from predict(..., type = "prob") corresponds to class 1.
# ====================================================================

df_train_fs$y_factor <- factor(
  df_train_fs$y,
  levels = c(0, 1),
  labels = c("Control", "Case")
)

df_test_fs$y_factor <- factor(
  df_test_fs$y,
  levels = c(0, 1),
  labels = c("Control", "Case")
)


# ====================================================================
# 5. KEEP ONLY VARIABLES USED BY THE TREE
# ====================================================================

tree_train <- df_train_fs[
  ,
  c("y_factor", features_selected)
]

tree_test <- df_test_fs[
  ,
  c("y_factor", features_selected)
]


# ====================================================================
# 6. DECISION TREE
# ====================================================================

cat("\n============================================\n")
cat("DECISION TREE ANALYSIS\n")
cat("============================================\n")


# --------------------------------------------------------------------
# Grid of tree parameters
# --------------------------------------------------------------------

cp_values <- c(
  0.00001,
  0.00005,
  0.0001,
  0.0005,
  0.001,
  0.005,
  0.01,
  0.05,
  0.1
)

maxdepth_values <- c(
  2,
  3,
  4,
  5,
  6,
  8
)

minsplit_values <- c(
  10,
  20,
  30,
  50,
  100
)

minbucket_values <- c(
  5,
  10,
  20,
  30,
  50
)


# ====================================================================
# 7. METRIC FUNCTIONS
# ====================================================================

# --------------------------------------------------------------------
# Brier score
# --------------------------------------------------------------------

brier_score <- function(obs, pred) {
  
  mean(
    (obs - pred)^2,
    na.rm = TRUE
  )
  
}


# --------------------------------------------------------------------
# Log loss
# --------------------------------------------------------------------

log_loss <- function(obs, pred) {
  
  eps <- 1e-15
  
  pred <- pmin(
    pmax(pred, eps),
    1 - eps
  )
  
  -mean(
    obs * log(pred) +
      (1 - obs) * log(1 - pred),
    na.rm = TRUE
  )
  
}


# ====================================================================
# 8. INTERNAL CROSS-VALIDATION FOR PARAMETER SELECTION
# ====================================================================
#
# IMPORTANT:
#
# The TEST set is NOT used here.
#
# We use 5-fold CV inside the TRAINING set to select
# the tree parameters.
# ====================================================================

set.seed(123)

folds <- sample(
  rep(
    1:5,
    length.out = nrow(tree_train)
  )
)

results_grid <- data.frame()


cat("\nStarting 5-fold cross-validation...\n")


for (cp_val in cp_values) {
  
  for (maxd in maxdepth_values) {
    
    for (mins in minsplit_values) {
      
      for (minb in minbucket_values) {
        
        # Skip impossible combinations
        if (minb >= mins) {
          next
        }
        
        fold_auc <- c()
        
        for (k in 1:5) {
          
          train_cv <- tree_train[folds != k, ]
          valid_cv <- tree_train[folds == k, ]
          
          # ----------------------------------------------------------
          # Fit tree
          # ----------------------------------------------------------
          
          model_cv <- tryCatch(
            
            rpart(
              y_factor ~ .,
              data = train_cv,
              method = "class",
              
              control = rpart.control(
                cp = cp_val,
                maxdepth = maxd,
                minsplit = mins,
                minbucket = minb,
                xval = 0
              )
            ),
            
            error = function(e) {
              NULL
            }
            
          )
          
          
          if (is.null(model_cv)) {
            next
          }
          
          
          # ----------------------------------------------------------
          # Predict probabilities
          # ----------------------------------------------------------
          
          pred_prob <- tryCatch(
            
            predict(
              model_cv,
              newdata = valid_cv,
              type = "prob"
            )[,"Case"],
            
            error = function(e) {
              NULL
            }
            
          )
          
          
          if (is.null(pred_prob)) {
            next
          }
          
          
          # ----------------------------------------------------------
          # AUC
          # ----------------------------------------------------------
          
          roc_cv <- tryCatch(
            
            roc(
              response = valid_cv$y_factor,
              predictor = pred_prob,
              levels = c("Control", "Case"),
              direction = "<",
              quiet = TRUE
            ),
            
            error = function(e) {
              NULL
            }
            
          )
          
          
          if (!is.null(roc_cv)) {
            
            fold_auc <- c(
              fold_auc,
              as.numeric(auc(roc_cv))
            )
            
          }
          
        }
        
        
        # ------------------------------------------------------------
        # Mean CV AUC
        # ------------------------------------------------------------
        
        if (length(fold_auc) > 0) {
          
          mean_cv_auc <- mean(
            fold_auc,
            na.rm = TRUE
          )
          
        } else {
          
          mean_cv_auc <- NA
          
        }
        
        
        # ------------------------------------------------------------
        # Store results
        # ------------------------------------------------------------
        
        # ------------------------------------------------------------
        # Fit the tree on the FULL training set for these parameters
        # and compute train/test metrics so we can build distributions
        # afterwards. NOTE: test metrics are recorded ONLY for
        # diagnostics; model selection still uses CV_AUC only.
        # ------------------------------------------------------------
        
        model_full <- tryCatch(
          
          rpart(
            y_factor ~ .,
            data = tree_train,
            method = "class",
            control = rpart.control(
              cp = cp_val,
              maxdepth = maxd,
              minsplit = mins,
              minbucket = minb,
              xval = 0
            )
          ),
          
          error = function(e) NULL
          
        )
        
        if (!is.null(model_full)) {
          
          p_tr <- tryCatch(
            predict(model_full, newdata = tree_train, type = "prob")[, "Case"],
            error = function(e) NULL
          )
          
          p_te <- tryCatch(
            predict(model_full, newdata = tree_test, type = "prob")[, "Case"],
            error = function(e) NULL
          )
          
          if (!is.null(p_tr) && !is.null(p_te)) {
            
            roc_tr <- tryCatch(
              roc(response = tree_train$y_factor, predictor = p_tr,
                  levels = c("Control", "Case"),
                  direction = "<", quiet = TRUE),
              error = function(e) NULL
            )
            
            roc_te <- tryCatch(
              roc(response = tree_test$y_factor, predictor = p_te,
                  levels = c("Control", "Case"),
                  direction = "<", quiet = TRUE),
              error = function(e) NULL
            )
            
            auc_tr     <- if (!is.null(roc_tr)) as.numeric(auc(roc_tr)) else NA
            auc_te     <- if (!is.null(roc_te)) as.numeric(auc(roc_te)) else NA
            brier_tr   <- brier_score(df_train_fs$y, p_tr)
            brier_te   <- brier_score(df_test_fs$y,  p_te)
            logloss_tr <- log_loss(df_train_fs$y, p_tr)
            logloss_te <- log_loss(df_test_fs$y,  p_te)
            
          } else {
            auc_tr <- auc_te <- brier_tr <- brier_te <-
              logloss_tr <- logloss_te <- NA
          }
          
        } else {
          auc_tr <- auc_te <- brier_tr <- brier_te <-
            logloss_tr <- logloss_te <- NA
        }
        
        
        # ------------------------------------------------------------
        # Store results
        # ------------------------------------------------------------
        
        results_grid <- rbind(
          
          results_grid,
          
          data.frame(
            cp = cp_val,
            maxdepth = maxd,
            minsplit = mins,
            minbucket = minb,
            CV_AUC = mean_cv_auc,
            Train_AUC        = auc_tr,       # <-- new
            Test_AUC         = auc_te,       # <-- new
            Train_Brier      = brier_tr,     # <-- new
            Test_Brier       = brier_te,     # <-- new
            Train_LogLoss    = logloss_tr,   # <-- new
            Test_LogLoss     = logloss_te    # <-- new
          )
          
        )
        
      }
    }
  }
}


# ====================================================================
# 9. SELECT BEST PARAMETERS
# ====================================================================

results_grid <- results_grid[
  !is.na(results_grid$CV_AUC),
]

if (nrow(results_grid) == 0) {
  
  stop(
    "No valid decision tree model was obtained."
  )
  
}


results_grid <- results_grid[
  order(-results_grid$CV_AUC),
]


best_params <- results_grid[1, ]


cat("\n============================================\n")
cat("BEST TREE PARAMETERS\n")
cat("============================================\n")

print(best_params)

# ====================================================================
# 9b. DISTRIBUTION OF METRICS ACROSS ALL GRID-SEARCH RUNS
# ====================================================================

cat("\n============================================\n")
cat("DISTRIBUTION OF METRICS ACROSS RUNS\n")
cat("============================================\n")

# ---- 1. Long-format data for plotting ------------------------------

metrics_to_plot <- c(
  "Train_AUC", "Test_AUC",
  "Train_Brier", "Test_Brier",
  "Train_LogLoss", "Test_LogLoss"
)

results_long <- data.frame(
  Metric = factor(
    rep(metrics_to_plot, each = nrow(results_grid)),
    levels = metrics_to_plot,
    labels = c("Train AUC", "Test AUC",
               "Train Brier", "Test Brier",
               "Train Log Loss", "Test Log Loss")
  ),
  Value = unlist(results_grid[, metrics_to_plot], use.names = FALSE)
)

# ---- 2. Violin + boxplot + mean diamond ----------------------------

p <- ggplot(results_long, aes(x = Metric, y = Value, fill = Metric)) +
  geom_violin(alpha = 0.55, trim = FALSE, colour = NA) +
  geom_boxplot(width = 0.18, alpha = 0.85, outlier.size = 0.8,
               colour = "grey20") +
  stat_summary(fun = mean, geom = "point", shape = 23,
               size = 2.6, fill = "white", colour = "black") +
  facet_wrap(~ Metric, scales = "free", ncol = 3) +
  theme_bw() +
  labs(
    title    = "Distribution of Metrics Across All Grid-Search Runs",
    subtitle = "White diamond = mean | Boxplot = median & quartiles | Violin = density",
    x = NULL, y = "Value"
  ) +
  theme(
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5),
    strip.text      = element_text(face = "bold")
  )

ggsave(
  "decision_tree_distribution.png",
  plot = p,
  width = 12, height = 8, dpi = 300
)

cat("Distribution plot saved to 'decision_tree_distribution.png'\n")

# ---- 3. Summary statistics -----------------------------------------

summary_stats <- data.frame(
  Metric = metrics_to_plot,
  Mean   = sapply(results_grid[, metrics_to_plot], mean,   na.rm = TRUE),
  Median = sapply(results_grid[, metrics_to_plot], median, na.rm = TRUE),
  SD     = sapply(results_grid[, metrics_to_plot], sd,     na.rm = TRUE),
  Q1     = sapply(results_grid[, metrics_to_plot],
                  quantile, probs = 0.25, na.rm = TRUE),
  Q3     = sapply(results_grid[, metrics_to_plot],
                  quantile, probs = 0.75, na.rm = TRUE)
)

summary_stats$IQR <- summary_stats$Q3 - summary_stats$Q1

print(summary_stats, row.names = FALSE)

write.csv(
  summary_stats,
  "decision_tree_distribution_summary.csv",
  row.names = FALSE
)

cat("Distribution summary saved to 'decision_tree_distribution_summary.csv'\n")

# ====================================================================
# 10. FIT FINAL TREE ON ALL TRAINING DATA
# ====================================================================

best_model <- rpart(
  
  y_factor ~ .,
  
  data = tree_train,
  
  method = "class",
  
  control = rpart.control(
    
    cp = best_params$cp,
    maxdepth = best_params$maxdepth,
    minsplit = best_params$minsplit,
    minbucket = best_params$minbucket,
    
    xval = 10
    
  )
  
)


cat("\nFinal decision tree fitted successfully.\n")


# ====================================================================
# 11. PREDICT TRAINING AND TESTING DATA
# ====================================================================

# Probability of CASE (y = 1)

p_train <- predict(
  best_model,
  newdata = tree_train,
  type = "prob"
)[,"Case"]


p_test <- predict(
  best_model,
  newdata = tree_test,
  type = "prob"
)[,"Case"]


# Class prediction

class_train <- predict(
  best_model,
  newdata = tree_train,
  type = "class"
)

class_test <- predict(
  best_model,
  newdata = tree_test,
  type = "class"
)


# ====================================================================
# 12. AUC-ROC
# ====================================================================

roc_train <- roc(
  response = tree_train$y_factor,
  predictor = p_train,
  levels = c("Control", "Case"),
  direction = "<",
  quiet = TRUE
)

roc_test <- roc(
  response = tree_test$y_factor,
  predictor = p_test,
  levels = c("Control", "Case"),
  direction = "<",
  quiet = TRUE
)


auc_train <- as.numeric(
  auc(roc_train)
)

auc_test <- as.numeric(
  auc(roc_test)
)


# ====================================================================
# 13. BRIER SCORE
# ====================================================================

brier_train <- brier_score(
  df_train_fs$y,
  p_train
)

brier_test <- brier_score(
  df_test_fs$y,
  p_test
)


# ====================================================================
# 14. LOG LOSS
# ====================================================================

logloss_train <- log_loss(
  df_train_fs$y,
  p_train
)

logloss_test <- log_loss(
  df_test_fs$y,
  p_test
)


# ====================================================================
# 15. CLASSIFICATION ACCURACY
# ====================================================================
#
# This is NOT one of your three requested metrics, but is useful
# for describing how the tree classifies samples.
# ====================================================================

pred_train_binary <- ifelse(
  class_train == "Case",
  1,
  0
)

pred_test_binary <- ifelse(
  class_test == "Case",
  1,
  0
)


accuracy_train <- mean(
  pred_train_binary == df_train_fs$y
)

accuracy_test <- mean(
  pred_test_binary == df_test_fs$y
)


# ====================================================================
# 16. CONFUSION MATRICES
# ====================================================================

confusion_train <- table(
  Actual = df_train_fs$y,
  Predicted = pred_train_binary
)

confusion_test <- table(
  Actual = df_test_fs$y,
  Predicted = pred_test_binary
)


# ====================================================================
# 17. FINAL PERFORMANCE
# ====================================================================

performance <- data.frame(
  
  Dataset = c(
    "Training",
    "Testing"
  ),
  
  AUC_ROC = c(
    auc_train,
    auc_test
  ),
  
  Brier_Score = c(
    brier_train,
    brier_test
  ),
  
  Log_Loss = c(
    logloss_train,
    logloss_test
  ),
  
  Accuracy = c(
    accuracy_train,
    accuracy_test
  )
  
)


cat("\n============================================\n")
cat("DECISION TREE PERFORMANCE\n")
cat("============================================\n")

print(
  performance,
  row.names = FALSE
)


# ====================================================================
# 18. CONFUSION MATRICES
# ====================================================================

cat("\n============================================\n")
cat("TRAINING CONFUSION MATRIX\n")
cat("============================================\n")

print(confusion_train)


cat("\n============================================\n")
cat("TESTING CONFUSION MATRIX\n")
cat("============================================\n")

print(confusion_test)


# ====================================================================
# 19. DISPLAY TREE
# ====================================================================

cat("\n============================================\n")
cat("DECISION TREE\n")
cat("============================================\n")

print(best_model)


# Plot tree

png(
  "decision_tree_plot.png",
  width = 2400,
  height = 1800,
  res = 250
)

rpart.plot(
  best_model,
  
  type = 4,
  extra = 104,
  
  fallen.leaves = TRUE,
  
  box.palette = "GnBu",
  
  shadow.col = "gray",
  
  nn = TRUE,
  
  main = "Decision Tree - Binary Phenotype"
)

dev.off()


# ====================================================================
# 20. EXTRACT TREE RULES
# ====================================================================

cat("\n============================================\n")
cat("DECISION TREE RULES\n")
cat("============================================\n")

tree_rules <- rpart.rules(
  best_model
)

print(tree_rules)


# Save rules

writeLines(
  capture.output(
    rpart.rules(best_model)
  ),
  "decision_tree_rules.txt"
)


# ====================================================================
# 21. SAVE FULL TREE SUMMARY
# ====================================================================

sink(
  "decision_tree_full_summary.txt"
)

cat("============================================\n")
cat("DECISION TREE - FULL SUMMARY\n")
cat("============================================\n\n")


cat("Number of selected SNPs:\n")
cat(length(features_selected), "\n\n")


cat("Selected SNPs:\n")
cat(paste(features_selected, collapse = "\n"))
cat("\n\n")


cat("Best parameters:\n")
print(best_params)


cat("\n\nTree:\n")
print(best_model)


cat("\n\nTree rules:\n")
print(
  rpart.rules(best_model)
)


cat("\n\nVariable importance:\n")

if (!is.null(best_model$variable.importance)) {
  
  importance_df <- data.frame(
    
    Variable =
      names(best_model$variable.importance),
    
    Importance =
      as.numeric(best_model$variable.importance)
    
  )
  
  importance_df <-
    importance_df[
      order(
        -importance_df$Importance
      ),
    ]
  
  print(
    importance_df,
    row.names = FALSE
  )
  
} else {
  
  cat(
    "No variable importance available.\n"
  )
  
}


sink()


# ====================================================================
# 22. VARIABLE IMPORTANCE
# ====================================================================

if (!is.null(best_model$variable.importance)) {
  
  importance_df <- data.frame(
    
    Variable =
      names(best_model$variable.importance),
    
    Importance =
      as.numeric(best_model$variable.importance)
    
  )
  
  importance_df <-
    importance_df[
      order(
        -importance_df$Importance
      ),
    ]
  
  write.csv(
    importance_df,
    "decision_tree_variable_importance.csv",
    row.names = FALSE
  )
  
}


# ====================================================================
# 23. SAVE PERFORMANCE METRICS
# ====================================================================

write.csv(
  performance,
  "decision_tree_metrics.csv",
  row.names = FALSE
)


write.table(
  performance,
  "decision_tree_metrics.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ====================================================================
# 24. SAVE GRID SEARCH RESULTS
# ====================================================================

write.csv(
  results_grid,
  "decision_tree_CV_grid_search.csv",
  row.names = FALSE
)


# ====================================================================
# 25. SAVE PREDICTIONS
# ====================================================================

predictions_train <- data.frame(
  
  Actual = df_train_fs$y,
  
  Predicted_Class =
    pred_train_binary,
  
  Predicted_Probability =
    p_train,
  
  Set = "Training"
  
)


predictions_test <- data.frame(
  
  Actual = df_test_fs$y,
  
  Predicted_Class =
    pred_test_binary,
  
  Predicted_Probability =
    p_test,
  
  Set = "Testing"
  
)


predictions_all <- rbind(
  predictions_train,
  predictions_test
)


write.csv(
  predictions_all,
  "decision_tree_predictions.csv",
  row.names = FALSE
)


# ====================================================================
# 26. SAVE MODEL
# ====================================================================

saveRDS(
  best_model,
  "best_decision_tree_model.rds"
)


# ====================================================================
# 27. PRINT FINAL SUMMARY
# ====================================================================

cat("\n============================================\n")
cat("FINAL SUMMARY\n")
cat("============================================\n")

cat(
  "\nNumber of selected SNPs:",
  length(features_selected),
  "\n"
)

cat(
  "\nBest CV AUC:",
  round(best_params$CV_AUC, 6),
  "\n"
)

cat(
  "\nFinal tree parameters:\n"
)

cat(
  "  cp = ",
  best_params$cp,
  "\n",
  sep = ""
)

cat(
  "  maxdepth = ",
  best_params$maxdepth,
  "\n",
  sep = ""
)

cat(
  "  minsplit = ",
  best_params$minsplit,
  "\n",
  sep = ""
)

cat(
  "  minbucket = ",
  best_params$minbucket,
  "\n",
  sep = ""
)


cat("\nPerformance:\n")
print(
  performance,
  row.names = FALSE
)


cat("\n\nTree rules:\n")
print(
  rpart.rules(best_model)
)


cat("\n============================================\n")
cat("FILES SAVED\n")
cat("============================================\n")

cat("1. decision_tree_metrics.csv\n")
cat("2. decision_tree_metrics.txt\n")
cat("3. decision_tree_rules.txt\n")
cat("4. decision_tree_full_summary.txt\n")
cat("5. decision_tree_variable_importance.csv\n")
cat("6. decision_tree_CV_grid_search.csv\n")
cat("7. decision_tree_predictions.csv\n")
cat("8. decision_tree_plot.png\n")
cat("9. best_decision_tree_model.rds\n")


cat("\nDecision tree analysis completed successfully.\n")


##
## PARTE RELATIVA A SVM BINARIO 
##


set.seed(123)


############################################################
# 2. LIBRARIES
############################################################

library(data.table)
library(foreach)
library(doParallel)
library(parallel)
library(caret)
library(kernlab)
library(pROC)


############################################################
# 3. PARALLEL COMPUTING
############################################################

num_cores <- detectCores() - 1

if (is.na(num_cores) || num_cores < 1) {
  num_cores <- 1
}

cat("\nNumber of CPU cores:", num_cores, "\n")

cl <- makeCluster(num_cores)
registerDoParallel(cl)


setwd(
  "C:\\Users\\Salvatore\\Desktop\\Dottorato\\richiesta_dati_reali\\EGAD00000000003_decriptati\\effettivamente_utili"
)


############################################################
# 5. INPUT FILES
############################################################

selected_file <- "selected_SNPs.txt"

train_file <- "train_selected_SNPs.csv"

test_file <- "test_selected_SNPs.csv"


############################################################
# 6. LOAD SELECTED SNP LIST
############################################################

feature_selected <- fread(
  selected_file,
  data.table = FALSE,
  header = FALSE
)


############################################################
# 7. CHECK SELECTED SNP COLUMN
############################################################

if (!"V1" %in% colnames(feature_selected)) {
  
  stop(
    "The file selected_SNPs.txt must contain a column called 'V1'."
  )
  
}


selected_snps <- as.character(
  feature_selected$V1
)


selected_snps <- selected_snps[
  !is.na(selected_snps) &
    selected_snps != ""
]


selected_snps <- unique(
  selected_snps
)


cat(
  "\nNumber of selected SNPs:",
  length(selected_snps),
  "\n"
)


############################################################
# 8. LOAD TRAINING AND TESTING DATA
############################################################

train_df <- fread(
  train_file,
  data.table = FALSE
)

test_df <- fread(
  test_file,
  data.table = FALSE
)

names(train_df) <- gsub(
  '"',
  '',
  names(train_df)
)

names(test_df) <- gsub(
  '"',
  '',
  names(test_df)
)


############################################################
# 9. CHECK PHENOTYPE COLUMN
############################################################

if (!"phenotype" %in% colnames(train_df)) {
  
  stop(
    "Column 'phenotype' not found in training dataset."
  )
  
}


if (!"phenotype" %in% colnames(test_df)) {
  
  stop(
    "Column 'phenotype' not found in testing dataset."
  )
  
}


############################################################
# 10. CHECK SELECTED SNPs
############################################################

missing_train <- setdiff(
  selected_snps,
  colnames(train_df)
)

missing_test <- setdiff(
  selected_snps,
  colnames(test_df)
)


if (length(missing_train) > 0) {
  
  cat(
    "\nThe following SNPs are missing from the training dataset:\n"
  )
  
  print(
    missing_train
  )
  
  stop(
    "Some selected SNPs are missing from train_selected_SNPs.csv."
  )
  
}


if (length(missing_test) > 0) {
  
  cat(
    "\nThe following SNPs are missing from the testing dataset:\n"
  )
  
  print(
    missing_test
  )
  
  stop(
    "Some selected SNPs are missing from test_selected_SNPs.csv."
  )
  
}


############################################################
# 11. SELECT PHENOTYPE + SNPs
############################################################

train_data <- train_df[
  ,
  c(
    "phenotype",
    selected_snps
  ),
  drop = FALSE
]

test_data <- test_df[
  ,
  c(
    "phenotype",
    selected_snps
  ),
  drop = FALSE
]


############################################################
# 12. CONVERT SNPs TO NUMERIC
############################################################

convert_snp <- function(x) {
  
  if (is.factor(x)) {
    
    return(
      as.numeric(
        as.character(x)
      )
    )
    
  }
  
  return(
    as.numeric(x)
  )
  
}


train_data[selected_snps] <- lapply(
  train_data[selected_snps],
  convert_snp
)

test_data[selected_snps] <- lapply(
  test_data[selected_snps],
  convert_snp
)


############################################################
# 13. CHECK FOR NA VALUES IN SNPs
############################################################

if (
  anyNA(
    train_data[
      ,
      selected_snps,
      drop = FALSE
    ]
  )
) {
  
  cat(
    "\nNA values detected in training SNPs.\n"
  )
  
  print(
    colSums(
      is.na(
        train_data[
          ,
          selected_snps,
          drop = FALSE
        ]
      )
    )
  )
  
  stop(
    "Training SNP matrix contains NA values."
  )
  
}


if (
  anyNA(
    test_data[
      ,
      selected_snps,
      drop = FALSE
    ]
  )
) {
  
  cat(
    "\nNA values detected in testing SNPs.\n"
  )
  
  print(
    colSums(
      is.na(
        test_data[
          ,
          selected_snps,
          drop = FALSE
        ]
      )
    )
  )
  
  stop(
    "Testing SNP matrix contains NA values."
  )
  
}


############################################################
# 14. CHECK PHENOTYPE
############################################################

if (anyNA(train_data$phenotype)) {
  
  stop(
    "NA values detected in training phenotype."
  )
  
}


if (anyNA(test_data$phenotype)) {
  
  stop(
    "NA values detected in testing phenotype."
  )
  
}


############################################################
# 15. CHECK BINARY PHENOTYPE
############################################################

train_unique <- sort(
  unique(
    train_data$phenotype
  )
)

test_unique <- sort(
  unique(
    test_data$phenotype
  )
)


cat(
  "\nTraining phenotype values:\n"
)

print(
  train_unique
)


cat(
  "\nTesting phenotype values:\n"
)

print(
  test_unique
)


if (length(train_unique) != 2) {
  
  stop(
    "Training phenotype must contain exactly two classes."
  )
  
}


if (length(test_unique) != 2) {
  
  stop(
    "Testing phenotype must contain exactly two classes."
  )
  
}


############################################################
# 16. FORCE PHENOTYPE TO BINARY FACTOR
#
# Assumes the phenotype is coded 0/1.
# Class1 = phenotype 1
# Class0 = phenotype 0
############################################################

if (
  !all(
    train_unique %in% c(
      0,
      1
    )
  )
) {
  
  stop(
    "The phenotype must be coded as 0 and 1."
  )
  
}


if (
  !all(
    test_unique %in% c(
      0,
      1
    )
  )
) {
  
  stop(
    "The phenotype must be coded as 0 and 1."
  )
  
}


train_data$phenotype <- factor(
  ifelse(
    train_data$phenotype == 1,
    "Class1",
    "Class0"
  ),
  levels = c(
    "Class1",
    "Class0"
  )
)


test_data$phenotype <- factor(
  ifelse(
    test_data$phenotype == 1,
    "Class1",
    "Class0"
  ),
  levels = c(
    "Class1",
    "Class0"
  )
)


############################################################
# 17. CHECK CLASS DISTRIBUTION
############################################################

cat(
  "\nTraining class distribution:\n"
)

print(
  table(
    train_data$phenotype
  )
)


cat(
  "\nTesting class distribution:\n"
)

print(
  table(
    test_data$phenotype
  )
)


############################################################
# 18. SAMPLE NUMBERS
############################################################

cat(
  "\nTraining samples:",
  nrow(train_data),
  "\n"
)

cat(
  "Testing samples :",
  nrow(test_data),
  "\n"
)

cat(
  "Number of SNPs  :",
  length(selected_snps),
  "\n"
)


############################################################
# 19. CROSS-VALIDATION SETTINGS
############################################################

ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  allowParallel = TRUE
)


############################################################
# 20. SVM HYPERPARAMETER GRID
############################################################

svm_grid <- expand.grid(
  cost = c(
    0.001,
    0.01,
    0.1,
    1
  )
)


############################################################
# 21. TRAIN LINEAR SVM
############################################################

cat(
  "\n========================================\n"
)

cat(
  "TRAINING BINARY LINEAR SVM\n"
)

cat(
  "========================================\n\n"
)


svm_model <- train(
  phenotype ~ .,
  data = train_data,
  method = "svmLinear2",
  trControl = ctrl,
  tuneGrid = svm_grid,
  metric = "ROC",
  preProcess = c(
    "center",
    "scale"
  )
)


############################################################
# 22. DISPLAY BEST MODEL
############################################################

cat(
  "\n========================================\n"
)

cat(
  "BEST MODEL\n"
)

cat(
  "========================================\n\n"
)


print(
  svm_model$bestTune
)


############################################################
# 23. EXTRACT FINAL e1071 MODEL
############################################################

final_model <- svm_model$finalModel


cat(
  "\n========================================\n"
)

cat(
  "FINAL MODEL CLASS\n"
)

cat(
  "========================================\n\n"
)

print(
  class(final_model)
)


############################################################
# 24. DISPLAY MODEL STRUCTURE
############################################################

cat(
  "\n========================================\n"
)

cat(
  "FINAL MODEL STRUCTURE\n"
)

cat(
  "========================================\n\n"
)

print(
  names(final_model)
)


############################################################
# 25. EXTRACT PREPROCESSING PARAMETERS
############################################################

preprocess_mean <- svm_model$preProcess$mean

preprocess_sd <- svm_model$preProcess$std


############################################################
# 26. CHECK PREPROCESSING PARAMETERS
############################################################

cat(
  "\n========================================\n"
)

cat(
  "CENTERING VALUES\n"
)

cat(
  "========================================\n\n"
)

print(
  preprocess_mean
)


cat(
  "\n========================================\n"
)

cat(
  "SCALING VALUES\n"
)

cat(
  "========================================\n\n"
)

print(
  preprocess_sd
)


############################################################
# 27. CHECK PREPROCESSING PARAMETERS
############################################################

if (
  !identical(
    names(preprocess_mean),
    names(preprocess_sd)
  )
) {
  
  stop(
    "The names/order of preprocessing means and SDs do not match."
  )
  
}


############################################################
# 28. REORDER PARAMETERS
############################################################

preprocess_mean <- preprocess_mean[
  selected_snps
]

preprocess_sd <- preprocess_sd[
  selected_snps
]


############################################################
# 29. CHECK ZERO SD
############################################################

if (
  any(
    preprocess_sd == 0
  )
) {
  
  cat(
    "\nVariables with zero standard deviation:\n"
  )
  
  print(
    names(preprocess_sd)[
      preprocess_sd == 0
    ]
  )
  
  stop(
    "At least one SNP has zero variance."
  )
  
}


############################################################
# 30. PREDICT TRAINING PROBABILITIES
############################################################

pred_train_prob <- predict(
  svm_model,
  newdata = train_data,
  type = "prob"
)


############################################################
# 31. PREDICT TESTING PROBABILITIES
############################################################

pred_test_prob <- predict(
  svm_model,
  newdata = test_data,
  type = "prob"
)


############################################################
# 32. EXTRACT PROBABILITY OF POSITIVE CLASS
############################################################

prob_train <- as.numeric(
  pred_train_prob$Class1
)

prob_test <- as.numeric(
  pred_test_prob$Class1
)


############################################################
# 33. CLASS PREDICTIONS
############################################################

pred_train_class <- predict(
  svm_model,
  newdata = train_data,
  type = "raw"
)

pred_test_class <- predict(
  svm_model,
  newdata = test_data,
  type = "raw"
)


############################################################
# 34. ORIGINAL PHENOTYPE AS 0/1
############################################################

y_train <- ifelse(
  train_data$phenotype == "Class1",
  1,
  0
)

y_test <- ifelse(
  test_data$phenotype == "Class1",
  1,
  0
)


############################################################
# 35. SAFE PROBABILITIES
############################################################

eps <- 1e-15


prob_train_safe <- pmin(
  pmax(
    prob_train,
    eps
  ),
  1 - eps
)


prob_test_safe <- pmin(
  pmax(
    prob_test,
    eps
  ),
  1 - eps
)


############################################################
# 36. AUC ROC
############################################################

train_roc <- roc(
  response = y_train,
  predictor = prob_train,
  levels = c(
    0,
    1
  ),
  direction = "<",
  quiet = TRUE
)


test_roc <- roc(
  response = y_test,
  predictor = prob_test,
  levels = c(
    0,
    1
  ),
  direction = "<",
  quiet = TRUE
)


train_AUC_ROC <- as.numeric(
  auc(
    train_roc
  )
)

test_AUC_ROC <- as.numeric(
  auc(
    test_roc
  )
)


############################################################
# 37. BRIER SCORE
############################################################

train_Brier_Score <- mean(
  (
    prob_train -
      y_train
  )^2
)


test_Brier_Score <- mean(
  (
    prob_test -
      y_test
  )^2
)


############################################################
# 38. LOG LOSS
############################################################

train_Log_Loss <- -mean(
  y_train *
    log(
      prob_train_safe
    ) +
    (
      1 - y_train
    ) *
    log(
      1 - prob_train_safe
    )
)


test_Log_Loss <- -mean(
  y_test *
    log(
      prob_test_safe
    ) +
    (
      1 - y_test
    ) *
    log(
      1 - prob_test_safe
    )
)


############################################################
# 39. CLASSIFICATION ACCURACY
############################################################

train_accuracy <- mean(
  pred_train_class ==
    train_data$phenotype
)

test_accuracy <- mean(
  pred_test_class ==
    test_data$phenotype
)


############################################################
# 40. PRINT PERFORMANCE
############################################################

cat(
  "\n========================================\n"
)

cat(
  "MODEL PERFORMANCE\n"
)

cat(
  "========================================\n\n"
)


cat(
  "Training AUC ROC   :",
  train_AUC_ROC,
  "\n"
)

cat(
  "Testing AUC ROC    :",
  test_AUC_ROC,
  "\n\n"
)


cat(
  "Training Brier     :",
  train_Brier_Score,
  "\n"
)

cat(
  "Testing Brier      :",
  test_Brier_Score,
  "\n\n"
)


cat(
  "Training Log Loss  :",
  train_Log_Loss,
  "\n"
)

cat(
  "Testing Log Loss   :",
  test_Log_Loss,
  "\n\n"
)


cat(
  "Training Accuracy  :",
  train_accuracy,
  "\n"
)

cat(
  "Testing Accuracy   :",
  test_accuracy,
  "\n"
)


############################################################
# 41. CREATE SNP MATRICES
############################################################

X_train <- train_data[
  ,
  selected_snps,
  drop = FALSE
]

X_test <- test_data[
  ,
  selected_snps,
  drop = FALSE
]


X_train_matrix <- as.matrix(
  X_train
)

storage.mode(
  X_train_matrix
) <- "numeric"


X_test_matrix <- as.matrix(
  X_test
)

storage.mode(
  X_test_matrix
) <- "numeric"


############################################################
# 42. STANDARDIZE TRAINING MATRIX
############################################################

X_train_standardized <- sweep(
  X_train_matrix,
  2,
  as.numeric(
    preprocess_mean
  ),
  FUN = "-"
)

X_train_standardized <- sweep(
  X_train_standardized,
  2,
  as.numeric(
    preprocess_sd
  ),
  FUN = "/"
)


############################################################
# 43. STANDARDIZE TEST MATRIX
############################################################

X_test_standardized <- sweep(
  X_test_matrix,
  2,
  as.numeric(
    preprocess_mean
  ),
  FUN = "-"
)

X_test_standardized <- sweep(
  X_test_standardized,
  2,
  as.numeric(
    preprocess_sd
  ),
  FUN = "/"
)


############################################################
# 44. EXTRACT LINEAR SVM SUPPORT VECTORS
############################################################

SV <- final_model$SV


############################################################
# 45. EXTRACT LINEAR SVM COEFFICIENTS
############################################################

coefs <- final_model$coefs


############################################################
# 46. EXTRACT SVM RHO
############################################################

rho <- as.numeric(
  final_model$rho
)


############################################################
# 47. EXTRACT STANDARDIZED SVM WEIGHTS
#
# For the linear e1071 SVM:
#
# f(x) = sum(w_j * x_j) - rho
#
# where:
#
# w = sum_i alpha_i * y_i * SV_i
############################################################

w_standardized <- as.numeric(
  crossprod(
    coefs[, 1],
    SV
  )
)

names(
  w_standardized
) <- selected_snps


############################################################
# 48. DISPLAY STANDARDIZED COEFFICIENTS
############################################################

cat(
  "\n========================================\n"
)

cat(
  "STANDARDIZED SVM DECISION-FUNCTION COEFFICIENTS\n"
)

cat(
  "========================================\n\n"
)

print(
  w_standardized
)


############################################################
# 49. ORIGINAL-SCALE COEFFICIENTS
############################################################

beta_original <-
  w_standardized /
  as.numeric(
    preprocess_sd
  )

names(
  beta_original
) <- selected_snps


############################################################
# 50. ORIGINAL-SCALE DECISION INTERCEPT
############################################################

intercept_original <-
  -rho -
  sum(
    w_standardized *
      as.numeric(
        preprocess_mean
      ) /
      as.numeric(
        preprocess_sd
      )
  )


############################################################
# 51. CREATE COEFFICIENT TABLE
############################################################

coefficient_table <- data.frame(
  
  SNP = selected_snps,
  
  Mean = as.numeric(
    preprocess_mean
  ),
  
  SD = as.numeric(
    preprocess_sd
  ),
  
  Coefficient_standardized =
    as.numeric(
      w_standardized
    ),
  
  Coefficient_original =
    as.numeric(
      beta_original
    ),
  
  stringsAsFactors = FALSE,
  
  check.names = FALSE
  
)


############################################################
# 52. ABSOLUTE COEFFICIENTS
############################################################

coefficient_table$Absolute_coefficient_standardized <-
  abs(
    coefficient_table$
      Coefficient_standardized
  )


coefficient_table$Absolute_coefficient_original <-
  abs(
    coefficient_table$
      Coefficient_original
  )


############################################################
# 53. SORT BY ABSOLUTE EFFECT
############################################################

coefficient_table_sorted <-
  coefficient_table[
    order(
      -coefficient_table$
        Absolute_coefficient_standardized
    ),
    ,
    drop = FALSE
  ]


############################################################
# 54. DISPLAY SNP COEFFICIENTS
############################################################

cat(
  "\n========================================\n"
)

cat(
  "SNP COEFFICIENTS\n"
)

cat(
  "========================================\n\n"
)

print(
  coefficient_table_sorted,
  row.names = FALSE
)


############################################################
# 55. BUILD STANDARDIZED DECISION-FUNCTION FORMULA
#
# f(z) = -rho + sum(w_j * z_j)
#
# z_j = (SNP_j - mean_j) / SD_j
############################################################

standardized_formula <- paste0(
  "Decision_score = ",
  format(
    -rho,
    digits = 15
  )
)


for (j in seq_along(w_standardized)) {
  
  coefficient <- w_standardized[j]
  
  sign_symbol <- ifelse(
    coefficient >= 0,
    " + ",
    " - "
  )
  
  standardized_formula <- paste0(
    standardized_formula,
    sign_symbol,
    format(
      abs(coefficient),
      digits = 15
    ),
    " * (( ",
    selected_snps[j],
    " - ",
    format(
      preprocess_mean[j],
      digits = 15
    ),
    " ) / ",
    format(
      preprocess_sd[j],
      digits = 15
    ),
    ")"
  )
  
}


############################################################
# 56. BUILD ORIGINAL-SCALE DECISION-FUNCTION FORMULA
############################################################

original_formula <- paste0(
  "Decision_score = ",
  format(
    intercept_original,
    digits = 15
  )
)


for (j in seq_along(beta_original)) {
  
  coefficient <- beta_original[j]
  
  sign_symbol <- ifelse(
    coefficient >= 0,
    " + ",
    " - "
  )
  
  original_formula <- paste0(
    original_formula,
    sign_symbol,
    format(
      abs(coefficient),
      digits = 15
    ),
    " * ",
    selected_snps[j]
  )
  
}


############################################################
# 57. PRINT FORMULAS
############################################################

cat(
  "\n========================================\n"
)

cat(
  "LINEAR SVM DECISION FUNCTION\n"
)

cat(
  "========================================\n\n"
)


cat(
  "STANDARDIZED FORM:\n\n"
)

cat(
  standardized_formula,
  "\n\n"
)


cat(
  "ORIGINAL SNP SCALE:\n\n"
)

cat(
  original_formula,
  "\n\n"
)


cat(
  "Standardized decision intercept:",
  -rho,
  "\n"
)


cat(
  "Original-scale intercept:",
  intercept_original,
  "\n"
)


############################################################
# 58. FORMULA-BASED TRAINING DECISION SCORES
############################################################

pred_train_decision_formula <-
  as.numeric(
    X_train_standardized %*%
      w_standardized
  ) -
  rho


############################################################
# 59. FORMULA-BASED TESTING DECISION SCORES
############################################################

pred_test_decision_formula <-
  as.numeric(
    X_test_standardized %*%
      w_standardized
  ) -
  rho


############################################################
# 60. DIRECT e1071 SVM DECISION SCORES
#
# e1071::predict.svm does not use type = "decision".
# Decision values are obtained with decision.values = TRUE.
############################################################

pred_train_svm <- predict(
  final_model,
  newdata = X_train_standardized,
  decision.values = TRUE
)

pred_test_svm <- predict(
  final_model,
  newdata = X_test_standardized,
  decision.values = TRUE
)


pred_train_decision <- as.numeric(
  attr(
    pred_train_svm,
    "decision.values"
  )
)

pred_test_decision <- as.numeric(
  attr(
    pred_test_svm,
    "decision.values"
  )
)


############################################################
# 61. FORMULA RECONSTRUCTION CHECK
############################################################

train_difference <-
  pred_train_decision -
  pred_train_decision_formula


test_difference <-
  pred_test_decision -
  pred_test_decision_formula


max_train_difference <-
  max(
    abs(
      train_difference
    )
  )


max_test_difference <-
  max(
    abs(
      test_difference
    )
  )


mean_train_difference <-
  mean(
    abs(
      train_difference
    )
  )


mean_test_difference <-
  mean(
    abs(
      test_difference
    )
  )


tolerance <- 1e-8


formula_matches_train <-
  max_train_difference <
  tolerance


formula_matches_test <-
  max_test_difference <
  tolerance


############################################################
# 62. PRINT RECONSTRUCTION CHECK
############################################################

cat(
  "\n========================================\n"
)

cat(
  "FORMULA RECONSTRUCTION CHECK\n"
)

cat(
  "========================================\n\n"
)


cat(
  "Maximum absolute difference - TRAIN:",
  max_train_difference,
  "\n"
)


cat(
  "Mean absolute difference - TRAIN:",
  mean_train_difference,
  "\n\n"
)


cat(
  "Maximum absolute difference - TEST:",
  max_test_difference,
  "\n"
)


cat(
  "Mean absolute difference - TEST:",
  mean_test_difference,
  "\n\n"
)


cat(
  "Formula reproduces TRAIN decision scores:",
  formula_matches_train,
  "\n"
)


cat(
  "Formula reproduces TEST decision scores:",
  formula_matches_test,
  "\n"
)


############################################################
# 63. CALCULATE TRAINING SNP CONTRIBUTIONS
############################################################

train_contributions <- sweep(
  X_train_standardized,
  2,
  w_standardized,
  FUN = "*"
)


train_contributions <- as.data.frame(
  train_contributions,
  check.names = FALSE
)


colnames(
  train_contributions
) <- selected_snps


############################################################
# 64. CALCULATE TESTING SNP CONTRIBUTIONS
############################################################

test_contributions <- sweep(
  X_test_standardized,
  2,
  w_standardized,
  FUN = "*"
)


test_contributions <- as.data.frame(
  test_contributions,
  check.names = FALSE
)


colnames(
  test_contributions
) <- selected_snps


############################################################
# 65. CHECK CONTRIBUTION SUMS
############################################################

train_contribution_sum <-
  -rho +
  rowSums(
    train_contributions
  )


test_contribution_sum <-
  -rho +
  rowSums(
    test_contributions
  )


cat(
  "\n========================================\n"
)

cat(
  "CONTRIBUTION CHECK\n"
)

cat(
  "========================================\n\n"
)


cat(
  "Maximum train difference:",
  max(
    abs(
      pred_train_decision -
        train_contribution_sum
    )
  ),
  "\n"
)


cat(
  "Maximum test difference:",
  max(
    abs(
      pred_test_decision -
        test_contribution_sum
    )
  ),
  "\n"
)


############################################################
# 66. TRAIN CONTRIBUTION TABLE
############################################################

train_contribution_table <- data.frame(
  
  phenotype_observed = y_train,
  
  Intercept_decision = -rho,
  
  train_contributions,
  
  Decision_score = pred_train_decision,
  
  Probability_Class1 = prob_train,
  
  Predicted_class =
    as.character(
      pred_train_class
    ),
  
  stringsAsFactors = FALSE,
  
  check.names = FALSE
  
)


############################################################
# 67. TEST CONTRIBUTION TABLE
############################################################

test_contribution_table <- data.frame(
  
  phenotype_observed = y_test,
  
  Intercept_decision = -rho,
  
  test_contributions,
  
  Decision_score = pred_test_decision,
  
  Probability_Class1 = prob_test,
  
  Predicted_class =
    as.character(
      pred_test_class
    ),
  
  stringsAsFactors = FALSE,
  
  check.names = FALSE
  
)


############################################################
# 68. ADD PREDICTION ERROR
############################################################

train_contribution_table$Error_probability <-
  y_train -
  prob_train


train_contribution_table$Absolute_error_probability <-
  abs(
    train_contribution_table$
      Error_probability
  )


test_contribution_table$Error_probability <-
  y_test -
  prob_test


test_contribution_table$Absolute_error_probability <-
  abs(
    test_contribution_table$
      Error_probability
  )


############################################################
# 69. STANDARDIZED GENOTYPE TABLES
############################################################

train_standardized_table <- data.frame(
  X_train_standardized,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


test_standardized_table <- data.frame(
  X_test_standardized,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


############################################################
# 70. ORIGINAL GENOTYPE TABLES
############################################################

train_genotype_table <- data.frame(
  X_train,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


test_genotype_table <- data.frame(
  X_test,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


############################################################
# 71. TRAINING PREDICTION TABLE
############################################################

train_prediction_table <- data.frame(
  
  Observed_phenotype = y_train,
  
  Predicted_probability_Class1 =
    prob_train,
  
  Predicted_class =
    as.character(
      pred_train_class
    ),
  
  Decision_score =
    pred_train_decision,
  
  Error_probability =
    y_train -
    prob_train,
  
  Absolute_error_probability =
    abs(
      y_train -
        prob_train
    ),
  
  stringsAsFactors = FALSE
  
)


############################################################
# 72. TESTING PREDICTION TABLE
############################################################

test_prediction_table <- data.frame(
  
  Observed_phenotype = y_test,
  
  Predicted_probability_Class1 =
    prob_test,
  
  Predicted_class =
    as.character(
      pred_test_class
    ),
  
  Decision_score =
    pred_test_decision,
  
  Error_probability =
    y_test -
    prob_test,
  
  Absolute_error_probability =
    abs(
      y_test -
        prob_test
    ),
  
  stringsAsFactors = FALSE
  
)


############################################################
# 73. PERFORMANCE TABLE
############################################################

performance_table <- data.frame(
  
  Dataset = c(
    "Training",
    "Testing"
  ),
  
  AUC_ROC = c(
    train_AUC_ROC,
    test_AUC_ROC
  ),
  
  Brier_Score = c(
    train_Brier_Score,
    test_Brier_Score
  ),
  
  Log_Loss = c(
    train_Log_Loss,
    test_Log_Loss
  ),
  
  Accuracy = c(
    train_accuracy,
    test_accuracy
  ),
  
  N = c(
    nrow(train_data),
    nrow(test_data)
  ),
  
  stringsAsFactors = FALSE
  
)


############################################################
# 74. MODEL INFORMATION TABLE
############################################################

model_information <- data.frame(
  
  Parameter = c(
    
    "Model",
    
    "Kernel",
    
    "Cost",
    
    "Number_of_SNPs",
    
    "Training_samples",
    
    "Testing_samples",
    
    "Positive_class",
    
    "Negative_class",
    
    "Cross_validation",
    
    "CV_folds",
    
    "CV_repeats",
    
    "Standardized_decision_intercept",
    
    "Original_scale_decision_intercept",
    
    "Number_of_support_vectors",
    
    "Formula_matches_training",
    
    "Formula_matches_testing",
    
    "Max_training_difference",
    
    "Max_testing_difference",
    
    "Training_AUC_ROC",
    
    "Testing_AUC_ROC",
    
    "Training_Brier_Score",
    
    "Testing_Brier_Score",
    
    "Training_Log_Loss",
    
    "Testing_Log_Loss"
    
  ),
  
  Value = c(
    
    "svmLinear2",
    
    "linear",
    
    as.character(
      svm_model$bestTune$cost
    ),
    
    as.character(
      length(selected_snps)
    ),
    
    as.character(
      nrow(train_data)
    ),
    
    as.character(
      nrow(test_data)
    ),
    
    "Class1 / phenotype = 1",
    
    "Class0 / phenotype = 0",
    
    "repeatedcv",
    
    "10",
    
    "5",
    
    as.character(
      -rho
    ),
    
    as.character(
      intercept_original
    ),
    
    as.character(
      nrow(final_model$SV)
    ),
    
    as.character(
      formula_matches_train
    ),
    
    as.character(
      formula_matches_test
    ),
    
    as.character(
      max_train_difference
    ),
    
    as.character(
      max_test_difference
    ),
    
    as.character(
      train_AUC_ROC
    ),
    
    as.character(
      test_AUC_ROC
    ),
    
    as.character(
      train_Brier_Score
    ),
    
    as.character(
      test_Brier_Score
    ),
    
    as.character(
      train_Log_Loss
    ),
    
    as.character(
      test_Log_Loss
    )
    
  ),
  
  stringsAsFactors = FALSE
  
)


############################################################
# 75. PREPROCESSING PARAMETERS
############################################################

preprocessing_table <- data.frame(
  
  SNP = selected_snps,
  
  Mean =
    as.numeric(
      preprocess_mean
    ),
  
  SD =
    as.numeric(
      preprocess_sd
    ),
  
  stringsAsFactors = FALSE,
  
  check.names = FALSE
  
)


############################################################
# 76. SAVE COEFFICIENT TABLE
############################################################

write.csv(
  coefficient_table,
  "SVM_binary_SNP_coefficients.csv",
  row.names = FALSE
)


write.csv(
  coefficient_table_sorted,
  "SVM_binary_SNP_coefficients_sorted_by_effect.csv",
  row.names = FALSE
)


############################################################
# 77. SAVE PREPROCESSING PARAMETERS
############################################################

write.csv(
  preprocessing_table,
  "SVM_binary_preprocessing_parameters.csv",
  row.names = FALSE
)


############################################################
# 78. SAVE STANDARDIZED GENOTYPES
############################################################

write.csv(
  train_standardized_table,
  "SVM_binary_train_standardized_SNPs.csv",
  row.names = FALSE
)


write.csv(
  test_standardized_table,
  "SVM_binary_test_standardized_SNPs.csv",
  row.names = FALSE
)


############################################################
# 79. SAVE ORIGINAL GENOTYPES
############################################################

write.csv(
  train_genotype_table,
  "SVM_binary_train_original_SNPs.csv",
  row.names = FALSE
)


write.csv(
  test_genotype_table,
  "SVM_binary_test_original_SNPs.csv",
  row.names = FALSE
)


############################################################
# 80. SAVE INDIVIDUAL SNP CONTRIBUTIONS
############################################################

write.csv(
  train_contribution_table,
  "SVM_binary_train_SNP_contributions.csv",
  row.names = FALSE
)


write.csv(
  test_contribution_table,
  "SVM_binary_test_SNP_contributions.csv",
  row.names = FALSE
)


############################################################
# 81. SAVE PREDICTIONS
############################################################

write.csv(
  train_prediction_table,
  "SVM_binary_train_predictions.csv",
  row.names = FALSE
)


write.csv(
  test_prediction_table,
  "SVM_binary_test_predictions.csv",
  row.names = FALSE
)


############################################################
# 82. SAVE PROBABILITIES
############################################################

train_probability_table <- data.frame(
  
  Observed_phenotype = y_train,
  
  Probability_Class0 =
    pred_train_prob$Class0,
  
  Probability_Class1 =
    pred_train_prob$Class1,
  
  Predicted_class =
    as.character(
      pred_train_class
    ),
  
  Decision_score =
    pred_train_decision,
  
  stringsAsFactors = FALSE
  
)


test_probability_table <- data.frame(
  
  Observed_phenotype = y_test,
  
  Probability_Class0 =
    pred_test_prob$Class0,
  
  Probability_Class1 =
    pred_test_prob$Class1,
  
  Predicted_class =
    as.character(
      pred_test_class
    ),
  
  Decision_score =
    pred_test_decision,
  
  stringsAsFactors = FALSE
  
)


write.csv(
  train_probability_table,
  "SVM_binary_train_probabilities.csv",
  row.names = FALSE
)


write.csv(
  test_probability_table,
  "SVM_binary_test_probabilities.csv",
  row.names = FALSE
)


############################################################
# 83. SAVE PERFORMANCE
############################################################

write.csv(
  performance_table,
  "SVM_binary_performance.csv",
  row.names = FALSE
)


############################################################
# 84. SAVE MODEL INFORMATION
############################################################

write.csv(
  model_information,
  "SVM_binary_model_information.csv",
  row.names = FALSE
)


############################################################
# 85. SAVE DECISION-FUNCTION FORMULAS
############################################################

formula_file <- file(
  "SVM_binary_formula.txt",
  open = "wt"
)


writeLines(
  "============================================================",
  formula_file
)

writeLines(
  "BINARY LINEAR SVM DECISION FUNCTION",
  formula_file
)

writeLines(
  "============================================================",
  formula_file
)

writeLines(
  "",
  formula_file
)

writeLines(
  "IMPORTANT:",
  formula_file
)

writeLines(
  "The formulas below describe the LINEAR SVM DECISION SCORE.",
  formula_file
)

writeLines(
  "They are not direct formulas for probability.",
  formula_file
)

writeLines(
  "",
  formula_file
)

writeLines(
  "STANDARDIZED FORM:",
  formula_file
)

writeLines(
  standardized_formula,
  formula_file
)

writeLines(
  "",
  formula_file
)

writeLines(
  "ORIGINAL SNP SCALE:",
  formula_file
)

writeLines(
  original_formula,
  formula_file
)

writeLines(
  "",
  formula_file
)

writeLines(
  paste0(
    "Standardized decision intercept = ",
    -rho
  ),
  formula_file
)

writeLines(
  paste0(
    "Original-scale intercept = ",
    intercept_original
  ),
  formula_file
)

writeLines(
  "",
  formula_file
)

writeLines(
  "============================================================",
  formula_file
)

writeLines(
  "FORMULA RECONSTRUCTION CHECK",
  formula_file
)

writeLines(
  "============================================================",
  formula_file
)

writeLines(
  paste0(
    "Maximum absolute difference, training = ",
    max_train_difference
  ),
  formula_file
)

writeLines(
  paste0(
    "Maximum absolute difference, testing = ",
    max_test_difference
  ),
  formula_file
)

writeLines(
  paste0(
    "Formula reproduces training decision scores = ",
    formula_matches_train
  ),
  formula_file
)

writeLines(
  paste0(
    "Formula reproduces testing decision scores = ",
    formula_matches_test
  ),
  formula_file
)


close(
  formula_file
)


############################################################
# 86. SAVE DETAILED COEFFICIENT FILE
############################################################

coefficient_file <- file(
  "SVM_binary_detailed_SNP_coefficients.txt",
  open = "wt"
)


writeLines(
  "============================================================",
  coefficient_file
)

writeLines(
  "DETAILED BINARY LINEAR SVM SNP COEFFICIENTS",
  coefficient_file
)

writeLines(
  "============================================================",
  coefficient_file
)

writeLines(
  "",
  coefficient_file
)

writeLines(
  "Coefficients refer to the linear SVM decision function.",
  coefficient_file
)

writeLines(
  "",
  coefficient_file
)

writeLines(
  paste0(
    "Standardized decision intercept = ",
    -rho
  ),
  coefficient_file
)

writeLines(
  paste0(
    "Original-scale decision intercept = ",
    intercept_original
  ),
  coefficient_file
)

writeLines(
  "",
  coefficient_file
)


for (j in seq_along(w_standardized)) {
  
  snp_name <- selected_snps[j]
  
  writeLines(
    paste0(
      "SNP: ",
      snp_name
    ),
    coefficient_file
  )
  
  writeLines(
    paste0(
      "  Mean = ",
      preprocess_mean[j]
    ),
    coefficient_file
  )
  
  writeLines(
    paste0(
      "  SD = ",
      preprocess_sd[j]
    ),
    coefficient_file
  )
  
  writeLines(
    paste0(
      "  Standardized decision coefficient = ",
      w_standardized[j]
    ),
    coefficient_file
  )
  
  writeLines(
    paste0(
      "  Original-scale decision coefficient = ",
      beta_original[j]
    ),
    coefficient_file
  )
  
  writeLines(
    "",
    coefficient_file
  )
  
}


close(
  coefficient_file
)


############################################################
# 87. SAVE TRAINED CARET MODEL
############################################################

saveRDS(
  svm_model,
  "SVM_binary_linear_model.rds"
)


############################################################
# 88. SAVE FINAL e1071 MODEL
############################################################

saveRDS(
  final_model,
  "SVM_binary_linear_final_e1071_model.rds"
)


############################################################
# 89. COMPLETE ANALYSIS OBJECT
############################################################

analysis_results <- list(
  
  model = svm_model,
  
  final_model = final_model,
  
  best_tune =
    svm_model$bestTune,
  
  selected_snps =
    selected_snps,
  
  preprocessing_mean =
    preprocess_mean,
  
  preprocessing_sd =
    preprocess_sd,
  
  standardized_coefficients =
    w_standardized,
  
  original_coefficients =
    beta_original,
  
  standardized_decision_intercept =
    -rho,
  
  original_decision_intercept =
    intercept_original,
  
  rho =
    rho,
  
  support_vectors =
    SV,
  
  standardized_formula =
    standardized_formula,
  
  original_formula =
    original_formula,
  
  train_probability =
    prob_train,
  
  test_probability =
    prob_test,
  
  train_decision_score =
    pred_train_decision,
  
  test_decision_score =
    pred_test_decision,
  
  train_predictions =
    pred_train_class,
  
  test_predictions =
    pred_test_class,
  
  train_AUC_ROC =
    train_AUC_ROC,
  
  test_AUC_ROC =
    test_AUC_ROC,
  
  train_Brier_Score =
    train_Brier_Score,
  
  test_Brier_Score =
    test_Brier_Score,
  
  train_Log_Loss =
    train_Log_Loss,
  
  test_Log_Loss =
    test_Log_Loss,
  
  train_accuracy =
    train_accuracy,
  
  test_accuracy =
    test_accuracy,
  
  coefficient_table =
    coefficient_table,
  
  train_contributions =
    train_contribution_table,
  
  test_contributions =
    test_contribution_table,
  
  formula_matches_train =
    formula_matches_train,
  
  formula_matches_test =
    formula_matches_test,
  
  max_train_difference =
    max_train_difference,
  
  max_test_difference =
    max_test_difference
  
)


saveRDS(
  analysis_results,
  "SVM_binary_complete_analysis_results.rds"
)


############################################################
# 90. PRINT TOP SNPs
############################################################

cat(
  "\n========================================\n"
)

cat(
  "TOP SNPs BY ABSOLUTE STANDARDIZED EFFECT\n"
)

cat(
  "========================================\n\n"
)


print(
  head(
    coefficient_table_sorted,
    21
  ),
  row.names = FALSE
)


############################################################
# 91. FINAL SUMMARY
############################################################

cat(
  "\n\n"
)

cat(
  "============================================================\n"
)

cat(
  "BINARY SVM ANALYSIS COMPLETE\n"
)

cat(
  "============================================================\n\n"
)


cat(
  "Number of training samples:",
  nrow(train_data),
  "\n"
)


cat(
  "Number of testing samples:",
  nrow(test_data),
  "\n"
)


cat(
  "Number of SNPs:",
  length(selected_snps),
  "\n\n"
)


cat(
  "Best SVM cost:",
  svm_model$bestTune$cost,
  "\n"
)


cat(
  "Number of support vectors:",
  nrow(final_model$SV),
  "\n\n"
)


cat(
  "Training AUC ROC:",
  train_AUC_ROC,
  "\n"
)


cat(
  "Testing AUC ROC:",
  test_AUC_ROC,
  "\n\n"
)


cat(
  "Training Brier Score:",
  train_Brier_Score,
  "\n"
)


cat(
  "Testing Brier Score:",
  test_Brier_Score,
  "\n\n"
)


cat(
  "Training Log Loss:",
  train_Log_Loss,
  "\n"
)


cat(
  "Testing Log Loss:",
  test_Log_Loss,
  "\n\n"
)


cat(
  "Training Accuracy:",
  train_accuracy,
  "\n"
)


cat(
  "Testing Accuracy:",
  test_accuracy,
  "\n\n"
)


cat(
  "Standardized decision intercept:",
  -rho,
  "\n"
)


cat(
  "Original-scale decision intercept:",
  intercept_original,
  "\n\n"
)


cat(
  "Formula reconstruction max error - TRAIN:",
  max_train_difference,
  "\n"
)


cat(
  "Formula reconstruction max error - TEST:",
  max_test_difference,
  "\n\n"
)


cat(
  "Formula reproduces TRAIN decision scores:",
  formula_matches_train,
  "\n"
)


cat(
  "Formula reproduces TEST decision scores:",
  formula_matches_test,
  "\n\n"
)


cat(
  "Formula saved to: SVM_binary_formula.txt\n"
)

cat(
  "Coefficients saved to: SVM_binary_SNP_coefficients.csv\n"
)

cat(
  "Individual contributions saved to:",
  "SVM_binary_test_SNP_contributions.csv\n"
)

cat(
  "Performance saved to: SVM_binary_performance.csv\n"
)

cat(
  "Complete model saved to:",
  "SVM_binary_complete_analysis_results.rds\n"
)

cat(
  "Prediction table saved to:",
  "SVM_binary_test_predictions.csv\n"
)


############################################################
# 92. STOP PARALLEL CLUSTER
############################################################

stopCluster(
  cl
)

registerDoSEQ()




#
A <- read.csv("symbolic_results_summary_binary_1_100.csv")
b <- read.csv("symbolic_results_summary_binary_101_500.csv")
c <- read.csv("symbolic_results_summary_binary_501_1000.csv")
unificati <- rbind(A,b,c)
write.table(unificati, file = "symbolic_results_summary_binary_1_1000.txt", sep = "\t", quote = F)
