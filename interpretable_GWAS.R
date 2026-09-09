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
best_test_r2 <- -Inf
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
              minsplit = mins, # how many samples in a node after each split to keep going
              minbucket = minb,
              cp = cp_val, # how much it has to improve after each split to keep going
              maxdepth = maxd,
              xval = 5
            )
          )
          
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
            train_r2 = train_r2,
            test_r2 = test_r2,
            train_mae = train_mae,
            test_mae = test_mae,
            train_rmse = train_rmse,
            test_rmse = test_rmse,
            tree_size = ifelse(is.null(dt_model$frame), 1, nrow(dt_model$frame))
          ))
          
          # Update best model (still based on R², but could be changed to MAE)
          if(!is.na(test_r2) && test_r2 > best_test_r2) {
            best_test_r2 <- test_r2
            best_params <- list(cp = cp_val, maxdepth = maxd, 
                                minsplit = mins, minbucket = minb)
            best_model <- dt_model
            
            cat(sprintf("New best: cp=%.6f, depth=%d, minsplit=%d, minbucket=%d\n",
                        cp_val, maxd, mins, minb))
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


# FINE PARTE FENOTIPO CONTINUO
# INIZIO PARTE FENOTIPO BINARIO
