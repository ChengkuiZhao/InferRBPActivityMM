# =============================================================================
# Survival Analysis and Random Forest Model Module
# =============================================================================

source("config/config.R")

# Load required libraries
suppressPackageStartupMessages({
  library(randomForestSRC)
  library(survival)
  library(survminer)
})

#' Perform hyperparameter tuning with cross-validation
#'
#' @param RBP_activity RBP activity matrix (samples x features)
#' @param survival_data Survival data with OS.time and OS columns
#' @param cv_folds Number of cross-validation folds
#' @param ntree_grid Vector of ntree values to test
#' @param mtry_grid Vector of mtry values to test
#' @param nodesize_grid Vector of nodesize values to test
#' @param seed Random seed for reproducibility
#' @return List with best parameters and error matrix
tune_hyperparameters <- function(RBP_activity, survival_data, 
                                cv_folds = 5,
                                ntree_grid = c(50, 100, 200, 500, 1000),
                                mtry_grid = c(10, 13, 15, 20, 30),
                                nodesize_grid = c(5, 10, 20, 50, 100),
                                seed = 123) {
  
  cat("Starting hyperparameter tuning...\n")
  
  # Prepare data
  data_train <- data.frame(t(RBP_activity), survival_data[, c("OS.time", "OS")])
  n_samples <- nrow(data_train)
  
  # Create cross-validation folds
  set.seed(seed)
  sample_indices <- sample(1:n_samples)
  fold_size <- floor(n_samples / cv_folds)
  
  cv_folds_list <- list()
  for (i in 1:cv_folds) {
    start_idx <- (i - 1) * fold_size + 1
    end_idx <- ifelse(i == cv_folds, n_samples, i * fold_size)
    cv_folds_list[[i]] <- sample_indices[start_idx:end_idx]
  }
  
  # Initialize error array
  n_ntree <- length(ntree_grid)
  n_mtry <- length(mtry_grid)
  n_nodesize <- length(nodesize_grid)
  
  error_array <- array(0, dim = c(n_ntree, n_mtry, n_nodesize, cv_folds))
  
  # Hyperparameter search
  total_combinations <- n_ntree * n_mtry * n_nodesize * cv_folds
  current_combination <- 0
  
  for (i in 1:n_ntree) {
    for (j in 1:n_mtry) {
      for (k in 1:n_nodesize) {
        for (fold in 1:cv_folds) {
          current_combination <- current_combination + 1
          
          # Split data
          test_indices <- cv_folds_list[[fold]]
          train_data <- data_train[-test_indices, ]
          test_data <- data_train[test_indices, ]
          
          # Train model
          tryCatch({
            model <- rfsrc(
              Surv(OS.time, OS) ~ ., 
              data = train_data,
              na.action = MODEL_CONFIG$missing_data_action,
              ntree = ntree_grid[i],
              mtry = mtry_grid[j],
              nodesize = nodesize_grid[k]
            )
            
            # Predict and calculate error
            prediction <- predict.rfsrc(
              model, 
              newdata = test_data,
              na.action = MODEL_CONFIG$missing_data_action
            )
            
            error_array[i, j, k, fold] <- prediction$err.rate[ntree_grid[i]]
            
          }, error = function(e) {
            warning("Error in combination ", current_combination, ": ", e$message)
            error_array[i, j, k, fold] <- Inf
          })
          
          if (current_combination %% 10 == 0) {
            cat("Completed", current_combination, "/", total_combinations, "combinations\n")
          }
        }
      }
    }
  }
  
  # Calculate mean errors
  mean_errors <- apply(error_array, c(1, 2, 3), mean)
  
  # Find best parameters
  best_indices <- which(mean_errors == min(mean_errors, na.rm = TRUE), arr.ind = TRUE)[1, ]
  
  best_params <- list(
    ntree = ntree_grid[best_indices[1]],
    mtry = mtry_grid[best_indices[2]],
    nodesize = nodesize_grid[best_indices[3]],
    error = mean_errors[best_indices[1], best_indices[2], best_indices[3]]
  )
  
  cat("Best parameters found:\n")
  cat("  ntree:", best_params$ntree, "\n")
  cat("  mtry:", best_params$mtry, "\n")  
  cat("  nodesize:", best_params$nodesize, "\n")
  cat("  CV error:", round(best_params$error, 4), "\n")
  
  return(list(
    best_params = best_params,
    error_array = error_array,
    mean_errors = mean_errors,
    cv_folds = cv_folds_list
  ))
}

#' Train final survival model with best parameters
#'
#' @param RBP_activity_train Training RBP activity matrix
#' @param survival_data_train Training survival data
#' @param best_params Best hyperparameters from tuning
#' @return Trained randomForestSRC model
train_survival_model <- function(RBP_activity_train, survival_data_train, best_params) {
  
  cat("Training final survival model...\n")
  
  # Prepare training data
  data_train <- data.frame(t(RBP_activity_train), survival_data_train[, c("OS.time", "OS")])
  
  # Train model with best parameters
  model <- rfsrc(
    Surv(OS.time, OS) ~ .,
    data = data_train,
    na.action = MODEL_CONFIG$missing_data_action,
    ntree = best_params$ntree,
    mtry = best_params$mtry,
    nodesize = best_params$nodesize
  )
  
  cat("Model training completed!\n")
  cat("  OOB error rate:", round(model$err.rate[model$ntree], 4), "\n")
  cat("  Variable importance available:", !is.null(model$importance), "\n")
  
  return(model)
}

#' Predict risk scores using trained model
#'
#' @param model Trained randomForestSRC model
#' @param RBP_activity_test Test RBP activity matrix
#' @return Vector of predicted risk scores
predict_risk_scores <- function(model, RBP_activity_test) {
  
  cat("Predicting risk scores...\n")
  
  # Prepare test data
  data_test <- data.frame(t(RBP_activity_test))
  
  # Make predictions
  predictions <- predict.rfsrc(
    model,
    newdata = data_test,
    na.action = MODEL_CONFIG$missing_data_action
  )
  
  risk_scores <- predictions$predicted
  names(risk_scores) <- colnames(RBP_activity_test)
  
  cat("Risk score prediction completed!\n")
  cat("  Number of predictions:", length(risk_scores), "\n")
  cat("  Risk score range:", paste(round(range(risk_scores), 4), collapse = " to "), "\n")
  
  return(risk_scores)
}

#' Classify samples into risk groups
#'
#' @param risk_scores Vector of risk scores
#' @param method Classification method ("median", "tertile", "custom")
#' @param threshold Custom threshold (if method = "custom")
#' @return Factor with risk group classifications
classify_risk_groups <- function(risk_scores, method = "median", threshold = NULL) {
  
  cat("Classifying samples into risk groups...\n")
  
  if (method == "median") {
    threshold <- median(risk_scores, na.rm = TRUE)
    risk_groups <- ifelse(risk_scores > threshold, "HighRisk", "LowRisk")
    
  } else if (method == "tertile") {
    tertiles <- quantile(risk_scores, c(1/3, 2/3), na.rm = TRUE)
    risk_groups <- cut(risk_scores, 
                      breaks = c(-Inf, tertiles[1], tertiles[2], Inf),
                      labels = c("LowRisk", "IntermediateRisk", "HighRisk"))
    
  } else if (method == "custom" && !is.null(threshold)) {
    risk_groups <- ifelse(risk_scores > threshold, "HighRisk", "LowRisk")
    
  } else {
    stop("Invalid classification method or missing threshold")
  }
  
  risk_groups <- factor(risk_groups)
  
  cat("Risk group classification completed!\n")
  cat("  Method:", method, "\n")
  if (method %in% c("median", "custom")) {
    cat("  Threshold:", round(threshold, 4), "\n")
  }
  cat("  Group distribution:\n")
  print(table(risk_groups))
  
  return(risk_groups)
}

#' Evaluate model performance
#'
#' @param risk_scores Vector of risk scores  
#' @param survival_data Survival data with OS.time and OS columns
#' @param risk_groups Risk group classifications
#' @return List with performance metrics
evaluate_model_performance <- function(risk_scores, survival_data, risk_groups) {
  
  cat("Evaluating model performance...\n")
  
  # Prepare data
  eval_data <- data.frame(
    risk_score = risk_scores,
    risk_group = risk_groups,
    OS.time = survival_data$OS.time,
    OS = survival_data$OS
  )
  
  # Survival analysis
  surv_object <- Surv(time = eval_data$OS.time, event = eval_data$OS)
  
  # Kaplan-Meier survival curves
  km_fit <- survfit(surv_object ~ risk_group, data = eval_data)
  
  # Log-rank test
  logrank_test <- survdiff(surv_object ~ risk_group, data = eval_data)
  logrank_pvalue <- pchisq(logrank_test$chisq, df = 1, lower.tail = FALSE)
  
  # Cox proportional hazards model
  cox_model <- coxph(surv_object ~ risk_score, data = eval_data)
  cox_summary <- summary(cox_model)
  
  # C-index (concordance index)
  c_index <- cox_summary$concordance[1]
  
  results <- list(
    km_fit = km_fit,
    logrank_pvalue = logrank_pvalue,
    cox_model = cox_model,
    cox_summary = cox_summary,
    c_index = c_index,
    eval_data = eval_data
  )
  
  cat("Model evaluation completed!\n")
  cat("  Log-rank p-value:", format.pval(logrank_pvalue), "\n")
  cat("  Cox HR (95% CI):", round(cox_summary$conf.int[1], 3), 
      "(", round(cox_summary$conf.int[3], 3), "-", round(cox_summary$conf.int[4], 3), ")\n")
  cat("  C-index:", round(c_index, 3), "\n")
  
  return(results)
}

#' Complete survival analysis pipeline
#'
#' @param RBP_activity_train Training RBP activity matrix
#' @param RBP_activity_validation Validation RBP activity matrix  
#' @param survival_data_train Training survival data
#' @param survival_data_validation Validation survival data
#' @param tune_hyperparams Whether to perform hyperparameter tuning
#' @param custom_params Custom parameters (if not tuning)
#' @return List with complete analysis results
run_survival_analysis <- function(RBP_activity_train, RBP_activity_validation,
                                 survival_data_train, survival_data_validation,
                                 tune_hyperparams = TRUE, custom_params = NULL) {
  
  cat("Starting complete survival analysis pipeline...\n")
  
  # Hyperparameter tuning (if requested)
  if (tune_hyperparams) {
    tuning_results <- tune_hyperparameters(
      RBP_activity = RBP_activity_train,
      survival_data = survival_data_train,
      cv_folds = MODEL_CONFIG$cv_folds,
      ntree_grid = MODEL_CONFIG$ntree_grid,
      mtry_grid = MODEL_CONFIG$mtry_grid,
      nodesize_grid = MODEL_CONFIG$nodesize_grid
    )
    best_params <- tuning_results$best_params
  } else {
    best_params <- custom_params %||% MODEL_CONFIG$rf_params
    tuning_results <- NULL
  }
  
  # Train final model
  model <- train_survival_model(RBP_activity_train, survival_data_train, best_params)
  
  # Predict on validation set
  risk_scores <- predict_risk_scores(model, RBP_activity_validation)
  risk_groups <- classify_risk_groups(risk_scores, method = "median")
  
  # Evaluate performance
  performance <- evaluate_model_performance(risk_scores, survival_data_validation, risk_groups)
  
  # Return complete results
  list(
    model = model,
    best_params = best_params,
    tuning_results = tuning_results,
    risk_scores = risk_scores,
    risk_groups = risk_groups,
    performance = performance
  )
}