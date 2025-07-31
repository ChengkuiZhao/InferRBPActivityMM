# =============================================================================
# Cofactors Analysis Module
# =============================================================================

source("config/config.R")

# Load required libraries
suppressPackageStartupMessages({
  library(survival)
  library(survminer)
})

#' Prepare cofactors data for Cox regression
#'
#' @param phenotype_data Phenotype data frame
#' @param risk_scores Vector of risk scores
#' @param survival_data Survival data
#' @return Cleaned data frame for Cox analysis
prepare_cofactors_data <- function(phenotype_data, risk_scores, survival_data) {
  
  cat("Preparing cofactors data for analysis...\n")
  
  # Extract relevant columns
  cofactor_cols <- c('demographic.gender', 'demographic.race', 
                    'demographic.ethnicity', 'diagnoses.age_at_diagnosis')
  
  # Create base data frame
  cox_data <- data.frame(
    phenotype_data[names(risk_scores), cofactor_cols, drop = FALSE]
  )
  
  # Add risk scores
  cox_data$risk_score <- risk_scores
  
  # Add survival data
  cox_data$OS.time <- survival_data[rownames(cox_data), "OS.time"]
  cox_data$OS <- survival_data[rownames(cox_data), "OS"]
  
  # Clean and encode categorical variables
  cox_data <- clean_cofactors(cox_data)
  
  cat("Cofactors data prepared with", nrow(cox_data), "samples\n")
  
  return(cox_data)
}

#' Clean and encode cofactors variables
#'
#' @param cox_data Raw Cox data frame
#' @return Cleaned data frame with encoded variables
clean_cofactors <- function(cox_data) {
  
  cat("Cleaning and encoding cofactors...\n")
  
  # Gender encoding
  cox_data$gender <- ifelse(cox_data$demographic.gender == "female", 0, 1)
  
  # Race encoding
  cox_data$race <- NA
  cox_data$race[cox_data$demographic.race == "white"] <- 0
  cox_data$race[cox_data$demographic.race == "black or african american"] <- 1
  cox_data$race[cox_data$demographic.race == "asian"] <- 2
  cox_data$race[cox_data$demographic.race == "other"] <- 3
  # Set problematic entries to NA
  cox_data$race[cox_data$demographic.race %in% c("no reported", "not allowed to collect", "not reported")] <- NA
  
  # Ethnicity encoding
  cox_data$ethnicity <- NA
  cox_data$ethnicity[cox_data$demographic.ethnicity == "not hispanic or latino"] <- 1
  cox_data$ethnicity[cox_data$demographic.ethnicity == "hispanic or latino"] <- 2
  cox_data$ethnicity[cox_data$demographic.ethnicity == "not reported"] <- NA
  
  # Age conversion (from days to years)
  cox_data$age_years <- cox_data$diagnoses.age_at_diagnosis / 365.25
  
  # Normalize risk score (0-1 scale)
  cox_data$risk_score_norm <- (cox_data$risk_score - min(cox_data$risk_score, na.rm = TRUE)) / 
                             (max(cox_data$risk_score, na.rm = TRUE) - min(cox_data$risk_score, na.rm = TRUE))
  
  # Select final columns
  final_data <- cox_data[, c("gender", "race", "ethnicity", "age_years", 
                            "risk_score", "risk_score_norm", "OS.time", "OS")]
  
  # Convert to appropriate types
  final_data$gender <- as.numeric(final_data$gender)
  final_data$race <- as.numeric(final_data$race)
  final_data$ethnicity <- as.numeric(final_data$ethnicity)
  
  cat("Cofactors cleaning completed\n")
  
  return(final_data)
}

#' Perform multivariate Cox regression analysis
#'
#' @param cox_data Prepared Cox data frame
#' @param include_normalized Whether to include normalized risk score
#' @return List with Cox model results
perform_cox_analysis <- function(cox_data, include_normalized = FALSE) {
  
  cat("Performing Cox regression analysis...\n")
  
  # Remove rows with missing survival data
  cox_data_clean <- cox_data[!is.na(cox_data$OS.time) & !is.na(cox_data$OS), ]
  
  # Model 1: Original risk score
  risk_var <- ifelse(include_normalized, "risk_score_norm", "risk_score")
  
  formula1 <- as.formula(paste("Surv(OS.time, OS) ~", risk_var))
  cox_model1 <- coxph(formula1, data = cox_data_clean)
  
  # Model 2: Full multivariate model
  formula2 <- as.formula(paste("Surv(OS.time, OS) ~ gender + race + ethnicity + age_years +", risk_var))
  cox_model2 <- coxph(formula2, data = cox_data_clean)
  
  # Model 3: Simplified model (excluding variables with many missing values)
  formula3 <- as.formula(paste("Surv(OS.time, OS) ~ gender + age_years +", risk_var))
  cox_model3 <- coxph(formula3, data = cox_data_clean)
  
  results <- list(
    data = cox_data_clean,
    model_univariate = cox_model1,
    model_full = cox_model2,
    model_simplified = cox_model3,
    summary_univariate = summary(cox_model1),
    summary_full = summary(cox_model2),
    summary_simplified = summary(cox_model3)
  )
  
  cat("Cox regression analysis completed\n")
  cat("  Univariate model - Risk score HR:", round(exp(coef(cox_model1)), 3), "\n")
  cat("  Full model - Risk score HR:", round(exp(coef(cox_model2)[risk_var]), 3), "\n")
  cat("  Simplified model - Risk score HR:", round(exp(coef(cox_model3)[risk_var]), 3), "\n")
  
  return(results)
}

#' Create Cox regression summary table
#'
#' @param cox_results Results from perform_cox_analysis
#' @return Data frame with formatted results
create_cox_summary_table <- function(cox_results) {
  
  cat("Creating Cox regression summary table...\n")
  
  # Extract results from full model
  full_summary <- cox_results$summary_full
  
  # Create summary table
  summary_table <- data.frame(
    Variable = rownames(full_summary$coefficients),
    HR = round(exp(full_summary$coefficients[, "coef"]), 3),
    CI_lower = round(full_summary$conf.int[, "lower .95"], 3),
    CI_upper = round(full_summary$conf.int[, "upper .95"], 3),
    P_value = format.pval(full_summary$coefficients[, "Pr(>|z|)"], digits = 3),
    stringsAsFactors = FALSE
  )
  
  # Format confidence intervals
  summary_table$CI_95 <- paste0("(", summary_table$CI_lower, "-", summary_table$CI_upper, ")")
  
  # Add significance stars
  p_vals <- as.numeric(summary_table$P_value)
  summary_table$Significance <- ""
  summary_table$Significance[p_vals < 0.001] <- "***"
  summary_table$Significance[p_vals < 0.01 & p_vals >= 0.001] <- "**"
  summary_table$Significance[p_vals < 0.05 & p_vals >= 0.01] <- "*"
  summary_table$Significance[p_vals < 0.1 & p_vals >= 0.05] <- "."
  
  # Clean up variable names
  summary_table$Variable <- gsub("_", " ", summary_table$Variable)
  summary_table$Variable <- tools::toTitleCase(summary_table$Variable)
  
  cat("Cox regression summary table created\n")
  
  return(summary_table[, c("Variable", "HR", "CI_95", "P_value", "Significance")])
}

#' Assess proportional hazards assumption
#'
#' @param cox_model Cox proportional hazards model
#' @return List with test results and plots
assess_proportional_hazards <- function(cox_model) {
  
  cat("Assessing proportional hazards assumption...\n")
  
  # Schoenfeld residuals test
  ph_test <- cox.zph(cox_model)
  
  # Create plots
  plot_data <- list()
  
  # Individual variable plots
  for (i in 1:ncol(ph_test$y)) {
    var_name <- colnames(ph_test$y)[i]
    plot_data[[var_name]] <- data.frame(
      time = ph_test$x,
      residual = ph_test$y[, i],
      variable = var_name
    )
  }
  
  results <- list(
    test = ph_test,
    plot_data = plot_data,
    global_p = ph_test$table["GLOBAL", "p"]
  )
  
  cat("Proportional hazards assessment completed\n")
  cat("  Global test p-value:", format.pval(results$global_p, digits = 3), "\n")
  
  return(results)
}

#' Compare models using likelihood ratio test
#'
#' @param cox_results Results from perform_cox_analysis
#' @return Data frame with model comparison results
compare_cox_models <- function(cox_results) {
  
  cat("Comparing Cox models...\n")
  
  # Likelihood ratio tests
  lr_test_1v2 <- anova(cox_results$model_univariate, cox_results$model_full)
  lr_test_1v3 <- anova(cox_results$model_univariate, cox_results$model_simplified)
  
  # AIC comparison
  aic_comparison <- data.frame(
    Model = c("Univariate", "Full", "Simplified"),
    AIC = c(AIC(cox_results$model_univariate),
            AIC(cox_results$model_full),
            AIC(cox_results$model_simplified)),
    Log_likelihood = c(cox_results$model_univariate$loglik[2],
                      cox_results$model_full$loglik[2],
                      cox_results$model_simplified$loglik[2])
  )
  
  # Calculate delta AIC
  aic_comparison$Delta_AIC <- aic_comparison$AIC - min(aic_comparison$AIC)
  
  results <- list(
    lr_test_full = lr_test_1v2,
    lr_test_simplified = lr_test_1v3,
    aic_comparison = aic_comparison
  )
  
  cat("Model comparison completed\n")
  cat("  Best model by AIC:", aic_comparison$Model[which.min(aic_comparison$AIC)], "\n")
  
  return(results)
}

#' Complete cofactors analysis pipeline
#'
#' @param phenotype_data Phenotype data frame
#' @param risk_scores Vector of risk scores
#' @param survival_data Survival data
#' @param output_dir Directory to save results (optional)
#' @return List with complete analysis results
run_cofactors_analysis <- function(phenotype_data, risk_scores, survival_data,
                                  output_dir = NULL) {
  
  cat("Starting complete cofactors analysis pipeline...\n")
  
  # Prepare data
  cox_data <- prepare_cofactors_data(phenotype_data, risk_scores, survival_data)
  
  # Perform Cox analysis
  cox_results <- perform_cox_analysis(cox_data, include_normalized = TRUE)
  
  # Create summary table
  summary_table <- create_cox_summary_table(cox_results)
  
  # Assess proportional hazards
  ph_assessment <- assess_proportional_hazards(cox_results$model_full)
  
  # Compare models
  model_comparison <- compare_cox_models(cox_results)
  
  # Compile results
  results <- list(
    cox_data = cox_data,
    cox_results = cox_results,
    summary_table = summary_table,
    ph_assessment = ph_assessment,
    model_comparison = model_comparison
  )
  
  # Save results if output directory provided
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Save summary table
    write.csv(summary_table, file.path(output_dir, "cox_regression_summary.csv"), 
              row.names = FALSE)
    
    # Save model comparison
    write.csv(model_comparison$aic_comparison, 
              file.path(output_dir, "model_comparison.csv"), row.names = FALSE)
    
    cat("Results saved to:", output_dir, "\n")
  }
  
  cat("Cofactors analysis pipeline completed!\n")
  
  return(results)
}