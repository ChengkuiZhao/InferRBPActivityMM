#!/usr/bin/env Rscript

# =============================================================================
# Main Analysis Script for InferRBPActivityMM Project
# RNA-binding proteins activity inferring model for multiple myeloma risk prediction
# =============================================================================

# Load configuration and all modules
source("config/config.R")
source("src/utils/data_utils.R")
source("src/data_processing/rbp_activity.R")
source("src/models/survival_model.R")
source("src/analysis/cofactors_analysis.R")
source("src/visualization/plots.R")

# Initialize project
cat("=== Initializing InferRBPActivityMM Analysis ===\n")
init_project()
setup_plotting_environment(ANALYSIS_CONFIG$plot_theme)

# =============================================================================
# STEP 1: Load and Prepare Data
# =============================================================================

cat("\n=== STEP 1: Loading and Preparing Data ===\n")

# Load raw data (adjust paths as needed)
raw_data_dir <- "code/MMRF Analysis/1.Survival"

# Load PSI matrix
cat("Loading PSI matrix...\n")
psi_env <- load_rdata_safely(file.path(raw_data_dir, "PSImatrix_selected.Rdata"))
PSI_matrix <- psi_env$PSImatrix_selected

# Load RBP data
cat("Loading RBP data...\n")
rbp_env <- load_rdata_safely(file.path(raw_data_dir, "RBP_selected.Rdata"))
RBP_matrix <- rbp_env$RBP_selected

# Load survival data
cat("Loading survival data...\n")
surv_env <- load_rdata_safely(file.path(raw_data_dir, "Survival_selected.Rdata"))
survival_data <- surv_env$Survival_selected

# Load phenotype data for cofactors analysis
cat("Loading phenotype data...\n")
pheno_env <- load_rdata_safely(file.path(raw_data_dir, "Phenotype_selected.Rdata"))
phenotype_data <- pheno_env$Phenotype_selected

# Load correlation matrix
cat("Loading correlation matrix...\n")
correlation_matrix <- read.csv(file.path(raw_data_dir, "CorrelationMatrix"))
rownames(correlation_matrix) <- correlation_matrix[, 1]
correlation_matrix <- t(as.matrix(correlation_matrix[, -1]))

# Load transcript names
cat("Loading transcript names...\n")
transcript_env <- load_rdata_safely(file.path(raw_data_dir, "Transcriptname.Rdata"))
transcript_names <- transcript_env$ENSGname

# Validate data consistency
validate_data_consistency(
  PSI_matrix = PSI_matrix,
  RBP_matrix = RBP_matrix,
  survival_data = survival_data,
  phenotype_data = phenotype_data
)

# =============================================================================
# STEP 2: RBP Activity Calculation
# =============================================================================

cat("\n=== STEP 2: RBP Activity Calculation ===\n")

# Split data into training and validation sets
data_split <- split_data(
  total_samples = DATA_CONFIG$total_samples,
  train_size = DATA_CONFIG$train_samples,
  seed = 123
)

# Process RBP activity pipeline
rbp_results <- process_rbp_activity_pipeline(
  data_split = data_split,
  PSI_matrix = PSI_matrix,
  RBP_matrix = RBP_matrix,
  correlation_matrix = correlation_matrix,
  transcript_names = transcript_names,
  correlation_threshold = DATA_CONFIG$correlation_threshold
)

# Validate RBP activity results
validate_rbp_activity(
  RBP_activity = rbp_results$RBP_train_activity,
  expected_samples = DATA_CONFIG$train_samples,
  expected_rbps = DATA_CONFIG$rbp_count
)

validate_rbp_activity(
  RBP_activity = rbp_results$RBP_validation_activity,
  expected_samples = DATA_CONFIG$validation_samples,
  expected_rbps = DATA_CONFIG$rbp_count
)

# =============================================================================
# STEP 3: Survival Analysis and Model Training
# =============================================================================

cat("\n=== STEP 3: Survival Analysis and Model Training ===\n")

# Prepare survival data for training and validation
train_indices <- rbp_results$data_indices$train
validation_indices <- rbp_results$data_indices$validation

survival_train <- survival_data[train_indices, ]
survival_validation <- survival_data[validation_indices, ]

# Run complete survival analysis
survival_results <- run_survival_analysis(
  RBP_activity_train = rbp_results$RBP_train_activity,
  RBP_activity_validation = rbp_results$RBP_validation_activity,
  survival_data_train = survival_train,
  survival_data_validation = survival_validation,
  tune_hyperparams = TRUE  # Set to FALSE for faster execution
)

# =============================================================================
# STEP 4: Cofactors Analysis
# =============================================================================

cat("\n=== STEP 4: Cofactors Analysis ===\n")

# Run cofactors analysis on validation set
validation_sample_names <- colnames(rbp_results$RBP_validation_activity)

cofactors_results <- run_cofactors_analysis(
  phenotype_data = phenotype_data,
  risk_scores = survival_results$risk_scores,
  survival_data = survival_validation,
  output_dir = PATHS$tables
)

# =============================================================================
# STEP 5: Visualization and Reporting
# =============================================================================

cat("\n=== STEP 5: Visualization and Reporting ===\n")

# Create comprehensive analysis report
analysis_plots <- create_analysis_report(
  analysis_results = c(survival_results, rbp_results),
  dataset_name = "MMRF_Validation",
  output_dir = PATHS$figures
)

# Create ISS stage analysis plot
if ("demographic.iss_stage" %in% colnames(phenotype_data)) {
  iss_stages <- phenotype_data[validation_sample_names, "demographic.iss_stage"]
  
  iss_plot <- plot_iss_stage_comparison(
    risk_scores = survival_results$risk_scores,
    iss_stages = iss_stages,
    dataset_name = "MMRF_Validation",
    save_path = file.path(PATHS$figures, "MMRF_ISS_stage_comparison.png")
  )
}

# =============================================================================
# STEP 6: External Validation (if data available)
# =============================================================================

cat("\n=== STEP 6: External Validation ===\n")

# Process external validation datasets
external_datasets <- c("GSE116324", "GSE110486")

for (dataset in external_datasets) {
  cat("\nProcessing external dataset:", dataset, "\n")
  
  external_dir <- file.path("code/External Validation", dataset)
  
  tryCatch({
    # Load external PSI data
    psi_file <- file.path(external_dir, paste0("PSImatrix_", dataset, ".Rdata"))
    if (file.exists(psi_file)) {
      external_psi_env <- load_rdata_safely(psi_file)
      external_psi <- external_psi_env$PSImatrix
      
      # Process external validation
      external_activity <- process_external_validation(
        PSI_external = external_psi,
        normalization_params = rbp_results$normalization_params,
        M_IQR = rbp_results$M_IQR
      )
      
      # Predict risk scores
      external_risk_scores <- predict_risk_scores(
        model = survival_results$model,
        RBP_activity_test = external_activity
      )
      
      # Save external results
      external_results <- list(
        dataset = dataset,
        risk_scores = external_risk_scores,
        RBP_activity = external_activity
      )
      
      save(external_results, 
           file = file.path(PATHS$models, paste0(dataset, "_results.Rdata")))
      
      cat("External validation completed for", dataset, "\n")
    } else {
      cat("PSI data not found for", dataset, "\n")
    }
    
  }, error = function(e) {
    cat("Error processing", dataset, ":", e$message, "\n")
  })
}

# =============================================================================
# STEP 7: Save Final Results
# =============================================================================

cat("\n=== STEP 7: Saving Final Results ===\n")

# Save main analysis results
main_results <- list(
  rbp_results = rbp_results,
  survival_results = survival_results,
  cofactors_results = cofactors_results,
  config = list(
    DATA_CONFIG = DATA_CONFIG,
    MODEL_CONFIG = MODEL_CONFIG,
    ANALYSIS_CONFIG = ANALYSIS_CONFIG
  )
)

save(main_results, file = file.path(PATHS$models, "main_analysis_results.Rdata"))

# Save risk scores for validation samples
risk_scores_df <- data.frame(
  sample_id = names(survival_results$risk_scores),
  risk_score = survival_results$risk_scores,
  risk_group = survival_results$risk_groups,
  stringsAsFactors = FALSE
)

write.csv(risk_scores_df, 
          file = file.path(PATHS$tables, "validation_risk_scores.csv"),
          row.names = FALSE)

# Create final summary report
cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Project:", PROJECT_NAME, "v", PROJECT_VERSION, "\n")
cat("Training samples:", DATA_CONFIG$train_samples, "\n")
cat("Validation samples:", DATA_CONFIG$validation_samples, "\n")
cat("RBP features:", DATA_CONFIG$rbp_count, "\n")
cat("Best model parameters:\n")
cat("  ntree:", survival_results$best_params$ntree, "\n")
cat("  mtry:", survival_results$best_params$mtry, "\n")
cat("  nodesize:", survival_results$best_params$nodesize, "\n")
cat("Model performance:\n")
cat("  Log-rank p-value:", format.pval(survival_results$performance$logrank_pvalue), "\n")
cat("  C-index:", round(survival_results$performance$c_index, 3), "\n")

# Print file locations
cat("\nOutput files saved to:\n")
cat("  Figures:", PATHS$figures, "\n")
cat("  Tables:", PATHS$tables, "\n")
cat("  Models:", PATHS$models, "\n")

cat("\n=== Analysis Completed Successfully! ===\n")

# Optional: Generate HTML report
if (requireNamespace("rmarkdown", quietly = TRUE)) {
  cat("Generating HTML report...\n")
  # Code to generate report would go here
}