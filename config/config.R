# =============================================================================
# Configuration File for InferRBPActivityMM Project
# RNA-binding proteins activity inferring model for multiple myeloma risk prediction
# =============================================================================

# Project Information
PROJECT_NAME <- "InferRBPActivityMM"
PROJECT_VERSION <- "1.0.0"
PROJECT_DESCRIPTION <- "Development of RNA-binding proteins activity inferring model for multiple myeloma risk prediction"

# Data Parameters
DATA_CONFIG <- list(
  # Sample sizes
  total_samples = 762,
  train_samples = 600,
  validation_samples = 162,
  
  # Feature dimensions
  rbp_count = 150,
  transcript_count = 49076,
  
  # Thresholds
  correlation_threshold = 0.3,
  psi_threshold = -2,
  
  # External validation datasets
  external_datasets = c("GSE116324", "GSE110486")
)

# Model Parameters
MODEL_CONFIG <- list(
  # Random Forest parameters
  rf_params = list(
    ntree = 50,
    mtry = 15,
    nodesize = 5
  ),
  
  # Cross-validation
  cv_folds = 5,
  
  # Hyperparameter grids
  ntree_grid = c(50, 100, 200, 500, 1000),
  mtry_grid = c(10, 13, 15, 20, 30),
  nodesize_grid = c(5, 10, 20, 50, 100)
)

# File Paths
PATHS <- list(
  # Data directories
  data_raw = "data/raw",
  data_processed = "data/processed", 
  data_results = "data/results",
  
  # Source directories
  src_data = "src/data_processing",
  src_models = "src/models",
  src_analysis = "src/analysis",
  src_viz = "src/visualization",
  src_utils = "src/utils",
  
  # Output directories
  figures = "data/results/figures",
  tables = "data/results/tables",
  models = "data/results/models"
)

# Analysis Parameters
ANALYSIS_CONFIG <- list(
  # Survival analysis
  risk_groups = c("HighRisk", "LowRisk"),
  
  # Statistical tests
  alpha_level = 0.05,
  
  # ISS staging
  iss_stages = c("I", "II", "III"),
  
  # Visualization
  plot_theme = "minimal",
  color_palette = c("lightblue", "green", "red")
)

# Quality Control Parameters
QC_CONFIG <- list(
  # Data filtering
  min_iqr_percentile = 0.25,
  max_outlier_threshold = 40,
  
  # Sample filtering
  min_survival_time = 0,
  
  # Missing data handling
  missing_data_action = "na.impute"
)

# Logging Configuration
LOG_CONFIG <- list(
  log_level = "INFO",
  log_file = "logs/analysis.log",
  verbose = TRUE
)

# Function to create directories if they don't exist
create_project_dirs <- function() {
  dirs_to_create <- c(
    PATHS$data_raw, PATHS$data_processed, PATHS$data_results,
    PATHS$figures, PATHS$tables, PATHS$models,
    "logs"
  )
  
  for (dir in dirs_to_create) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("Created directory:", dir, "\n")
    }
  }
}

# Function to validate configuration
validate_config <- function() {
  # Check if required packages are installed
  required_packages <- c(
    "survival", "survminer", "randomForestSRC", 
    "ggplot2", "ArrayTools"
  )
  
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    warning("Missing required packages: ", paste(missing_packages, collapse = ", "))
    return(FALSE)
  }
  
  return(TRUE)
}

# Initialize project
init_project <- function() {
  cat("Initializing InferRBPActivityMM project...\n")
  create_project_dirs()
  
  if (validate_config()) {
    cat("Project initialized successfully!\n")
  } else {
    stop("Project initialization failed. Please install missing packages.")
  }
}