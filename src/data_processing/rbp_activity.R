# =============================================================================
# RBP Activity Calculation Module
# =============================================================================

source("config/config.R")
source("src/utils/data_utils.R")

#' Calculate RBP activity scores
#'
#' @param target_matrix Target matrix (M_IQR)
#' @param PSI_scores Centered PSI scores
#' @param psi_threshold Threshold for valid PSI values
#' @return Matrix of RBP activity scores
calculate_rbp_activity <- function(target_matrix, PSI_scores, psi_threshold = -2) {
  
  cat("Calculating RBP activity scores...\n")
  
  rbp_count <- nrow(target_matrix)
  sample_count <- ncol(PSI_scores)
  
  # Initialize activity matrix
  RBP_activity <- matrix(0, rbp_count, sample_count)
  rownames(RBP_activity) <- rownames(target_matrix)
  colnames(RBP_activity) <- colnames(PSI_scores)
  
  for (i in 1:rbp_count) {
    for (j in 1:sample_count) {
      # Find positive and negative targets with valid PSI values
      positive_targets <- intersect(
        which(target_matrix[i, ] == 1),
        which(PSI_scores[, j] > psi_threshold)
      )
      
      negative_targets <- intersect(
        which(target_matrix[i, ] == -1),
        which(PSI_scores[, j] > psi_threshold)
      )
      
      # Calculate activity score
      if (length(positive_targets) + length(negative_targets) > 0) {
        positive_score <- ifelse(length(positive_targets) > 0, 
                               sum(PSI_scores[positive_targets, j]), 0)
        negative_score <- ifelse(length(negative_targets) > 0,
                               sum(PSI_scores[negative_targets, j]), 0)
        
        RBP_activity[i, j] <- (positive_score - negative_score) / 
                             (length(positive_targets) + length(negative_targets))
      }
    }
    
    if (i %% 10 == 0) cat("Calculated activity for RBP", i, "/", rbp_count, "\n")
  }
  
  cat("RBP activity calculation completed!\n")
  return(RBP_activity)
}

#' Process complete RBP activity pipeline for training data
#'
#' @param data_split Result from split_data function
#' @param PSI_matrix Raw PSI matrix
#' @param RBP_matrix RBP expression matrix
#' @param correlation_matrix ENCODE correlation matrix
#' @param transcript_names Vector of transcript names
#' @param correlation_threshold Threshold for correlation filtering
#' @return List with processed training and validation data
process_rbp_activity_pipeline <- function(data_split, PSI_matrix, RBP_matrix, 
                                         correlation_matrix, transcript_names,
                                         correlation_threshold = 0.3) {
  
  cat("Starting RBP activity processing pipeline...\n")
  
  # Split data
  train_indices <- data_split$train
  validation_indices <- data_split$validation
  
  RBP_train <- RBP_matrix[, train_indices]
  PSI_train <- PSI_matrix[, train_indices]
  RBP_validation <- RBP_matrix[, validation_indices]
  PSI_validation <- PSI_matrix[, validation_indices]
  
  # Create CLIP target matrix
  M_CLIP <- create_clip_matrix(
    correlation_matrix = correlation_matrix,
    transcript_names = transcript_names,
    rbp_count = DATA_CONFIG$rbp_count,
    transcript_count = DATA_CONFIG$transcript_count
  )
  
  # Create correlation-based target matrix
  M_corr <- create_correlation_matrix(
    M_CLIP = M_CLIP,
    RBP_data = RBP_train,
    PSI_data = PSI_train,
    threshold = correlation_threshold
  )
  
  # Calculate PSI metrics and apply IQR filter
  psi_metrics <- calculate_psi_metrics(PSI_train)
  iqr_threshold <- quantile(psi_metrics$iqr, QC_CONFIG$min_iqr_percentile)
  
  M_IQR <- M_corr
  M_IQR[, psi_metrics$iqr < iqr_threshold] <- 0
  
  cat("Applied IQR filter. Retained", sum(apply(abs(M_IQR), 2, sum) > 0), 
      "transcripts out of", ncol(M_IQR), "\n")
  
  # Normalize and center training PSI data
  norm_result_train <- normalize_psi_data(PSI_train)
  center_result_train <- center_psi_scores(norm_result_train$normalized_data)
  
  # Normalize and center validation PSI data using training parameters
  norm_result_validation <- normalize_psi_data(
    PSI_validation,
    min_values = norm_result_train$min_values,
    range_diff = norm_result_train$range_diff
  )
  center_result_validation <- center_psi_scores(
    norm_result_validation$normalized_data,
    median_values = center_result_train$median_values
  )
  
  # Calculate RBP activities
  RBP_train_activity <- calculate_rbp_activity(
    target_matrix = M_IQR,
    PSI_scores = center_result_train$centered_scores,
    psi_threshold = DATA_CONFIG$psi_threshold
  )
  
  RBP_validation_activity <- calculate_rbp_activity(
    target_matrix = M_IQR,
    PSI_scores = center_result_validation$centered_scores,
    psi_threshold = DATA_CONFIG$psi_threshold
  )
  
  # Return comprehensive results
  list(
    # Processed data
    RBP_train_activity = RBP_train_activity,
    RBP_validation_activity = RBP_validation_activity,
    
    # Intermediate matrices
    M_CLIP = M_CLIP,
    M_correlation = M_corr,
    M_IQR = M_IQR,
    
    # Normalization parameters
    normalization_params = list(
      min_values = norm_result_train$min_values,
      range_diff = norm_result_train$range_diff,
      median_values = center_result_train$median_values,
      iqr_threshold = iqr_threshold
    ),
    
    # Data splits
    data_indices = list(
      train = train_indices,
      validation = validation_indices
    )
  )
}

#' Process external validation data
#'
#' @param PSI_external External PSI matrix
#' @param normalization_params Parameters from training data
#' @param M_IQR Target matrix from training
#' @return RBP activity matrix for external data
process_external_validation <- function(PSI_external, normalization_params, M_IQR) {
  
  cat("Processing external validation data...\n")
  
  # Normalize using training parameters
  norm_result <- normalize_psi_data(
    PSI_external,
    min_values = normalization_params$min_values,
    range_diff = normalization_params$range_diff
  )
  
  # Center using training parameters
  center_result <- center_psi_scores(
    norm_result$normalized_data,
    median_values = normalization_params$median_values
  )
  
  # Calculate RBP activity
  RBP_external_activity <- calculate_rbp_activity(
    target_matrix = M_IQR,
    PSI_scores = center_result$centered_scores,
    psi_threshold = DATA_CONFIG$psi_threshold
  )
  
  return(RBP_external_activity)
}

#' Validate RBP activity results
#'
#' @param RBP_activity RBP activity matrix
#' @param expected_samples Expected number of samples
#' @param expected_rbps Expected number of RBPs
#' @return Boolean indicating validation success
validate_rbp_activity <- function(RBP_activity, expected_samples = NULL, expected_rbps = NULL) {
  
  cat("Validating RBP activity results...\n")
  
  # Check dimensions
  if (!is.null(expected_rbps) && nrow(RBP_activity) != expected_rbps) {
    warning("Unexpected number of RBPs: ", nrow(RBP_activity), " vs ", expected_rbps)
    return(FALSE)
  }
  
  if (!is.null(expected_samples) && ncol(RBP_activity) != expected_samples) {
    warning("Unexpected number of samples: ", ncol(RBP_activity), " vs ", expected_samples)
    return(FALSE)
  }
  
  # Check for excessive missing values
  na_proportion <- sum(is.na(RBP_activity)) / length(RBP_activity)
  if (na_proportion > 0.5) {
    warning("High proportion of missing values: ", round(na_proportion * 100, 2), "%")
    return(FALSE)
  }
  
  # Check value ranges
  finite_values <- RBP_activity[is.finite(RBP_activity)]
  if (length(finite_values) == 0) {
    warning("No finite values in RBP activity matrix")
    return(FALSE)
  }
  
  value_range <- range(finite_values)
  cat("RBP activity range:", paste(round(value_range, 4), collapse = " to "), "\n")
  cat("Missing values:", sum(is.na(RBP_activity)), "/", length(RBP_activity), "\n")
  
  cat("RBP activity validation passed!\n")
  return(TRUE)
}