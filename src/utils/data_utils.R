# =============================================================================
# Data Utility Functions for InferRBPActivityMM Project
# =============================================================================

#' Load and validate R data file
#'
#' @param file_path Path to the .Rdata file
#' @param expected_vars Expected variable names (optional)
#' @return Environment with loaded variables
load_rdata_safely <- function(file_path, expected_vars = NULL) {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  env <- new.env()
  load(file_path, envir = env)
  
  if (!is.null(expected_vars)) {
    missing_vars <- setdiff(expected_vars, ls(env))
    if (length(missing_vars) > 0) {
      warning("Missing expected variables: ", paste(missing_vars, collapse = ", "))
    }
  }
  
  cat("Loaded data from:", file_path, "\n")
  cat("Variables:", paste(ls(env), collapse = ", "), "\n")
  
  return(env)
}

#' Split data into training and validation sets
#'
#' @param total_samples Total number of samples
#' @param train_size Number of training samples
#' @param seed Random seed for reproducibility
#' @return List with train and validation indices
split_data <- function(total_samples, train_size, seed = 123) {
  set.seed(seed)
  train_indices <- sample(1:total_samples, train_size)
  validation_indices <- setdiff(1:total_samples, train_indices)
  
  list(
    train = train_indices,
    validation = validation_indices
  )
}

#' Create target matrix from correlation data
#'
#' @param correlation_matrix Matrix of correlations
#' @param transcript_names Vector of transcript names
#' @param rbp_count Number of RBPs
#' @param transcript_count Number of transcripts
#' @return Binary target matrix
create_clip_matrix <- function(correlation_matrix, transcript_names, 
                              rbp_count, transcript_count) {
  
  cat("Creating CLIP target matrix...\n")
  
  # Initialize matrix
  M_CLIP <- matrix(0, rbp_count, transcript_count)
  rownames(M_CLIP) <- rownames(correlation_matrix)
  colnames(M_CLIP) <- transcript_names
  
  # Fill matrix based on correlations
  for (i in 1:rbp_count) {
    encode_targets <- colnames(correlation_matrix)[which(correlation_matrix[i, ] != 0)]
    
    for (j in 1:transcript_count) {
      if (transcript_names[j] %in% encode_targets) {
        M_CLIP[i, j] <- 1
      }
    }
    
    if (i %% 10 == 0) cat("Processed RBP", i, "/", rbp_count, "\n")
  }
  
  target_numbers <- apply(M_CLIP, 1, sum)
  cat("Target numbers - Median:", median(target_numbers), 
      "Range:", range(target_numbers), "\n")
  
  return(M_CLIP)
}

#' Calculate Spearman correlations and create final target matrix
#'
#' @param M_CLIP Binary CLIP matrix
#' @param RBP_data RBP expression data
#' @param PSI_data PSI data
#' @param threshold Correlation threshold
#' @return Final target matrix with correlation-based filtering
create_correlation_matrix <- function(M_CLIP, RBP_data, PSI_data, threshold = 0.3) {
  
  cat("Creating correlation-based target matrix...\n")
  
  rbp_count <- nrow(M_CLIP)
  transcript_count <- ncol(M_CLIP)
  
  # Initialize final matrix
  M_final <- matrix(0, rbp_count, transcript_count)
  rownames(M_final) <- rownames(M_CLIP)
  colnames(M_final) <- colnames(M_CLIP)
  
  for (i in 1:rbp_count) {
    target_indices <- which(M_CLIP[i, ] == 1)
    
    if (length(target_indices) > 0) {
      correlations <- numeric(length(target_indices))
      
      for (j in seq_along(target_indices)) {
        idx <- target_indices[j]
        valid_samples <- which(PSI_data[idx, ] >= 0)
        
        if (length(valid_samples) > 3) {  # Need at least 4 samples for correlation
          correlation <- cor(RBP_data[i, valid_samples], 
                           PSI_data[idx, valid_samples], 
                           method = "spearman")
          correlations[j] <- ifelse(is.na(correlation), 0, correlation)
        }
      }
      
      # Set targets based on correlation threshold
      positive_targets <- target_indices[correlations >= threshold]
      negative_targets <- target_indices[correlations <= -threshold]
      
      M_final[i, positive_targets] <- 1
      M_final[i, negative_targets] <- -1
    }
    
    if (i %% 10 == 0) cat("Processed correlations for RBP", i, "/", rbp_count, "\n")
  }
  
  final_targets <- apply(abs(M_final), 1, sum)
  cat("Final target numbers - Median:", median(final_targets), 
      "Range:", range(final_targets), "\n")
  
  return(M_final)
}

#' Calculate statistical metrics for PSI data
#'
#' @param PSI_data PSI matrix
#' @return List with range differences and IQR values
calculate_psi_metrics <- function(PSI_data) {
  
  cat("Calculating PSI metrics...\n")
  
  transcript_count <- nrow(PSI_data)
  D <- numeric(transcript_count)
  IQR_values <- numeric(transcript_count)
  
  for (i in 1:transcript_count) {
    valid_values <- PSI_data[i, ][PSI_data[i, ] >= 0]
    
    if (length(valid_values) > 0) {
      D[i] <- diff(range(valid_values))
      IQR_values[i] <- IQR(valid_values)
    }
    
    if (i %% 5000 == 0) cat("Processed", i, "/", transcript_count, "transcripts\n")
  }
  
  list(
    range_diff = D,
    iqr = IQR_values
  )
}

#' Normalize PSI data
#'
#' @param PSI_data Raw PSI data
#' @param min_values Minimum values for normalization
#' @param range_diff Range differences for normalization
#' @return Normalized PSI data
normalize_psi_data <- function(PSI_data, min_values = NULL, range_diff = NULL) {
  
  cat("Normalizing PSI data...\n")
  
  PSI_norm <- apply(PSI_data, 2, as.numeric)
  rownames(PSI_norm) <- rownames(PSI_data)
  transcript_count <- nrow(PSI_norm)
  
  # Calculate normalization parameters if not provided
  if (is.null(min_values) || is.null(range_diff)) {
    metrics <- calculate_psi_metrics(PSI_data)
    
    if (is.null(min_values)) {
      min_values <- numeric(transcript_count)
      for (i in 1:transcript_count) {
        valid_values <- PSI_norm[i, ][PSI_norm[i, ] >= 0]
        if (length(valid_values) > 0) {
          min_values[i] <- min(valid_values, na.rm = TRUE)
        }
      }
    }
    
    if (is.null(range_diff)) {
      range_diff <- metrics$range_diff
    }
  }
  
  # Normalize
  for (i in 1:transcript_count) {
    valid_indices <- which(PSI_norm[i, ] >= 0)
    if (length(valid_indices) > 0 && range_diff[i] > 0) {
      PSI_norm[i, valid_indices] <- (PSI_norm[i, valid_indices] - min_values[i]) / range_diff[i]
    }
    
    if (i %% 5000 == 0) cat("Normalized", i, "/", transcript_count, "transcripts\n")
  }
  
  list(
    normalized_data = PSI_norm,
    min_values = min_values,
    range_diff = range_diff
  )
}

#' Center PSI scores using median
#'
#' @param PSI_norm Normalized PSI data
#' @param median_values Pre-calculated median values (optional)
#' @return List with centered scores and median values
center_psi_scores <- function(PSI_norm, median_values = NULL) {
  
  cat("Centering PSI scores...\n")
  
  PSI_score <- PSI_norm
  transcript_count <- nrow(PSI_score)
  
  if (is.null(median_values)) {
    median_values <- numeric(transcript_count)
    
    for (i in 1:transcript_count) {
      valid_scores <- PSI_score[i, ][!is.na(PSI_score[i, ])]
      if (length(valid_scores) > 0) {
        median_values[i] <- median(valid_scores)
      }
    }
  }
  
  # Center scores
  for (i in 1:transcript_count) {
    valid_indices <- which(!is.na(PSI_score[i, ]))
    if (length(valid_indices) > 0) {
      PSI_score[i, valid_indices] <- PSI_score[i, valid_indices] - median_values[i]
    }
    
    if (i %% 5000 == 0) cat("Centered", i, "/", transcript_count, "transcripts\n")
  }
  
  list(
    centered_scores = PSI_score,
    median_values = median_values
  )
}

#' Validate data dimensions and consistency
#'
#' @param ... Named data objects to validate
#' @return Boolean indicating if validation passed
validate_data_consistency <- function(...) {
  data_objects <- list(...)
  
  cat("Validating data consistency...\n")
  
  # Check if all objects exist
  for (name in names(data_objects)) {
    if (is.null(data_objects[[name]])) {
      warning("Data object '", name, "' is NULL")
      return(FALSE)
    }
  }
  
  # Add more specific validation logic as needed
  cat("Data validation passed!\n")
  return(TRUE)
}