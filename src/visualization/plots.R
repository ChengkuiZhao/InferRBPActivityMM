# =============================================================================
# Visualization Module for InferRBPActivityMM Project
# =============================================================================

source("config/config.R")

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(survival)
  library(survminer)
  library(gridExtra)
  library(RColorBrewer)
})

#' Create survival curve plot
#'
#' @param survival_fit survfit object from survival analysis
#' @param data_frame Data frame with survival data
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @param save_path Path to save the plot (optional)
#' @return ggplot object
plot_survival_curves <- function(survival_fit, data_frame, 
                                title = "Kaplan-Meier Survival Curves",
                                subtitle = NULL,
                                save_path = NULL) {
  
  cat("Creating survival curve plot...\n")
  
  # Create survival plot
  survival_plot <- ggsurvplot(
    survival_fit,
    data = data_frame,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.col = "strata",
    linetype = "strata",
    surv.median.line = "hv",
    ggtheme = theme_minimal(),
    palette = ANALYSIS_CONFIG$color_palette[1:length(levels(data_frame$risk_group))],
    title = title,
    subtitle = subtitle,
    xlab = "Time (days)",
    ylab = "Overall Survival Probability",
    legend.title = "Risk Group",
    legend.labs = levels(data_frame$risk_group)
  )
  
  # Save plot if path provided
  if (!is.null(save_path)) {
    ggsave(save_path, plot = survival_plot$plot, 
           width = 10, height = 8, dpi = 300)
    cat("Survival plot saved to:", save_path, "\n")
  }
  
  return(survival_plot)
}

#' Create risk score distribution plot
#'
#' @param risk_scores Vector of risk scores
#' @param risk_groups Factor of risk group classifications
#' @param title Plot title
#' @param save_path Path to save the plot (optional)
#' @return ggplot object
plot_risk_distribution <- function(risk_scores, risk_groups,
                                  title = "Risk Score Distribution",
                                  save_path = NULL) {
  
  cat("Creating risk score distribution plot...\n")
  
  # Prepare data
  plot_data <- data.frame(
    risk_score = risk_scores,
    risk_group = risk_groups
  )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = risk_score, fill = risk_group)) +
    geom_histogram(alpha = 0.7, bins = 30, position = "identity") +
    geom_vline(xintercept = median(risk_scores), 
               linetype = "dashed", color = "red", size = 1) +
    scale_fill_manual(values = ANALYSIS_CONFIG$color_palette) +
    theme_minimal() +
    labs(
      title = title,
      x = "Risk Score",
      y = "Frequency",
      fill = "Risk Group"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "top"
    )
  
  # Save plot if path provided
  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = 10, height = 6, dpi = 300)
    cat("Risk distribution plot saved to:", save_path, "\n")
  }
  
  return(p)
}

#' Create ISS stage comparison boxplot
#'
#' @param risk_scores Vector of risk scores
#' @param iss_stages Vector of ISS stages
#' @param dataset_name Name of the dataset
#' @param save_path Path to save the plot (optional)
#' @return ggplot object
plot_iss_stage_comparison <- function(risk_scores, iss_stages, 
                                     dataset_name = "Dataset",
                                     save_path = NULL) {
  
  cat("Creating ISS stage comparison plot...\n")
  
  # Clean and prepare data
  valid_indices <- !is.na(risk_scores) & !is.na(iss_stages) & iss_stages != ""
  
  plot_data <- data.frame(
    risk_score = risk_scores[valid_indices],
    iss_stage = factor(iss_stages[valid_indices], levels = c("I", "II", "III"))
  )
  
  # Remove outliers for better visualization
  outlier_threshold <- quantile(plot_data$risk_score, 0.95, na.rm = TRUE)
  plot_data_clean <- plot_data[plot_data$risk_score <= outlier_threshold, ]
  
  # Statistical test
  stage_12 <- plot_data_clean$risk_score[plot_data_clean$iss_stage %in% c("I", "II")]
  stage_3 <- plot_data_clean$risk_score[plot_data_clean$iss_stage == "III"]
  
  if (length(stage_12) > 0 && length(stage_3) > 0) {
    t_test_result <- t.test(stage_12, stage_3, alternative = "less")
    p_value <- t_test_result$p.value
  } else {
    p_value <- NA
  }
  
  # Create main plot
  p1 <- ggplot(plot_data_clean, aes(x = iss_stage, y = risk_score, fill = iss_stage)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    scale_fill_manual(values = ANALYSIS_CONFIG$color_palette) +
    scale_x_discrete(labels = paste("Stage", c("I", "II", "III"))) +
    theme_minimal() +
    labs(
      title = paste("Risk Score by ISS Stage -", dataset_name),
      x = "ISS Stage",
      y = "Risk Score",
      fill = "ISS Stage"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "none"
    )
  
  # Add p-value annotation if available
  if (!is.na(p_value)) {
    p1 <- p1 + annotate("text", x = 2, y = max(plot_data_clean$risk_score) * 0.9,
                       label = paste("p =", format.pval(p_value, digits = 3)),
                       size = 4, fontface = "bold")
  }
  
  # Create binary comparison plot (Stage I&II vs III)
  plot_data_binary <- plot_data_clean
  plot_data_binary$stage_binary <- ifelse(plot_data_binary$iss_stage == "III", 
                                         "Stage III", "Stage I & II")
  
  p2 <- ggplot(plot_data_binary, aes(x = stage_binary, y = risk_score, fill = stage_binary)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    scale_fill_manual(values = ANALYSIS_CONFIG$color_palette[c(1, 3)]) +
    theme_minimal() +
    labs(
      title = paste("Risk Score Binary Comparison -", dataset_name),
      subtitle = ifelse(!is.na(p_value), paste("t-test p-value =", format.pval(p_value)), ""),
      x = "ISS Stage Group",
      y = "Risk Score",
      fill = "Stage Group"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "none"
    )
  
  # Combine plots
  combined_plot <- grid.arrange(p1, p2, ncol = 2, 
                               top = paste("ISS Stage Analysis -", dataset_name))
  
  # Save plot if path provided
  if (!is.null(save_path)) {
    ggsave(save_path, plot = combined_plot, width = 14, height = 6, dpi = 300)
    cat("ISS stage comparison plot saved to:", save_path, "\n")
  }
  
  return(list(individual = p1, binary = p2, combined = combined_plot))
}

#' Create target number visualization
#'
#' @param target_matrix Target matrix (e.g., M_CLIP, M_correlation)
#' @param title Plot title
#' @param save_path Path to save the plot (optional)
#' @return ggplot object
plot_target_numbers <- function(target_matrix, title = "RBP Target Numbers",
                                save_path = NULL) {
  
  cat("Creating target numbers plot...\n")
  
  # Calculate target numbers
  target_numbers <- apply(abs(target_matrix), 1, sum)
  
  # Create data frame
  plot_data <- data.frame(
    RBP = factor(1:length(target_numbers)),
    targets = target_numbers
  )
  
  # Create boxplot
  p <- ggplot(plot_data, aes(y = targets)) +
    geom_boxplot(fill = ANALYSIS_CONFIG$color_palette[1], alpha = 0.7, width = 0.5) +
    geom_text(aes(x = 0, y = median(targets)), 
              label = round(median(target_numbers), 0),
              hjust = -0.5, vjust = -0.5, size = 4, fontface = "bold") +
    theme_minimal() +
    labs(
      title = title,
      y = "Number of Targets",
      x = ""
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  # Save plot if path provided
  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = 8, height = 6, dpi = 300)
    cat("Target numbers plot saved to:", save_path, "\n")
  }
  
  return(p)
}

#' Create correlation heatmap
#'
#' @param correlation_matrix Matrix of correlations
#' @param title Plot title
#' @param max_features Maximum number of features to display
#' @param save_path Path to save the plot (optional)
#' @return ggplot object
plot_correlation_heatmap <- function(correlation_matrix, 
                                    title = "RBP-Target Correlations",
                                    max_features = 50,
                                    save_path = NULL) {
  
  cat("Creating correlation heatmap...\n")
  
  # Sample features if matrix is too large
  if (ncol(correlation_matrix) > max_features) {
    selected_cols <- sample(ncol(correlation_matrix), max_features)
    correlation_matrix <- correlation_matrix[, selected_cols]
  }
  
  if (nrow(correlation_matrix) > max_features) {
    selected_rows <- sample(nrow(correlation_matrix), max_features)
    correlation_matrix <- correlation_matrix[selected_rows, ]
  }
  
  # Convert to long format
  library(reshape2)
  corr_long <- melt(correlation_matrix)
  colnames(corr_long) <- c("RBP", "Target", "Correlation")
  
  # Create heatmap
  p <- ggplot(corr_long, aes(x = Target, y = RBP, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                        midpoint = 0, limit = c(-1, 1)) +
    theme_minimal() +
    labs(
      title = title,
      x = "Targets",
      y = "RBPs"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
    )
  
  # Save plot if path provided
  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = 10, height = 8, dpi = 300)
    cat("Correlation heatmap saved to:", save_path, "\n")
  }
  
  return(p)
}

#' Create comprehensive analysis report plots
#'
#' @param analysis_results Complete analysis results
#' @param dataset_name Name of the dataset
#' @param output_dir Directory to save plots
#' @return List of all generated plots
create_analysis_report <- function(analysis_results, dataset_name = "Analysis",
                                  output_dir = NULL) {
  
  cat("Creating comprehensive analysis report...\n")
  
  if (!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  plots <- list()
  
  # Survival curves
  if (!is.null(analysis_results$performance$km_fit)) {
    survival_path <- if (!is.null(output_dir)) {
      file.path(output_dir, paste0(dataset_name, "_survival_curves.png"))
    } else NULL
    
    plots$survival <- plot_survival_curves(
      survival_fit = analysis_results$performance$km_fit,
      data_frame = analysis_results$performance$eval_data,
      title = paste("Survival Analysis -", dataset_name),
      save_path = survival_path
    )
  }
  
  # Risk distribution
  if (!is.null(analysis_results$risk_scores)) {
    risk_dist_path <- if (!is.null(output_dir)) {
      file.path(output_dir, paste0(dataset_name, "_risk_distribution.png"))
    } else NULL
    
    plots$risk_distribution <- plot_risk_distribution(
      risk_scores = analysis_results$risk_scores,
      risk_groups = analysis_results$risk_groups,
      title = paste("Risk Score Distribution -", dataset_name),
      save_path = risk_dist_path
    )
  }
  
  # Target numbers (if available)
  if (!is.null(analysis_results$M_IQR)) {
    target_path <- if (!is.null(output_dir)) {
      file.path(output_dir, paste0(dataset_name, "_target_numbers.png"))
    } else NULL
    
    plots$target_numbers <- plot_target_numbers(
      target_matrix = analysis_results$M_IQR,
      title = paste("RBP Target Numbers -", dataset_name),
      save_path = target_path
    )
  }
  
  cat("Analysis report plots created successfully!\n")
  
  return(plots)
}

#' Setup plotting environment
#'
#' @param theme_name ggplot theme name
setup_plotting_environment <- function(theme_name = "minimal") {
  
  # Set global ggplot theme
  if (theme_name == "minimal") {
    theme_set(theme_minimal())
  } else if (theme_name == "classic") {
    theme_set(theme_classic())
  } else if (theme_name == "bw") {
    theme_set(theme_bw())
  }
  
  # Set color palette
  options(ggplot2.discrete.colour = ANALYSIS_CONFIG$color_palette)
  options(ggplot2.discrete.fill = ANALYSIS_CONFIG$color_palette)
  
  cat("Plotting environment configured with theme:", theme_name, "\n")
}