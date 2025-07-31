# InferRBPActivityMM

Development of RNA-binding proteins activity inferring model for multiple myeloma risk prediction

## 📋 Overview

**InferRBPActivityMM** is a comprehensive bioinformatics pipeline that infers RNA-binding protein (RBP) activity patterns to predict risk in multiple myeloma (MM) patients. The project combines machine learning approaches with survival analysis to develop robust prognostic models.

### Key Features

- 🧬 **RBP Activity Inference**: Calculate RBP activity scores from RNA splicing data
- 🤖 **Machine Learning Models**: Random Forest Survival models with hyperparameter optimization  
- 📊 **Survival Analysis**: Kaplan-Meier curves and Cox proportional hazards modeling
- 🔍 **External Validation**: Multi-dataset validation (GSE116324, GSE110486)
- 📈 **Comprehensive Visualization**: Professional plots and analysis reports
- ⚙️ **Modular Design**: Clean, well-documented, and extensible codebase

## 🏗️ Project Structure

```
InferRBPActivityMM/
├── README.md                      # Project documentation
├── config/
│   └── config.R                   # Global configuration and parameters
├── src/                           # Source code modules
│   ├── utils/
│   │   └── data_utils.R          # Data processing utilities
│   ├── data_processing/
│   │   └── rbp_activity.R        # RBP activity calculation
│   ├── models/
│   │   └── survival_model.R      # Survival modeling functions
│   ├── analysis/
│   │   └── cofactors_analysis.R  # Cofactors and Cox regression
│   └── visualization/
│       └── plots.R               # Plotting and visualization
├── scripts/
│   └── run_main_analysis.R       # Main analysis pipeline
├── data/                         # Data directories
│   ├── raw/                      # Raw input data
│   ├── processed/                # Processed intermediate data
│   └── results/                  # Analysis results
│       ├── figures/              # Generated plots
│       ├── tables/               # Summary tables
│       └── models/               # Trained models
├── code/                         # Original code (legacy)
│   ├── MMRF Analysis/
│   └── External Validation/
└── logs/                         # Analysis logs
```

## 🚀 Quick Start

### Prerequisites

**R version**: 4.0.2 or higher

**Required R packages**:
```r
# Core packages
install.packages(c("survival", "survminer", "randomForestSRC", 
                   "ggplot2", "gridExtra", "RColorBrewer", "reshape2"))

# Additional packages for extended functionality
install.packages(c("ArrayTools", "rmarkdown"))
```

### Installation

1. **Clone the repository**:
```bash
git clone https://github.com/yourusername/InferRBPActivityMM.git
cd InferRBPActivityMM
```

2. **Download required data files**:
   - MMRF gene expression files (TPM.Rdata, TPM_selected.Rdata) from:
     - http://predscmirna.chengqi.site/RBPData/TPM.Rdata
     - http://predscmirna.chengqi.site/RBPData/TPM_selected.Rdata
   - Place in appropriate data directories

3. **Run the analysis**:
```bash
Rscript scripts/run_main_analysis.R
```

## 📖 Usage

### Basic Analysis

The main analysis pipeline can be executed with a single command:

```r
# In R console
source("scripts/run_main_analysis.R")
```

### Customized Analysis

For more control over the analysis parameters:

```r
# Load configuration
source("config/config.R")

# Modify parameters as needed
DATA_CONFIG$correlation_threshold <- 0.35
MODEL_CONFIG$rf_params$ntree <- 100

# Run specific modules
source("src/data_processing/rbp_activity.R")
# ... continue with analysis
```

### Configuration Options

Key parameters can be modified in `config/config.R`:

```r
# Data parameters
DATA_CONFIG$correlation_threshold <- 0.3  # RBP-target correlation threshold
DATA_CONFIG$train_samples <- 600          # Number of training samples

# Model parameters  
MODEL_CONFIG$rf_params$ntree <- 50        # Number of trees in Random Forest
MODEL_CONFIG$rf_params$mtry <- 15         # Number of variables per split

# Analysis parameters
ANALYSIS_CONFIG$alpha_level <- 0.05       # Statistical significance level
```

## 🔬 Methodology

### 1. Data Preprocessing
- **PSI Normalization**: Range normalization and median centering
- **Quality Control**: IQR-based filtering of low-variance features
- **Data Splitting**: Reproducible train/validation splits

### 2. RBP Activity Inference
- **CLIP Target Identification**: Integration of ENCODE CLIP-seq data
- **Correlation Analysis**: Spearman correlation between RBP expression and target splicing
- **Activity Scoring**: Weighted combination of positive/negative target effects

### 3. Survival Modeling
- **Random Forest Survival**: Non-parametric ensemble method
- **Hyperparameter Optimization**: 5-fold cross-validation grid search
- **Risk Stratification**: Median-based high/low risk classification

### 4. Statistical Analysis
- **Kaplan-Meier Curves**: Survival probability estimation
- **Log-rank Tests**: Group comparison statistics
- **Cox Regression**: Multivariate hazard modeling with cofactors

## 📊 Results

### Primary Analysis (MMRF Dataset)
- **Training Set**: 600 multiple myeloma patients
- **Validation Set**: 162 patients  
- **Features**: 150 RBP activity scores
- **Endpoint**: Overall survival

### External Validation
- **GSE116324**: 44 MM samples with ISS staging
- **GSE110486**: 17 MM samples vs normal controls

### Key Findings
- RBP activity signatures effectively stratify MM patients into risk groups
- Model shows significant survival differences (p < 0.05)
- External validation confirms robustness across datasets
- ISS stage correlation validates clinical relevance

## 📁 Output Files

The analysis generates several types of output:

### Figures (`data/results/figures/`)
- `*_survival_curves.png`: Kaplan-Meier survival plots
- `*_risk_distribution.png`: Risk score distributions  
- `*_ISS_stage_comparison.png`: ISS stage analysis plots
- `*_target_numbers.png`: RBP target number visualizations

### Tables (`data/results/tables/`)
- `validation_risk_scores.csv`: Individual patient risk scores
- `cox_regression_summary.csv`: Cofactors analysis results
- `model_comparison.csv`: Model performance metrics

### Models (`data/results/models/`)
- `main_analysis_results.Rdata`: Complete analysis objects
- `*_results.Rdata`: External validation results

## 🔧 Advanced Usage

### Custom Data Analysis

To analyze your own dataset:

```r
# Load your data
your_psi_matrix <- load_your_psi_data()
your_rbp_matrix <- load_your_rbp_data()
your_survival_data <- load_your_survival_data()

# Run RBP activity calculation
rbp_results <- process_rbp_activity_pipeline(
  data_split = your_data_split,
  PSI_matrix = your_psi_matrix,
  RBP_matrix = your_rbp_matrix,
  correlation_matrix = your_correlation_matrix,
  transcript_names = your_transcript_names
)

# Continue with survival analysis...
```

### Extending the Pipeline

Add new analysis modules:

```r
# Create new analysis function
your_analysis_function <- function(data, params) {
  # Your analysis code here
  return(results)
}

# Integrate into main pipeline
source("src/your_analysis/your_module.R")
```

## 🐛 Troubleshooting

### Common Issues

1. **Missing dependencies**:
   ```r
   # Check and install missing packages
   source("config/config.R")
   validate_config()
   ```

2. **Memory issues with large datasets**:
   ```r
   # Reduce dataset size or increase memory
   options(java.parameters = "-Xmx8g")
   ```

3. **File path issues**:
   ```r
   # Check working directory
   getwd()
   
   # Verify file paths in config
   source("config/config.R")
   ```

### Error Messages

| Error | Solution |
|-------|----------|
| "File not found" | Check data file paths in config |
| "Package not installed" | Install missing R packages |
| "Insufficient memory" | Reduce data size or increase RAM |

## 🤝 Contributing

We welcome contributions! Please follow these steps:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Development Guidelines

- Follow R coding standards
- Add comprehensive documentation
- Include unit tests for new functions
- Update README for new features

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📚 Citation

If you use this software in your research, please cite:

```bibtex
@software{inferrbpactivitymm,
  title={InferRBPActivityMM: RNA-binding protein activity inference for multiple myeloma risk prediction},
  author={Your Name},
  year={2024},
  url={https://github.com/yourusername/InferRBPActivityMM}
}
```

## 📞 Support

- **Issues**: [GitHub Issues](https://github.com/yourusername/InferRBPActivityMM/issues)
- **Documentation**: This README and inline code documentation
- **Contact**: your.email@institution.edu

## 🙏 Acknowledgments

- MMRF (Multiple Myeloma Research Foundation) for providing the dataset
- TCGA consortium for genomic data
- GEO database for external validation datasets
- R community for excellent statistical computing tools

---

*Last updated: December 2024*
