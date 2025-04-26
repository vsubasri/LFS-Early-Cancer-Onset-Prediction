# DNA Methylation Predicts Early Onset of Primary Tumor in Patients with Li-Fraumeni Syndrome

## Project Overview

Li-Fraumeni syndrome (LFS) is an autosomal dominant cancer predisposition syndrome. Approximately 80% of LFS patients harbor a germline TP53 mutation, rendering them susceptible to a wide spectrum of early onset malignancies. A comprehensive surveillance regimen termed the 'Toronto Protocol' has recently been adopted for early tumor detection, demonstrating significant improvement in survival among these patients. However, the protocol's "one-size-fits-all" approach fails to consider an individual patient's risk of cancer.

This project has developed a machine learning model that predicts early onset of primary tumors in LFS patients by estimating the probability of cancer onset before the age of six, leveraging a patient's peripheral blood leukocyte methylation profile.

## Repository Structure

```
/
├── src/                   # Source code
│   ├── preprocessing/     # Data preprocessing scripts
│   ├── model/             # Model training and evaluation
│   ├── analyses/          # Analysis scripts
│   ├── utils/             # Utility functions
│   └── scripts/           # Shell scripts for running pipelines
├── data/                  # Data directory (not committed)
│   ├── raw/               # Raw data
│   ├── processed/         # Processed data
│   └── results/           # Model results
├── figures/               # Generated figures and plots
├── checkpoint/            # Model checkpoints
├── idol/                  # Immune cell deconvolution
└── sampling/              # Sampling analysis scripts
```

### Preprocessing (src/preprocessing)

Contains scripts for preprocessing DNA methylation data, including optimal preprocessing and normalization methods, outlier detection, and confounder removal.

### Model Development (src/model)

Contains scripts for feature selection, cross-validation, and model selection. Includes utilities for running models with different seeds and parameters, and combining results.

### Analyses (src/analyses)

Scripts for analyzing model results, including SHAP value analysis, covariate comparisons, age-related probe analysis, model comparisons, and family-based analyses.

### Utility Functions (src/utils)

Common utility functions used across different parts of the codebase.

### Shell Scripts (src/scripts)

Shell scripts for running complete pipelines and workflows.

### Data (data/)

Data directory (not committed to version control) containing raw data, processed data, and model results.

### Figures (figures/)

Contains generated figures and visualization scripts for the project.

### Checkpoints (checkpoint/)

Model checkpoints and saved model states.

### Immune Cell Deconvolution (idol/)

Scripts for estimating immune cell proportions from methylation data using the IDOL (Identifying Optimal Libraries) method.

### Sampling Analysis (sampling/)

Contains scripts related to resampling techniques used to address class imbalance and validate model robustness.

## Usage

Scripts in this repository are primarily R and shell scripts. Main model execution can be performed using:

```bash
./src/scripts/run_model.sh
```

For parallel execution:

```bash
./src/scripts/run_parallel.sh
```

## Citation

If you use this code or model in your research, please cite our work (citation details to be added upon publication).


