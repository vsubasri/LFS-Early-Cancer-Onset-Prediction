# DNA Methylation Predicts Early Onset of Primary Tumor in Patients with Li-Fraumeni Syndrome

## Project Overview

Li-Fraumeni syndrome (LFS) is an autosomal dominant cancer predisposition syndrome. Approximately 80% of LFS patients harbor a germline TP53 mutation, rendering them susceptible to a wide spectrum of early onset malignancies. A comprehensive surveillance regimen termed the 'Toronto Protocol' has recently been adopted for early tumor detection, demonstrating significant improvement in survival among these patients. However, the protocol's "one-size-fits-all" approach fails to consider an individual patient's risk of cancer.

This project has developed a machine learning model that predicts early onset of primary tumors in LFS patients by estimating the probability of cancer onset before the age of six, leveraging a patient's peripheral blood leukocyte methylation profile.

## Repository Structure

### Preprocessing

```
/preprocessing
```

Contains scripts for preprocessing DNA methylation data, including optimal preprocessing and normalization methods, outlier detection, and confounder removal.

### Model Development

```
/model
```

Contains scripts for feature selection, cross-validation, and model selection. Includes utilities for running models with different seeds and parameters, and combining results.

### Analyses

```
/analyses
```

Scripts for analyzing model results, including SHAP value analysis, covariate comparisons, age-related probe analysis, model comparisons, and family-based analyses.

### Sampling Analysis

```
/sampling
```

Contains scripts related to resampling techniques used to address class imbalance and validate model robustness.

### Immune Cell Deconvolution

```
/idol
```

Scripts for estimating immune cell proportions from methylation data using the IDOL (Identifying Optimal Libraries) method.

### Figures

```
/figures
```

Contains generated figures and visualization scripts for the project.

### Checkpoints

```
/checkpoint
```

Model checkpoints and saved model states.

### Data Directory

```
/data
```

This project uses a centralized data directory structure with the following organization:

- `/data/rds/` - For RDS data files
- `/data/Output/` - For output files
- `/data/Plots/` - For plot files  
- `/data/Resources/` - For resource files like probe lists
- `/data/Data/LFS_450/` - For 450k array data
- `/data/Data/LFS_850/` - For 850k array data

## Setting up the Data Directory

The project has been updated to use relative paths based on a centralized data directory. Run the setup script to create the necessary directory structure:

```bash
./scripts/setup_data_dirs.sh
```

This script will create the required folder structure and provide instructions for migrating existing data.

## Usage

Scripts in this repository are primarily R and shell scripts. Main model execution can be performed using:

```bash
./model/runModel.sh
```

For parallel execution:

```bash
./model/runAllModelinParallel.sh
```

## Citation

If you use this code or model in your research, please cite our work (citation details to be added upon publication).


