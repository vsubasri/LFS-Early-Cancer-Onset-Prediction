library(dplyr)
library(stringr)
library(ggplot2)
library(ROCR)
library(data.table)
library(caret)
library(gridExtra)
library(grid)
library(argparse)

# Use relative paths
script_dir <- dirname(sys.frame(1)$ofile)
source(file.path(script_dir, 'modelUtils.R'))

parser <- ArgumentParser()
parser$add_argument("--dir", action="store")
args <- parser$parse_args()
dir <- args$dir

get_avg_results(dir)
