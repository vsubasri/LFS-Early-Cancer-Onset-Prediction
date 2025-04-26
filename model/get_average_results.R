library(dplyr)
library(stringr)
library(ggplot2)
library(ROCR)
library(data.table)
library(caret)
library(gridExtra)
library(grid)
library(argparse)

source('modelUtils.R')

parser <- ArgumentParser()
parser$add_argument("--dir", action="store")
args <- parser$parse_args()
dir <- args$dir

get_avg_results(dir)
