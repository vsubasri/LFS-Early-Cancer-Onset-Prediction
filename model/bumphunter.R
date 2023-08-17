source('Scripts/bumphunter_utils.R')
library(bumphunter)
library(dplyr)
library(doParallel)
library(minfi)
library(qvalue)
library(argparse)

setwd('/hpf/largeprojects/davidm/vsubasri/methyl_data/')

## set up parser
parser <- ArgumentParser()
parser$add_argument("--datafile", action="store")

args <- parser$parse_args()
data <- readRDS(args$datafile)

## set up cores
num_cores = detectCores()
registerDoParallel(cores = num_cores-1)

# the desgin matrix =  with rows representing samples and columns representing covariates.
# regression is applied to each row of mat. basically a column for the intercept (= 1) 
# and a column for the variable of interest (cases vs controls).
designMatrix <- data.frame(intcpt = 1,
                           cancer = data$cancerstatus)

ratio_set <- readRDS('Data_objects/raw_ratio_set.rds')
# get granges object

params <- process_params(data,ratio_set)

bumps <- bump_hunter(params, designMatrix, "bumps")

g_ranges <- as.data.frame(getLocations(ratio_set))
# get probes from rownames
g_ranges$probe <- rownames(g_ranges)
# remove ch and duplicatee
g_ranges <- g_ranges[!duplicated(g_ranges$start),]
g_ranges <- g_ranges[!grepl('ch', g_ranges$probe),]

#bumps <- get(load("Data_objects/bumps.rda"))
bumps_table <- bumps$table
bumps_sig <- bumps_table[bumps_table$fwer < 0.05,]
bumps_probes <- inner_join(bumps_sig, g_ranges)$probe


