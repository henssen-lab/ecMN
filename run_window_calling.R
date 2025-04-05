library(tidyr)
library(stringr)
library(ggplot2)
library(fitdistrplus)
library(MASS)
library(ggpubr)

args = commandArgs(trailingOnly=TRUE)

meta_table = args[1]
reconstruction = args[2]
output_dir = args[3]



source("./R/data_handler.R")
source("./R/compute_thresholds.R")
source("./R/classify_windows.R")
source("./R/QC_tools.R")
source("./R/analysis_tools.R")


dir.create(paste0(output_dir,"/theshold_plots"))


#0. load meta table and reconstructions
meta <- read.csv(meta_table, sep="\t", header = F, col.names = c("path","sample","class"))
reconstrucion_bed <- read.csv(reconstruction, header = 1, sep = "\t")
#1. read files
bed_file_list <- load_bed_files(meta_table)
#2. compute thresholds
sample_threshold_df <- generate_threshold_df(bed_file_list,P=0.01,percentile=0.99,plot=T,out_dir_plots=paste0(output_dir,"/theshold_plots"))
#3. classify windows
window_calls <- classify_windows(bed_file_list,sample_threshold_df,plot = T, plot_outdir =paste0(output_dir,"/theshold_plots"),window_proximity_n=4)

write.table(window_calls,paste0(output_dir,"/window_calls.txt"), append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE)


