#!/usr/bin/env Rscript
.libPaths("/ssd/sda1/sbeatty/software/miniconda3/lib/R/library")
options(verbose=FALSE)
args = commandArgs(trailingOnly=TRUE)
require("stringr", quietly=TRUE)
require("readr", quietly=TRUE)
require("dplyr", quietly=TRUE)




if(length(args) == 0){
	args <- "Rscript //scratch/shahlab_tmp/sbeatty/IND-15/scripts/subset_annotation.R --annotation_file_in=//scratch/shahlab_tmp/sbeatty/ind231/ind231_summary_V14.csv --annotation_file_out=test_annotation_file_out.tsv --task=1 --target_sample=SA1228N"
	args <- str_split(args," ") %>% unlist
}

annotation_file_in <- args[str_detect(args,"--annotation_file_in")] %>% str_replace("--annotation_file_in=","")
annotation_file_out <- args[str_detect(args,"--annotation_file_out")] %>% str_replace("--annotation_file_out=","")
target_task <- args[str_detect(args,"--task")] %>% str_replace("--task=","")
target_sample <- args[str_detect(args,"--target_sample")] %>% str_replace("--target_sample=","")

annotation_file <- read_csv(annotation_file_in, col_types=cols(chr=col_character(), SIFT_effect_prediction=col_character(), SIFT_substitution_type=col_character()))

annotation_file <- annotation_file %>% filter(SA_ID == target_sample & task == target_task) 


annotation_file$PhenotypeList <- annotation_file$PhenotypeList %>%  str_replace_all(";",",")

write_tsv(annotation_file, paste(annotation_file_out))