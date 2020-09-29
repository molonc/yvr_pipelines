#!/usr/bin/env Rscript
.libPaths("/ssd/sda1/sbeatty/software/miniconda3/lib/R/library")
options(verbose=FALSE)
args = commandArgs(trailingOnly=TRUE)
require("stringr", quietly=TRUE)
require("readr", quietly=TRUE)
require("dplyr", quietly=TRUE)




if(length(args) == 0){
	args <- "Rscript /scratch/shahlab_tmp/sbeatty/yvr_pipelines/vcf_manipulation/reorder_annoation_columns.R --annotation_file_in=/scratch/shahlab_tmp/sbeatty/IND-15/shrunk_annotation_file/SA1284T.annotation_file.tsv --annotation_file_out=test_annotation_file_out.tsv"
	args <- str_split(args," ") %>% unlist
}

annotation_file_in <- args[str_detect(args,"--annotation_file_in")] %>% str_replace("--annotation_file_in=","")
annotation_file_out <- args[str_detect(args,"--annotation_file_out")] %>% str_replace("--annotation_file_out=","")

annotation_file <- read_tsv(annotation_file_in, col_types=cols(chr=col_character(), SIFT_effect_prediction=col_character(), SIFT_substitution_type=col_character()))

annotation_file <- annotation_file %>% rename(CHROM=chr, POS=start, END=end)

annotation_file <- annotation_file %>% relocate("END") %>% relocate("POS") %>% relocate("CHROM")

names(annotation_file)[1] <- "#CHROM"
write_tsv(annotation_file, annotation_file_out, col_names=TRUE)