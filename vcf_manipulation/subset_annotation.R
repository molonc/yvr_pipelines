#!/usr/bin/env Rscript
#.libPaths("/ssd/sda1/sbeatty/software/miniconda3/lib/R/library")
options(verbose=FALSE)
args = commandArgs(trailingOnly=TRUE)
require("stringr", quietly=TRUE)
require("readr", quietly=TRUE)
require("dplyr", quietly=TRUE)


if(length(args) == 0){
	args <- "Rscript /projects/molonc/aparicio_lab/sbeatty/yvr_pipelines/vcf_manipulation/subset_annotation.R --target_file=/projects/molonc/aparicio_lab/sbeatty/archive/IND-231/variant_calls_May2020/ind231_summary_V14.csv --annotation_file_in=/projects/molonc/aparicio_lab/sbeatty/archive/IND-231/variant_calls_May2020/ind231_summary_V14.csv --annotation_file_out=shrunk_annotation_file/SA1259T.SING_STRE.annotation_file.tsv --task=1,2 --target_sample=SA1259T.SING_STRE"
	args <- str_split(args," ") %>% unlist
}



task_file_mapping <- read_csv("/projects/molonc/aparicio_lab/sbeatty/IND/IND-31/sample_task_mapping.csv")
annotation_file_in <- args[str_detect(args,"--annotation_file_in")] %>% str_replace("--annotation_file_in=","")
annotation_file_out <- args[str_detect(args,"--annotation_file_out")] %>% str_replace("--annotation_file_out=","")
target_file <- args[str_detect(args,"--target_file")] %>% str_replace("--target_file=","")

target_sample <- args[str_detect(args,"--target_sample")] %>% str_replace("--target_sample=","")

target_task <- task_file_mapping %>% filter(sample == target_sample) %>% dplyr::select(tasks) %>% unlist %>% first %>% str_replace(";",",")
#target_task <- args[str_detect(args,"--task")] %>% str_replace("--task=","")
target_task <- target_task %>% str_split(",") %>% unlist

target_sample <- str_replace(target_sample, '.PAIR',"")  %>% str_replace(".SING","") %>% str_replace("_MUTA","") %>% str_replace("_STRE","")




annotation_file <- read_csv(annotation_file_in, col_types=cols(chr=col_character(), SIFT_effect_prediction=col_character(), SIFT_substitution_type=col_character()))
dummy_row <- annotation_file[1,]
for(col_n in names(dummy_row)){
	dummy_row[1,paste(col_n)] <- NA
}
dummy_row[1,'chr'] <- "1"
dummy_row[1,"start"] <- 1 
dummy_row[1,'end'] <- 1
dummy_row[1,"strand"] <- "+"
dummy_row[1,"ref"] <- "A"
dummy_row[1,"alt"] <- "T"
dummy_row <- rbind(dummy_row,dummy_row)
annotation_file <- annotation_file[annotation_file$task %in% target_task,]
annotation_file <- annotation_file %>% filter(SA_ID == target_sample)
if(nrow(annotation_file) == 0){
	annotation_file <- dummy_row
}
write_tsv(annotation_file, paste(annotation_file_out))
annotation_file <- annotation_file %>% rename(CHROM=chr, POS=start, END=end)
annotation_file$PhenotypeList <- annotation_file$PhenotypeList %>% str_replace_all(";",",")
write_tsv(annotation_file[,c("CHROM","POS","END")], paste0(annotation_file_out,".regions.tsv"), col_names=FALSE)