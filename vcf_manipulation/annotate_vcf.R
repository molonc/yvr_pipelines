#!/usr/bin/env Rscript
.libPaths("/ssd/sda1/sbeatty/software/miniconda3/lib/R/library")
options(verbose=FALSE)
args = commandArgs(trailingOnly=TRUE)

require("stringr", quietly=TRUE)
require("readr", quietly=TRUE)
require("dplyr", quietly=TRUE)
#require("parallel", quietly=TRUE)
#require("VariantAnnotation", quietly=TRUE)
#require("data.table", quietly=TRUE)
#options(echo=TRUE)

args = commandArgs(trailingOnly=TRUE)

#target_vcf_path <- "/scratch/shahlab_tmp/sbeatty/ind231/data_export/vcf_subset.vcf"
#output_file_name <- "vcf_output.vcf"
#temporary_header_path <- paste(target_vcf_path,"temporary_header.hdr", sep=".")
##line_to_add <- "##INFO=<ID=FATHMM,Number=.,Type=String,Description=\"FATHMM score\">"
#new_header_line_ID <- "FATHMM,Number"
#new_header_line_description <- "FATHMM score"
#new_annotation_column <- "FATHMM"


if(length(args) == 0){
args <- "Rscript /scratch/shahlab_tmp/sbeatty/yvr_pipelines/vcf_manipulation/annotate_vcf.R --target_file=/scratch/shahlab_tmp/sbeatty/ind231/data_export/vcf_subset.vcf --output_file=/scratch/shahlab_tmp/sbeatty/ind231/data_export/vcf_output.vcf --new_header_line_ID=FATHMM,Number --new_header_line_description=FATHMM_score --new_annotation_column=FATHMM"
args <- str_split(args," ") %>% unlist
}


target_vcf_path <- args[str_detect(args,"--target_file")] %>% str_replace("--target_file=","")
output_file_name <- args[str_detect(args,"--output_file")] %>% str_replace("--output_file=","")
new_header_line_info <- args[str_detect(args,"--new_header_line_info")] %>% str_replace("--new_header_line_info=","")
new_annotation_column <- args[str_detect(args,"--outputfile")] %>% str_replace("--outputfile=","")


new_header_line_ID <- args[str_detect(args,"--new_header_line_ID")] %>% str_replace("--new_header_line_ID=","")
new_header_line_description <- args[str_detect(args,"--new_header_line_description")] %>% str_replace("--new_header_line_description=","") %>% str_replace("_", " ")

line_to_add <- paste0("##INFO=<ID=", new_header_line_ID, "=.,Type=String,Description=",'"', new_header_line_description, '"')

annotation_file_path <- args[str_detect(args,"--annotation_path")] %>% str_replace("--annotation_path=","")
new_annotation_column <- args[str_detect(args,"--outputfile")] %>% str_replace("--outputfile=","")

temporary_header_path <- paste(target_vcf_path,"temporary_header.hdr", sep=".")

#add header 
new_annotation_column <- "FATHMM"
annotation_file_path <- "annots.tsv.gz"
#final_vcf_path <- "newest_vcf.vcf"

#annote
#command <- paste("bcftools annotate -a", annotation_file_path, "-h", new_header_path, "-c", paste("CHROM,FROM,TO", new_annotation_column, sep=","), vcf_with_new_header_path,">", final_vcf_path)
#system(command)


# create header file
command <- paste("bcftools view -h", target_vcf_path, ">", temporary_header_path)
system(command)


header <- read_lines(temporary_header_path)
last_info_line <- which(str_detect(header, "##INFO") == TRUE) %>% max
header_front <- header[1:last_info_line]
header_back <- header[seq(last_info_line+1,length(header))]
header_revised <- c(header_front, line_to_add, header_back)
write_lines(header_revised,temporary_header_path, na = "NA", append = FALSE)

command <- paste("bcftools reheader -h", temporary_header_path, "-o", paste0(output_file_name,".reheader.vcf"), target_vcf_path)
system(command)

command <- paste("bcftools annotate -a", annotation_file_path, "-c", paste("CHROM,FROM,TO", new_annotation_column, sep=","), paste0(output_file_name,".reheader.vcf"),">", output_file_name)
system(command)
system(paste("rm", temporary_header_path))
system(paste0("rm ",output_file_name,".reheader.vcf"))

