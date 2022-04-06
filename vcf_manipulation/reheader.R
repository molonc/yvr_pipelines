#!/usr/bin/env Rscript
#.libPaths("/ssd/sda1/sbeatty/software/miniconda3/lib/R/library")
options(verbose=FALSE)
args = commandArgs(trailingOnly=TRUE)
require("stringr", quietly=TRUE)
require("readr", quietly=TRUE)
require("dplyr", quietly=TRUE)

args = commandArgs(trailingOnly=TRUE)

if(length(args) == 0){
	args <- "Rscript /scratch/shahlab_tmp/sbeatty/yvr_pipelines/vcf_manipulation/annotate_vcf.R --target_file=/scratch/shahlab_tmp/sbeatty/IND-15/shrunk_vcf/SA1228N.shrunk.vcf --output_file=/scratch/shahlab_tmp/sbeatty/ind231/data_export/SA1228test.rehead.vcf --new_header_lines_dataframe=/scratch/shahlab_tmp/sbeatty/ind231/data_export/new_info_tags_dataframe.csv --header_lines_out=test.hdr"
	args <- str_split(args," ") %>% unlist
}


target_vcf_path <- args[str_detect(args,"--target_file")] %>% str_replace("--target_file=","")
output_file_name <- args[str_detect(args,"--output_file")] %>% str_replace("--output_file=","")
header_lines_out <- args[str_detect(args,"--header_lines_out")] %>% str_replace("--header_lines_out=","")

annotation_file_path <- args[str_detect(args,"--annotation_path")] %>% str_replace("--annotation_path=","")


new_header_lines_dataframe <- args[str_detect(args,"--new_header_lines_dataframe")] %>% str_replace("--new_header_lines_dataframe=","")

new_header_lines <- read_csv(new_header_lines_dataframe)
new_header_lines_df <- new_header_lines
temp_header_fix <- "##FILTER=<ID=INDL,Description=“INDEL filter from Kronos”>"
new_header_lines <- mapply(1:nrow(new_header_lines), USE.NAMES=FALSE, FUN=function(x) {unlist(paste0("##INFO=<ID=", new_header_lines$ID[x], ",Number=", new_header_lines$Number[x], ",Type=",new_header_lines$Type[x], ",Description=\"", new_header_lines$Description[x],"\">"))})

#temporary_header_path <- paste(target_vcf_path,"temporary_header.hdr", sep=".")
temporary_header_path <-paste(str_replace_all(target_vcf_path, '\\.vcf', ".temp_header.hdr"))



# create header file
command <- paste("bcftools view -h", target_vcf_path, ">", temporary_header_path)
system(command)


header <- read_lines(temporary_header_path)
last_info_line <- which(str_detect(header, "##INFO") == TRUE) %>% max
header_front <- header[1:last_info_line]
format_fix <- 
header_back <- header[seq(last_info_line+1,length(header))]
header_revised <- c(header_front, new_header_lines, temp_header_fix, header_back)
write_lines(new_header_lines, paste(header_lines_out), na="NA", append = FALSE)
write_lines(header_revised,temporary_header_path, na = "NA", append = FALSE)

#command <- paste("bcftools reheader -h", temporary_header_path, "-o", paste(str_replace_all(target_vcf_path, '\\.vcf', ".reheader.vcf")), target_vcf_path)

command <- paste("bcftools reheader -h", temporary_header_path, "-o", paste(output_file_name), target_vcf_path)

system(command)

system(paste("rm",paste(temporary_header_path)))

