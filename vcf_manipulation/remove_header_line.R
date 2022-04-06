#!/usr/bin/env Rscript
#.libPaths("/ssd/sda1/sbeatty/software/miniconda3/lib/R/library")
options(verbose=FALSE)
args = commandArgs(trailingOnly=TRUE)
require("stringr", quietly=TRUE)
require("readr", quietly=TRUE)
require("dplyr", quietly=TRUE)

args = commandArgs(trailingOnly=TRUE)

if(length(args) == 0){
	args <- "Rscript /projects/molonc/aparicio_lab/sbeatty/yvr_pipelines/vcf_manipulation/remove_header_line.R --target_file=/projects/molonc/aparicio_lab/sbeatty/IND/IND-31/removed_tags/SA1266T.PAIR_MUTA.shrunk.rehead.annotated.removed_tags.vcf --output_file=test.vcf"
	args <- str_split(args," ") %>% unlist
}

target_vcf_path <- args[str_detect(args,"--target_file")] %>% str_replace("--target_file=","")
output_file_name <- args[str_detect(args,"--output_file")] %>% str_replace("--output_file=","")
#IDs_to_remove <- args[str_detect(args,"--IDs_to_remove")] %>% str_replace("--IDs_to_remove=","")



#temporary_header_path <- paste(target_vcf_path,"temporary_header.hdr", sep=".")
temporary_header_path <-paste(str_replace_all(target_vcf_path, '\\.vcf', ".temp_header.hdr"))



# create header file
command <- paste("bcftools view -h", target_vcf_path, ">", temporary_header_path)
system(command)


header <- read_lines(temporary_header_path)
#ids_to_remove <- str_split(IDs_to_remove, ",") %>% unlist

#for(id in ids_to_remove){
#	match_row <- which(str_detect(header, paste0("INFO=<ID=",id)) == TRUE)
#	header <- header[-c(match_row)]
#}

write_lines(header,temporary_header_path, na = "NA", append = FALSE)

command <- paste("bcftools reheader -h", temporary_header_path, "-o", paste(output_file_name), target_vcf_path)
system(command)

system(paste("rm",paste(temporary_header_path)))

