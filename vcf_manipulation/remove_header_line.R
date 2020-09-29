#!/usr/bin/env Rscript
.libPaths("/ssd/sda1/sbeatty/software/miniconda3/lib/R/library")
options(verbose=FALSE)
args = commandArgs(trailingOnly=TRUE)
require("stringr", quietly=TRUE)
require("readr", quietly=TRUE)
require("dplyr", quietly=TRUE)

args = commandArgs(trailingOnly=TRUE)

if(length(args) == 0){
	args <- "Rscript /scratch/shahlab_tmp/sbeatty/yvr_pipelines/vcf_manipulation/annotate_vcf.R --target_file=/scratch/shahlab_tmp/sbeatty/IND-15/reheaded_vcf/SA1284AN.shrunk.rehead.vcf --output_file=test.vcf --IDs_to_remove=ANN,chr,start,end,SA_ID,strand,ref_allele,alt,normal_id,tumour_id,mut_pr,tr,ta,nr,na,dbsnp,thousand_genomes,PhenotypeList"
	args <- str_split(args," ") %>% unlist
}

target_vcf_path <- args[str_detect(args,"--target_file")] %>% str_replace("--target_file=","")
output_file_name <- args[str_detect(args,"--output_file")] %>% str_replace("--output_file=","")
IDs_to_remove <- args[str_detect(args,"--IDs_to_remove")] %>% str_replace("--IDs_to_remove=","")



#temporary_header_path <- paste(target_vcf_path,"temporary_header.hdr", sep=".")
temporary_header_path <-paste(str_replace_all(target_vcf_path, '\\.vcf', ".temp_header.hdr"))



# create header file
command <- paste("bcftools view -h", target_vcf_path, ">", temporary_header_path)
system(command)


header <- read_lines(temporary_header_path)
ids_to_remove <- str_split(IDs_to_remove, ",") %>% unlist

for(id in ids_to_remove){
	match_row <- which(str_detect(header, paste0("INFO=<ID=",id)) == TRUE)
	header <- header[-c(match_row)]
}

write_lines(header,temporary_header_path, na = "NA", append = FALSE)

command <- paste("bcftools reheader -h", temporary_header_path, "-o", paste(output_file_name), target_vcf_path)
system(command)

system(paste("rm",paste(temporary_header_path)))

