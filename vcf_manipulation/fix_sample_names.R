#!/usr/bin/env Rscript
#.libPaths("/ssd/sda1/sbeatty/software/miniconda3/lib/R/library")
options(verbose=FALSE)
args = commandArgs(trailingOnly=TRUE)
require("stringr", quietly=TRUE)
require("readr", quietly=TRUE)
require("dplyr", quietly=TRUE)

args = commandArgs(trailingOnly=TRUE)

if(length(args) == 0){
	args <- "Rscript /projects/molonc/aparicio_lab/sbeatty/yvr_pipelines/vcf_manipulation/fix_sample_names.R --target_file=/projects/molonc/aparicio_lab/sbeatty/IND/IND-31/removed_tags/SA1266T.PAIR_MUTA.shrunk.rehead.annotated.removed_tags.vcf --output_file=test.vcf --sample=SA1266T.PAIR_MUTA"
	args <- str_split(args," ") %>% unlist
}

target_vcf_path <- args[str_detect(args,"--target_file")] %>% str_replace("--target_file=","")
output_file_name <- args[str_detect(args,"--output_file")] %>% str_replace("--output_file=","")
sample_in <- args[str_detect(args,"--sample")] %>% str_replace("--sample=","")
sample_in <- sample_in %>% str_split("\\.") %>% unlist %>% first
temp_header_path <- target_vcf_path %>% str_replace("\\.vcf", ".temp.hdr")

sample_sheet <- read_csv('//projects/molonc/aparicio_lab/sbeatty/IND/germline_snv_sample_sheet.csv', col_types=cols(TB_ID=col_character()))



command <- paste("bcftools view -h", target_vcf_path, "> ", temp_header_path)
system(command)
header <- read_lines(temp_header_path)

tb_id <- filter(sample_sheet, sample_name == sample_in) %>% dplyr::select(TB_ID) %>% unlist
header <- str_replace_all(header,sample_in, tb_id) %>% str_replace_all("##FILTER=<ID=INDL,Description=\"Missing description: this FILTER line was added by Picard's FixVCFHeader\">","##FILTER=<ID=INDL,Description=\"Kronos INDEL filter\">")

header <- header[!str_detect(header,"bcftools")]

write_lines(header,temp_header_path, na = "NA", append = FALSE)


#command <- paste("bcftools reheader -h", temp_header_path, "-o", output_file_name, target_vcf_path)

command <- paste("bcftools reheader -h", temp_header_path,target_vcf_path, "> ", output_file_name )

system(paste(command))
system(paste("rm", temp_header_path))