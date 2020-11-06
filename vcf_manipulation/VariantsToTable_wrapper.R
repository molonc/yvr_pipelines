#!/usr/bin/env Rscript

options(verbose=FALSE)
args = commandArgs(trailingOnly=TRUE)
require("stringr", quietly=TRUE)
require("readr", quietly=TRUE)
require("dplyr", quietly=TRUE)

args = commandArgs(trailingOnly=TRUE)

if(length(args) == 0){
	args <- "Rscript /projects/molonc/aparicio_lab/sbeatty/yvr_pipelines/vcf_manipulation/VariantsToTable_wrapper.R --target_file=input_vcfs/SA1284T.SING_STRE.vcf.gz --output_file=test.tsv -—gatk_command_path=/projects/molonc/aparicio_lab/sbeatty/software/gatk-4.1.8.1/gatk"
	args <- str_split(args," ") %>% unlist
}

target_file_path <- args[str_detect(args,"--target_file")] %>% str_replace("--target_file=","")
output_file_path <- args[str_detect(args,"--output_file")] %>% str_replace("--output_file=","")
gatk_command_path <- args[str_detect(args,"-—gatk_command_path")] %>% str_replace("-—gatk_command_path=","")


temp_header_path <- target_file_path %>% str_split("/") %>% unlist %>% last %>% paste0("./",.,".temp.hdr")

command <- paste("bcftools view -h",target_file_path, ">",temp_header_path)

print(command)

system(command)

header <- readLines(temp_header_path)
info_lines <- header[str_detect(header,"##INFO=<")]

parse_info_line <- function(info_line_string){
	output <- info_line_string %>% str_split(",") %>% unlist %>% first %>% str_split("=") %>% unlist %>% last
	output
}



info_tags <- mapply(info_lines, USE.NAMES=FALSE, FUN=parse_info_line)
standard_tags <- header %>% last %>% str_split("\t") %>% unlist %>% str_replace("#","")
GATK_tags <- c("EVENTLENGTH", "TRANSITION", "HET", "HOM-REF", "HOM-VAR", "NO-CALL", "TYPE", "VAR", "NSAMPLES", "NCALLED", "MULTI-ALLELIC")
tags <- c(standard_tags, info_tags, GATK_tags)
tags <- tags  %>% paste("-F",.) %>% paste(collapse=" ")

GATK_command_call <- paste("/projects/molonc/aparicio_lab/sbeatty/software/gatk-4.1.8.1/gatk VariantsToTable -V", target_file_path,"--show-filtered", tags,"-O", output_file_path)
print(GATK_command_call)
system(paste(GATK_command_call))
system(paste("/projects/molonc/aparicio_lab/sbeatty/software/gatk-4.1.8.1/gatk VariantsToTable -V", target_file_path,"--show-filtered", tags,"-O test_vanilla.tsv"))

#system("/projects/molonc/aparicio_lab/sbeatty/software/gatk-4.1.8.1/gatk VariantsToTable -V /projects/molonc/aparicio_lab/sbeatty/IND/IND-31/annotate_task6/filter_low_qual/SA1284T.SING_STRE.fixed_header.remove.low.qual.vcf --show-filtered -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F INFO -F FORMAT -F SA1284T -F END -F BLOCKAVG_min30p3a -F SNVHPOL -F CIGAR -F RU -F REFREP -F IDREP -F MQ -F EVENTLENGTH -F TRANSITION -F HET -F HOM-REF -F HOM-VAR -F NO-CALL -F TYPE -F VAR -F NSAMPLES -F NCALLED -F MULTI-ALLELIC -O vanilla_test.tsv")

#system("/projects/molonc/aparicio_lab/sbeatty/software/gatk-4.1.8.1/gatk VariantsToTable -V /projects/molonc/aparicio_lab/sbeatty/IND/IND-31/annotate_task6/filter_low_qual/SA1284T.SING_STRE.fixed_header.remove.low.qual.vcf -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F INFO -F FORMAT -F SA1284T -F END -F BLOCKAVG_min30p3a -F SNVHPOL -F CIGAR -F RU -F REFREP -F IDREP -GF MQ -F EVENTLENGTH -F TRANSITION -F HET -F HOM-REF -F HOM-VAR -F NO-CALL -F TYPE -F VAR -F NSAMPLES -F NCALLED -F MULTI-ALLELIC --show-filtered -O test_move_show_filetered.tsv")

#/projects/molonc/aparicio_lab/sbeatty/software/gatk-4.1.8.1/gatk VariantsToTable -V /projects/molonc/aparicio_lab/sbeatty/IND/IND-31/annotate_task6/filter_low_qual/SA1284T.SING_STRE.fixed_header.remove.low.qual.vcf -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F INFO -F FORMAT -F SA1284T -F END -F BLOCKAVG_min30p3a -F SNVHPOL -F CIGAR -F RU -F REFREP -F IDREP -GF MQ -F EVENTLENGTH -F TRANSITION -F HET -F HOM-REF -F HOM-VAR -F NO-CALL -F TYPE -F VAR -F NSAMPLES -F NCALLED -F MULTI-ALLELIC --show-filtered -O test_csv.csv