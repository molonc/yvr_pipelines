#!/usr/bin/env Rscript
.libPaths("/ssd/sda1/sbeatty/software/miniconda3/lib/R/library")
options(verbose=FALSE)
args = commandArgs(trailingOnly=TRUE)
require("stringr", quietly=TRUE)
require("readr", quietly=TRUE)
require("dplyr", quietly=TRUE)

args = commandArgs(trailingOnly=TRUE)

#target_vcf_path <- "/scratch/shahlab_tmp/sbeatty/ind231/data_export/vcf_subset.vcf"
#output_file_name <- "vcf_output.vcf"
#temporary_header_path <- paste(target_vcf_path,"temporary_header.hdr", sep=".")
##line_to_add <- "##INFO=<ID=FATHMM,Number=.,Type=String,Description=\"FATHMM score\">"
#new_header_line_ID <- "FATHMM,Number"
#new_header_line_description <- "FATHMM score"
#new_annotation_column <- "FATHMM"
/scratch/shahlab_tmp/danlai/APARICIO-590/SING_MUTA/mutationseq_pipeline/OUTPUT/RUN/SA1228N_museq/outputs/results/TASK_6_COSMIC_SA1228N_MutationSeq.annotSnpEff.annotMA.flagDBsnp.flag1000gen.flagCosmic.vcf

if(length(args) == 0){
	args <- "Rscript /scratch/shahlab_tmp/sbeatty/yvr_pipelines/vcf_manipulation/annotate_vcf.R --target_file=TASK_6_STRELKA_FLAG_COSMIC_SA1228T_strelka.passed.somatic.snvs.annotSnpEff.annotMA.flagDBsnp.flag1000gen.flagCosmic.vcf --output_file=/scratch/shahlab_tmp/sbeatty/ind231/data_export/SA1228test.rehead.vcf --new_header_lines_dataframe=/scratch/shahlab_tmp/sbeatty/IND-15/metadata/new_header_lines.csv --annotation_path=/scratch/shahlab_tmp/sbeatty/ind231/ind231_summary_V14.csv --target_sample=SA1228T --task=2"
	args <- str_split(args," ") %>% unlist
}



target_vcf_path <- args[str_detect(args,"--target_file")] %>% str_replace("--target_file=","")
target_task <- args[str_detect(args,"--task")] %>% str_replace("--task=","")
target_task <- target_task %>% str_split(",") %>% unlist 
target_sample <- args[str_detect(args,"--target_sample")] %>% str_replace("--target_sample=","")
target_vcf_path <- args[str_detect(args,"--target_file")] %>% str_replace("--target_file=","")
output_file_name <- args[str_detect(args,"--output_file")] %>% str_replace("--output_file=","")
new_header_lines_dataframe <- args[str_detect(args,"--new_header_lines_dataframe")] %>% str_replace("--new_header_lines_dataframe=","")
annotation_file_path <- args[str_detect(args,"--annotation_path")] %>% str_replace("--annotation_path=","")

new_header_lines <- read_csv("new_info_tags_dataframe.csv")
new_header_lines_df <- new_header_lines
new_header_lines <- mapply(1:nrow(new_header_lines), USE.NAMES=FALSE, FUN=function(x) {unlist(paste0("##INFO=<ID=", new_header_lines$ID[x], ",Number=", new_header_lines$Number[x], ",Type=",new_header_lines$Type[x], ",Description=\'", new_header_lines$Description[x],"\'>"))})

temporary_header_path <- paste(target_vcf_path,"temporary_header.hdr", sep=".")


# create header file
command <- paste("bcftools view -h", target_vcf_path, ">", temporary_header_path)
system(command)



header <- read_lines(temporary_header_path)
last_info_line <- which(str_detect(header, "##INFO") == TRUE) %>% max
header_front <- header[1:last_info_line]
header_back <- header[seq(last_info_line+1,length(header))]
header_revised <- c(header_front, new_header_lines, header_back)
write_lines(header_revised,temporary_header_path, na = "NA", append = FALSE)

command <- paste("bcftools reheader -h", temporary_header_path, "-o", paste0(target_vcf_path,".reheader.vcf"), target_vcf_path)
system(command)


# shrink the annoation file 
annotation_file <- read_csv(annotation_file_path, col_types=cols(chr=col_character(), SIFT_effect_prediction=col_character(), SIFT_substitution_type=col_character()))
annotation_file <- annotation_file[annotation_file$task %in% target_task,]
annotation_file <- annotation_file %>% filter(SA_ID == target_sample)
annotation_file <- annotation_file %>% rename(ref_allele=ref)
annotation_file <- annotation_file[,c(new_header_lines_df$ID)]
annotation_file <- annotation_file %>% rename(CHROM = chr, FROM = start, TO=end)
annotation_file <- annotation_file %>% relocate(TO) %>% relocate(FROM) %>% relocate(CHROM)   
annotation_file_output <- annotation_file %>% as.data.frame

write_tsv(data.frame(annotation_file_output), path=paste(output_file_name,"annotation_reference.tsv", sep="."), col_names=FALSE)

system(paste("rm", paste(output_file_name,"annotation_reference.tsv.gz", sep=".")))

command <- paste("bgzip", paste(output_file_name,"annotation_reference.tsv", sep="."))

system(command)
s <- which(names(annotation_file) == "CHROM")
b <- which(names(annotation_file) == "FROM")
e <- which(names(annotation_file) == "TO")


command <- paste0("tabix ",paste(output_file_name,"annotation_reference.tsv", "gz",sep="."), " -s", s, " -b", b, " -e", e)
system(command)
new_annoation_columns <-  new_header_lines_df$ID %>% str_replace("chr", "CHROM") %>% str_replace("start", 'FROM') %>% str_replace("end", 'TO')
new_annoation_columns <- c("CHROM",'FROM','TO',new_annoation_columns[!c(new_annoation_columns %in% c("CHROM","FROM","TO"))])

command <- paste("bcftools annotate -a", paste(output_file_name,"annotation_reference.tsv", "gz",sep="."), "-c",paste( new_annoation_columns, collapse=","),paste0(target_vcf_path,".reheader.vcf"),">", output_file_name)

write_tsv(annotation_file[1:3,], "regions.tsv", col_names=FALSE)

output_file_name

system(command)
a <- readVcf("small.vcf")

system(paste("rm", temporary_header_path))
system(paste0("rm ",output_file_name,".reheader.vcf"))

