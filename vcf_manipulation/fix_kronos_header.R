#!/usr/bin/env Rscript

options(verbose=FALSE)
args = commandArgs(trailingOnly=TRUE)
require("stringr", quietly=TRUE)
require("readr", quietly=TRUE)
require("dplyr", quietly=TRUE)

args = commandArgs(trailingOnly=TRUE)

if(length(args) == 0){
	args <- "Rscript /projects/molonc/aparicio_lab/sbeatty/yvr_pipelines/vcf_manipulation/fix_kronos_header.R --target_file=/projects/molonc/aparicio_lab/sbeatty/archive/BXE/kronos_task6/SA1139T.PAIR_MUTA.TASK_6_COSMIC_MutationSeq.annotSnpEff.annotMA.flagDBsnp.flag1000gen.flagCosmic.vcf --output_file=fixed_header/SA1139T.PAIR_MUTA.fixed_header.vcf"
	args <- str_split(args," ") %>% unlist
}


target_vcf_path <- args[str_detect(args,"--target_file")] %>% str_replace("--target_file=","")
output_file_name <- args[str_detect(args,"--output_file")] %>% str_replace("--output_file=","")
add_contig_lines <- args[str_detect(args,"--add_contig_lines")] %>% str_replace("--add_contig_lines=","")

temporary_header_path <- paste0(target_vcf_path,".temp_header.hdr")

print(target_vcf_path)
print(output_file_name)

# create header file
command <- paste("bcftools view -h", target_vcf_path, ">", temporary_header_path)
system(command)

fields_to_fix <- list(
	c("Phred-scaled genotype likelihoods>","Phred-scaled genotype likelihoods"),
	c("#INFO=<ID=GT,Number=1","#FORMAT=<ID=GT,Number=G"),
	#c("#INFO=<ID=PL","#FORMAT=<ID=PL,Number=G"),
	c("#INFO=<ID=DP","#FORMAT=<ID=DP"),
	c("#INFO=<ID=FT","#FORMAT=<ID=FT"),
	c("#INFO=<ID=GL,","#FORMAT=<ID=GL,"),
	c("#INFO=<ID=GLE","#FORMAT=<ID=GLE"),
	c("#INFO=<ID=GP","#FORMAT=<ID=GP"),
	c("#INFO=<ID=GQ","#FORMAT=<ID=GQ"),
	c("#INFO=<ID=HQ","#FORMAT=<ID=GP"),
	c("#INFO=<ID=GP","#FORMAT=<ID=HQ"),
	c("#INFO=<ID=PS","#FORMAT=<ID=PS"),
	c("#INFO=<ID=PQ","#FORMAT=<ID=PQ"),
	c("#INFO=<ID=EC","#FORMAT=<ID=EC"),
	c("#INFO=<ID=MQ","#FORMAT=<ID=MQ"),
	c("##FORMAT=<ID=PL,Number=G,Number=.,Type=String,Description=\"Missing description: this INFO line was added by Picard's FixVCFHeader\">", "##FORMAT=<ID=PL,Number=G,Number=.,Type=String,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">"),
	c("##INFO=<ID=GT,Number=.,Type=String,Description=\"Missing description: this INFO line was added by Picard's FixVCFHeader\">",
"##INFO=<ID=GT,Number=.,Type=String,Description=\"Genotype\">"),
	c("##FILTER=<ID=INDL,Description=\"Missing description: this FILTER line was added by Picard's FixVCFHeader\">","##FILTER=<ID=INDL,Description=\"Kronos INDEL filter\">")
	)


kronos_contig_header_lines <- read_lines("/projects/molonc/aparicio_lab/sbeatty/IND/IND-31/contig_header_lines.tsv")



header <- read_lines(temporary_header_path)
for(field in fields_to_fix){
	header <-  header %>% str_replace(field[1],field[2])
}


if(str_detect(target_vcf_path, "_MUTA") & add_contig_lines == "TRUE"){
	if(length(which(str_detect(header, "##contig") == TRUE)) == 0){
		header_front <- header[1:which(str_detect(header,"##reference") == TRUE)]
		header_back <- header[-c(1:which(str_detect(header,"##reference") == TRUE))]
		header <- c(header_front,kronos_contig_header_lines,header_back)
	}
}

#header <-  header %>% str_replace("Phred-scaled genotype likelihoods>","Phred-scaled genotype likelihoods")

#header <-  header %>% str_replace("#INFO=<ID=GT,Number=1","#FORMAT=<ID=GT,Number=G")
#header <-  header %>% str_replace("#INFO=<ID=PL","#FORMAT=<ID=PL,Number=3")


write_lines(header,temporary_header_path, na = "NA", append = FALSE)

#command <- paste("bcftools reheader -h", temporary_header_path, "-o", paste(str_replace_all(target_vcf_path, '\\.vcf', ".reheader.vcf")), target_vcf_path)

command <- paste("bcftools reheader -h", temporary_header_path, "-o", output_file_name, target_vcf_path)

#standard_tags <- header %>% last %>% str_split("\t") %>% unlist %>% str_replace("#","") %>% paste("-F",., collapse=" ")
system(command)

system(paste("rm",paste(temporary_header_path)))

