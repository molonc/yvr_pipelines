#!/usr/bin/env Rscript

options(verbose=FALSE)
args = commandArgs(trailingOnly=TRUE)
require("stringr", quietly=TRUE)
require("readr", quietly=TRUE)
require("dplyr", quietly=TRUE)

args = commandArgs(trailingOnly=TRUE)

annotation_col_specs <- cols(
  `#CHROM` = col_character(),
  FROM = col_double(),
  TO = col_double(),
  ID = col_character(),
  ref = col_character(),
  alt = col_character(),
  QUAL = col_double(),
  FILTER = col_character(),
  INFO = col_character(),
  PR = col_double(),
  TC = col_character(),
  TR = col_character(),
  TA = col_character(),
  NR = col_character(),
  V14 = col_character(),
  ND = col_character(),
  NI = col_character(),
  Allele = col_character(),
  Annotation = col_character(),
  Annotation_Impact = col_character(),
  Gene_Name = col_character(),
  Gene_ID = col_character(),
  Feature_Type = col_character(),
  Feature_ID = col_character(),
  Transcript_BioType = col_character(),
  Rank = col_character(),
  HGVS.c = col_character(),
  HGVS.p = col_character(),
  `cDNA.pos / cDNA.length` = col_character(),
  `CDS.pos / CDS.length` = col_character(),
  `AA.pos / AA.length` = col_character(),
  Distance = col_character(),
  `ERRORS / WARNINGS / INFO` = col_character(),
  LOF = col_character(),
  NMD = col_character(),
  MA = col_character(),
  DBSNP = col_character(),
  X1000Gen = col_character(),
  Cosmic = col_character(),
  EVENTLENGTH = col_character(),
  TRANSITION = col_character(),
  HET = col_character(),
  HOM.REF = col_character(),
  HOM.VAR = col_character(),
  NO.CALL = col_character(),
  TYPE = col_character(),
  VAR = col_character(),
  NSAMPLES = col_character(),
  NCALLED = col_character(),
  MULTI.ALLELIC = col_character(),
  V35 = col_character(),
  genecode_gene_name = col_character(),
  genecode_strand = col_character(),
  unique_id = col_character(),
  X.AlleleID = col_character(),
  Type = col_character(),
  Name = col_character(),
  GeneID = col_character(),
  GeneSymbol = col_character(),
  HGNC_ID = col_character(),
  ClinicalSignificance = col_character(),
  ClinSigSimple = col_character(),
  LastEvaluated = col_character(),
  RS...dbSNP. = col_character(),
  nsv.esv..dbVar. = col_character(),
  RCVaccession = col_character(),
  PhenotypeIDS = col_character(),
  PhenotypeList = col_character(),
  Origin = col_character(),
  OriginSimple = col_character(),
  Assembly = col_character(),
  ChromosomeAccession = col_character(),
  chr.clinvar = col_character(),
  start.clinvar = col_character(),
  stop.clinvar = col_character(),
  ReferenceAllele = col_character(),
  AlternateAllele = col_character(),
  Cytogenetic = col_character(),
  ReviewStatus = col_character(),
  NumberSubmitters = col_character(),
  Guidelines = col_character(),
  TestedInGTR = col_character(),
  OtherIDs = col_character(),
  SubmitterCategories = col_character(),
  VariationID = col_character()
)
columns_to_pull <-  "#CHROM,FROM,TO,Annotation_Impact,Gene_ID,Feature_Type,Feature_ID,Transcript_BioType,Rank,HGVS.c,HGVS.p,LOF,NMD,MA,DBSNP,X1000Gen,Cosmic,genecode_gene_name,genecode_strand,Type,Name,GeneID,GeneSymbol,HGNC_ID,ClinicalSignificance,ClinSigSimple,LastEvaluated,RCVaccession,PhenotypeIDS,PhenotypeList,Origin,OriginSimple,Assembly,ChromosomeAccession,Cytogenetic,ReviewStatus,NumberSubmitters,Guidelines,TestedInGTR,OtherIDs,SubmitterCategories,VariationID"

columns_to_pull <-  "#CHROM,FROM,TO,LOF"

#columns_to_pull <-  "#CHROM,FROM,TO,Annotation_Impact"

columns_to_pull <- columns_to_pull %>% str_split(",") %>% unlist 



if(length(args) == 0){
	args <- "Rscript projects/molonc/aparicio_lab/sbeatty/yvr_pipelines/vcf_manipulation/VariantsToTable_wrapper.R --target_file=filter_low_qual/SA1145N.SING_MUTA.fixed_header_gatk.fixed_header_twice.remove.low.qual.vcf --output_file=test.tsv --annotation_file=annotation_file_fixed_columns/SA1145N.SING_MUTA.annotation_file.fixed_columns.noheader.tsv.gz"
	args <- str_split(args," ") %>% unlist
}

target_file_path <- args[str_detect(args,"--target_file")] %>% str_replace("--target_file=","")
output_file_path <- args[str_detect(args,"--output_file")] %>% str_replace("--output_file=","")
annotation_file <- args[str_detect(args,"--annotation_file")] %>% str_replace("--annotation_file=","")


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


annotation_file_df <- read_tsv(annotation_file, col_types = annotation_col_specs )

names(annotation_file_df)

annotation_column_vector <- names(annotation_file_df)

annotation_column_vector[!c(annotation_column_vector %in% columns_to_pull)] <- "-"

paste("bcftools annotate -a",annotation_file,"-h", "test_header_lines.hdr", "-c",paste(str_replace(annotation_column_vector,"#CHROM","CHROM"), collapse=","), 
target_file_path,">",output_file_path)

paste("bcftools annotate -a",annotation_file,"-h", "test_header_lines.hdr", "-c",paste(str_replace(names(annotation_file_df),"#CHROM","CHROM"), collapse=","), 
target_file_path,">",output_file_path)

info_tags <- mapply(info_lines, USE.NAMES=FALSE, FUN=parse_info_line)
standard_tags <- header %>% last %>% str_split("\t") %>% unlist %>% str_replace("#","")
GATK_tags <- c("EVENTLENGTH", "TRANSITION", "HET", "HOM-REF", "HOM-VAR", "NO-CALL", "TYPE", "VAR", "NSAMPLES", "NCALLED", "MULTI-ALLELIC")
tags <- c(standard_tags, info_tags, GATK_tags)
tags <- tags  %>% paste("-F",.) %>% paste(collapse=" ")

GATK_command_call <- paste("/projects/molonc/aparicio_lab/sbeatty/software/gatk-4.1.8.1/gatk VariantsToTable -V", target_file_path,"--show-filtered", tags,"-O", output_file_path)
print(GATK_command_call)
system(paste(GATK_command_call))


 shell("bcftools annotate -a {input.annotation_in} -h {input.new_header_lines} -c CHROM,FROM,TO,TB_ID,-,gencode_gene_name,-,-,-,germline_or_somatic,task,transcript_status,-,-,-,-,-,-,-,-,-,gencode_feature_type,Clinvar_AlleleID,Type,Name,GeneID,HGNC_ID,ClinicalSignificance,-,Assembly,Cytogenetic,NumberSubmitters,Clinvar_VariationID,cosmic_mutation_strand,COSMIC_GENOMIC_MUTATION_ID,cosmic_patients_total,cosmic_patients_gwas,snpeff_impact,stop_codon_gain_potential,stop_codon_loss_potential,CCDS,FATHMM.prediction,FATHMM.score,SIFT_effect_prediction,SIFT_substitution_type {input.vcf_in} > {output}")

bcftools annotate -a annotation_file_fixed_columns/SA1145N.SING_MUTA.annotation_file.fixed_columns.noheader.tsv.gz -c CHROM,FROM,TO,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,Gene_ID,Feature_Type,Feature_ID,Transcript_BioType,Rank,HGVS.c,HGVS.p,-,-,-,-,-,LOF,NMD,MA,DBSNP,X1000Gen,Cosmic,-,-,-,-,-,-,-,-,-,-,-,-,genecode_gene_name,genecode_strand,-,-,Type,Name,GeneID,GeneSymbol,HGNC_ID,ClinicalSignificance,ClinSigSimple,LastEvaluated,-,-,RCVaccession,PhenotypeIDS,PhenotypeList,Origin,OriginSimple,Assembly,ChromosomeAccession,-,-,-,-,-,Cytogenetic,ReviewStatus,NumberSubmitters,Guidelines,TestedInGTR,OtherIDs,SubmitterCategories,VariationID filter_low_qual/SA1145N.SING_MUTA.fixed_header_gatk.fixed_header_twice.remove.low.qual.vcf > test.tsv
#system(paste("/projects/molonc/aparicio_lab/sbeatty/software/gatk-4.1.8.1/gatk VariantsToTable -V", target_file_path,"--show-filtered", tags,"-O test_vanilla.tsv"))

#system("/projects/molonc/aparicio_lab/sbeatty/software/gatk-4.1.8.1/gatk VariantsToTable -V /projects/molonc/aparicio_lab/sbeatty/IND/IND-31/annotate_task6/filter_low_qual/SA1284T.SING_STRE.fixed_header.remove.low.qual.vcf --show-filtered -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F INFO -F FORMAT -F SA1284T -F END -F BLOCKAVG_min30p3a -F SNVHPOL -F CIGAR -F RU -F REFREP -F IDREP -F MQ -F EVENTLENGTH -F TRANSITION -F HET -F HOM-REF -F HOM-VAR -F NO-CALL -F TYPE -F VAR -F NSAMPLES -F NCALLED -F MULTI-ALLELIC -O vanilla_test.tsv")

#system("/projects/molonc/aparicio_lab/sbeatty/software/gatk-4.1.8.1/gatk VariantsToTable -V /projects/molonc/aparicio_lab/sbeatty/IND/IND-31/annotate_task6/filter_low_qual/SA1284T.SING_STRE.fixed_header.remove.low.qual.vcf -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F INFO -F FORMAT -F SA1284T -F END -F BLOCKAVG_min30p3a -F SNVHPOL -F CIGAR -F RU -F REFREP -F IDREP -GF MQ -F EVENTLENGTH -F TRANSITION -F HET -F HOM-REF -F HOM-VAR -F NO-CALL -F TYPE -F VAR -F NSAMPLES -F NCALLED -F MULTI-ALLELIC --show-filtered -O test_move_show_filetered.tsv")

#/projects/molonc/aparicio_lab/sbeatty/software/gatk-4.1.8.1/gatk VariantsToTable -V /projects/molonc/aparicio_lab/sbeatty/IND/IND-31/annotate_task6/filter_low_qual/SA1284T.SING_STRE.fixed_header.remove.low.qual.vcf -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F INFO -F FORMAT -F SA1284T -F END -F BLOCKAVG_min30p3a -F SNVHPOL -F CIGAR -F RU -F REFREP -F IDREP -GF MQ -F EVENTLENGTH -F TRANSITION -F HET -F HOM-REF -F HOM-VAR -F NO-CALL -F TYPE -F VAR -F NSAMPLES -F NCALLED -F MULTI-ALLELIC --show-filtered -O test_csv.csv