#!/usr/bin/env Rscript
options(verbose=FALSE)
args = commandArgs(trailingOnly=TRUE)
require("stringr", quietly=TRUE)
require("readr", quietly=TRUE)
require("dplyr", quietly=TRUE)


column_specs <- cols(
  chr = col_character(),
  start = col_double(),
  end = col_double(),
  ID = col_character(),
  ref = col_character(),
  alt = col_character(),
  QUAL = col_double(),
  FILTER = col_character(),
  INFO = col_character(),
  PR = col_double(),
  TC = col_character(),
  TR = col_double(),
  TA = col_double(),
  NR = col_double(),
  V14 = col_double(),
  ND = col_double(),
  NI = col_double(),
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
  Distance = col_double(),
  `ERRORS / WARNINGS / INFO` = col_character(),
  LOF = col_character(),
  NMD = col_character(),
  MA = col_character(),
  DBSNP = col_character(),
  X1000Gen = col_character(),
  Cosmic = col_character(),
  EVENTLENGTH = col_double(),
  TRANSITION = col_double(),
  HET = col_double(),
  HOM.REF = col_double(),
  HOM.VAR = col_double(),
  NO.CALL = col_double(),
  TYPE = col_character(),
  VAR = col_double(),
  NSAMPLES = col_double(),
  NCALLED = col_double(),
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
  VariationID = col_character(),
  SIFT_effect_prediction= col_character(),
  SIFT_substitution_type=  col_character()
)

if(length(args) == 0){
	args <- "Rscript /projects/molonc/aparicio_lab/sbeatty/yvr_pipelines/vcf_manipulation/sort_annotation.R --annotation_file_in=//projects/molonc/aparicio_lab/sbeatty/BXE/BXE-264/annotation_task6/simplified_snpeff/SA1019X1_1L.SING_MUTA.clinvar.PR85.simple_snpeff.csv --annotation_file_out=test.tsv"
	args <- str_split(args," ") %>% unlist
}

annotation_file_in <- args[str_detect(args,"--annotation_file_in")] %>% str_replace("--annotation_file_in=","")
annotation_file_out <- args[str_detect(args,"--annotation_file_out")] %>% str_replace("--annotation_file_out=","")

annotation_file_ext <- annotation_file_in %>% str_split("/") %>% unlist %>%  last %>%  str_split("\\.") %>% unlist %>% last

#annotation_file[,"chr_modified"] <- annotation_file$chr %>% str_replace_all("Y","24") %>% str_replace_all("X","23")


if(annotation_file_ext == "tsv")
	{
		annotation_file <- read_tsv(annotation_file_in, col_types=column_specs)
    annotation_file[,"chr_modified"] <- annotation_file$chr %>% str_replace_all("Y","24") %>% str_replace_all("X","23")
    annotation_file <- annotation_file %>% arrange(.,chr_modified, start, end) %>% dplyr::select(-chr_modified)
	} else if(annotation_file_ext == "csv"){
		annotation_file <- read_csv(annotation_file_in, col_types=column_specs)
    annotation_file[,"chr_modified"] <- annotation_file$chr %>% str_replace_all("Y","24") %>% str_replace_all("X","23")
    annotation_file <- annotation_file %>% arrange(.,chr_modified, start, end) %>% dplyr::select(-chr_modified)
	} else {
		print("error reading annotation file extension")
		q()
	}




for(col_i in names(annotation_file)){
  annotation_file[,paste(col_i)][[1]] <- annotation_file[,paste(col_i)][[1]] %>% str_replace_all(" ","_")
}

annotation_file <- annotation_file %>% rename(CHROM=chr, FROM=start, TO=end)

annotation_file <- annotation_file %>% relocate("TO") %>% relocate("FROM") %>% relocate("CHROM")

names(annotation_file)[1] <- "#CHROM"
write_csv(annotation_file, annotation_file_out, col_names=TRUE)