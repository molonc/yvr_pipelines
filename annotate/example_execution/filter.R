#!/usr/bin/env Rscript
.libPaths("/ssd/sda1/sbeatty/software/miniconda3/lib/R/library")


##filtering for task 3 

options(echo=TRUE)
options(verbose=FALSE)
args = commandArgs(trailingOnly=TRUE)

library("stringr")
library("dplyr")
library("readr")
library("VariantAnnotation")
library("data.table")



cols_tsv <- 
cols(
  .default = col_character(),
  case_id = col_character(),
  normal_id = col_character(),
  tumour_id = col_character(),
  chromosome = col_character(),
  start = col_double(),
  stop = col_double(),
  gene = col_character(),
  gene_id = col_character(),
  type = col_character(),
  filter = col_character(),
  ref = col_character(),
  alt = col_character(),
  gt = col_character(),
  pl = col_character(),
  mut_pr = col_double(),
  tr = col_double(),
  ta = col_double(),
  nr = col_double(),
  na = col_double(),
  dbsnp = col_character(),
  thousand_genomes = col_character(),
  cosmic = col_character(),
  caller = col_character(),
  trinucleotide_ref = col_character(),
  trinucleotide_alt = col_character(),
  amino_acid_change = col_character(),
  functional_class = col_character(),
  gene_coding = col_character(),
  project = col_character(),
  substitution = col_character(),
  impact = col_character(),
  dna_change = col_character(),
  biotype = col_character()
)

get_anno_columns <- function(header){
  for(i in 1:nrow(header)){
    if(i == 1){
      anno_columns <- NA
    }
    if(str_detect(header[i,], "##INFO=<ID=ANN") == TRUE){
      anno_columns <- header[i,]
    }
  }
  anno_columns <- gsub('##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: \'',"", anno_columns)
  anno_columns <- anno_columns %>% str_split('\\|') %>% unlist %>% trimws
  anno_columns
}

extract_anno_results <- function(info_result){
  anno_info_i_split <- info_result %>% str_split(";") ; info_result <- info_result[[1]]
  ANN_match <- mapply(anno_info_i_split, USE.NAMES=FALSE, FUN=function(x){str_detect(x,"ANN=")})
  ANN_match <- mapply(anno_info_i_split, FUN=function(x){str_detect(x,"ANN=")})[,1]
  anno_result <- unlist(anno_info_i_split)[ANN_match]
  anno_result
}


impact_from_anno_table <- function(anno_result_compact){
  impact <- NA
  if(!is.na(anno_result_compact)){
  info_columns <- c('Allele', 'Annotation', 'Annotation_Impact', 'Gene_Name', 'Gene_ID', 'Feature_Type', 'Feature_ID', 'Transcript_BioType', 'Rank', 'HGVS.c', 'HGVS.p', 'cDNA.pos / cDNA.length', 'CDS.pos / CDS.length', 'AA.pos / AA.length', 'Distance', 'ERRORS / WARNINGS / INFO')
  anno_result_compact <- extract_anno_results(anno_result_compact)
  anno_result_split <- anno_result_compact %>% str_split(",") %>% unlist()
  for(i in 1:length(anno_result_split)){
    if(i == 1){
      out <- str_split(anno_result_split[i],'\\|')[[1]]
    } else if (i > 1){
      out_i <- str_split(anno_result_split[i],'\\|')[[1]]
      out <- rbind(out, out_i)
    }
  }

  anno_result_split_df <- data.frame(out)
  if(min(dim(anno_result_split_df)) > 1){
    names(anno_result_split_df) <- info_columns
    impact <- unique(anno_result_split_df$Annotation_Impact)
    impact <- paste(impact, collapse = ",")
  } else {
    anno_result_split_df <- anno_result_split_df %>% t %>% data.frame
    names(anno_result_split_df) <- info_columns
    impact <- unique(anno_result_split_df$Annotation_Impact)
    impact <- paste(impact, collapse = ",")
  }
}
  impact
}



calculate_vaf <- function(sample_df){
  tumour_id <- sample_df$tumour_id %>% unique
  normal_id <- sample_df$normal_id %>% unique
  if(normal_id[1] == "None" | tumour_id[1] == "None"){
    if(normal_id[1] == "None"){
      vaf <- sample_df$ta/c(sample_df$ta + sample_df$tr)
      } else if (tumour_id[1] == "None"){
        vaf <- sample_df$na/c(sample_df$na + sample_df$nr)
      } 
      } else {
        vaf <- rep(NA,nrow(sample_df))
        }
        sample_df[,"vaf"] <- vaf
        sample_df
      }

convert_snpeff_impact <- function(impact_i){
  if(!is.na(impact_i)){
    if(str_detect(impact_i,"HIGH")){
    impact_output <- "HIGH"
  } else if (str_detect(impact_i,"MODERATE")){
    impact_output <- "MODERATE"
  } else if (str_detect(impact_i,"LOW")){
    impact_output <- "LOW"
  } else if (str_detect(impact_i,"MODIFIER")){
    impact_output <- "MODIFIER"
  } else {
    impact_output <- NA
  }
  } else {
    impact_output <- NA
  }
  
  impact_output
}


source("/scratch/shahlab_tmp/sbeatty/yvr_pipelines/annotate/core_functions.R")
target_file <- args[str_detect(args,"--targetfile")] %>% str_replace("--targetfile=","")

output_file <- args[str_detect(args,"--outputfile")] %>% str_replace("--outputfile=","")
ncpus <- args[str_detect(args,"--cpu_count")] %>% str_replace("--cpu_count=","")

if(length(args) == 0){
  target_file <- "exons_flagged/SA609.clinvar.cosmic.exons.csv"
  output_file <- "test.csv"
  ncpus <- 20
}
target_file_type <-gsub( "vcf.gz","vcf",basename(target_file)) %>% str_split("\\.") %>% unlist() %>% last()
target_genes <- c("TP53", "PTEN", "MYC", "BRCA1", "BRCA2","PIK3CA", "PTEN")

#import referenc sequences
gencode_reference_col_types = cols(
  X1 = col_character(),
  X2 = col_character(),
  X3 = col_character(),
  X4 = col_double(),
  X5 = col_double(),
  X6 = col_character(),
  X7 = col_character(),
  X8 = col_character(),
  X9 = col_character()
)


Pathogenic_classes <- c("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic","risk factor","Pathogenic/Likely pathogenic, drug response"," Pathogenic, risk factor","Likely pathogenic, drug response","Pathogenic, drug response","Likely pathogenic, risk factor","Pathogenic/Likely pathogenic, other","Pathogenic, Affects","Pathogenic, association","Pathogenic, other","Likely pathogenic, other, risk factor","Pathogenic/Likely pathogenic, risk factor","Likely pathogenic, other","Likely pathogenic, association")


df <- read_csv(target_file, col_types=cols_tsv)  %>% data.frame
df[,"target_gene_flag"] <- FALSE
df[c(df$genecode_gene_name %in% target_genes),"target_gene_flag"] <- TRUE
names(df) <- substitute_column_name("seqnames", "chr", df)
names(df) <- gsub(".variantcaller","",names(df))

if(length(df) > 1){
matches <- which(c(c(!is.na(df$ClinicalSignificance) & c(df$ClinicalSignificance %in% Pathogenic_classes)) | c(df$GENOMIC_MUTATION_ID != "")) == TRUE)
df <- df[matches, ]

if(length(matches) == 0){
   #print("error 1")
   fwrite(data.frame("no valid results"), file=output_file)
   quit(save="no", status=0)
} 
fwrite(data.frame(df), file=output_file)
} else {
  #print("error 2")
  fwrite(data.frame("no valid results"), file=output_file)
  quit(save="no", status=0)
}




