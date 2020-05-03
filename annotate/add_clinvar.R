# script is setup for GRCh37
# caveats will match clinvar indels that are the same length as the sample indels. SNV matches are exact. 
# in cases where clinvar has a different reference variant than the sample 
# reference, and the variants are 1 bp, the clinvar match will be ignored
.libPaths("/ssd/sda1/sbeatty/software/miniconda3/lib/R/library")
options(echo=TRUE)
options(verbose=FALSE)
require("stringr", quietly=TRUE)
require("readr", quietly=TRUE)
require("dplyr", quietly=TRUE)
require("parallel", quietly=TRUE)
require("VariantAnnotation", quietly=TRUE)
require("data.table", quietly=TRUE)

args = commandArgs(trailingOnly=TRUE)

# supply column specifications to readr greatly speeds up file reading
tsv_column_spec <- 
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

clinvar_column_specs <- cols(
  `#AlleleID` = col_double(),
  Type = col_character(),
  Name = col_character(),
  GeneID = col_character(),
  GeneSymbol = col_character(),
  HGNC_ID = col_character(),
  ClinicalSignificance = col_character(),
  ClinSigSimple = col_double(),
  LastEvaluated = col_character(),
  `RS# (dbSNP)` = col_double(),
  `nsv/esv (dbVar)` = col_character(),
  RCVaccession = col_character(),
  PhenotypeIDS = col_character(),
  PhenotypeList = col_character(),
  Origin = col_character(),
  OriginSimple = col_character(),
  Assembly = col_character(),
  ChromosomeAccession = col_character(),
  Chromosome = col_character(),
  Start = col_double(),
  Stop = col_double(),
  ReferenceAllele = col_character(),
  AlternateAllele = col_character(),
  Cytogenetic = col_character(),
  ReviewStatus = col_character(),
  NumberSubmitters = col_double(),
  Guidelines = col_character(),
  TestedInGTR = col_character(),
  OtherIDs = col_character(),
  SubmitterCategories = col_double(),
  VariationID = col_double()
)

if(length(args) == 0){
args <- "Rscript //scratch/shahlab_tmp/sbeatty/yvr_pipelines/annotate/add_clinvar.R --targetfile=/scratch/shahlab_tmp/sbeatty/ind231/clinvar_in/GERM_STRE_INDEL_SA1228T.vcf.gz --outputfile=clinvar_added/GERM_STRE_INDEL_SA1228T.clinvar.csv --loadreference=TRUE"
args <- str_split(args," ") %>% unlist
}
cores <- args[str_detect(args,"--cores")] %>% str_replace("--cores=","")
if(length(cores) == 0){
  cores <- 1
} else {
  cores <- as.numeric(cores)
}

rangeordf <- args[str_detect(args,"--rangeordf")] %>% str_replace("--rangeordf=","")
if(length(rangeordf) == 0){
  rangeordf <- "df"
}

target_file <- args[str_detect(args,"--targetfile")] %>% str_replace("--targetfile=","")
output_file <- args[str_detect(args,"--outputfile")] %>% str_replace("--outputfile=","")
load_parsed_reference_data <- args[str_detect(args,"--loadreference")] %>% str_replace("--loadreference=","")



target_file_type <-gsub( "vcf.gz","vcf",basename(target_file)) %>% str_split("\\.") %>% unlist() %>% last()

#load_parsed_reference_data <- TRUE
if(length(load_parsed_reference_data) == 0){
  load_parsed_reference_data <- FALSE
}

if(load_parsed_reference_data == TRUE & c( !file.exists("/scratch/shahlab_tmp/sbeatty/reference/clinvar_reference_dat_ranges.Rdata") |
  !file.exists("/scratch/shahlab_tmp/sbeatty/reference/clinvar_reference_dat.Rdata") |
  !file.exists("/scratch/shahlab_tmp/sbeatty/reference/genecode_reference.Rdata") |
  !file.exists("/scratch/shahlab_tmp/sbeatty/reference/gencode_reference_dat_ranges.Rdata") )){
  load_parsed_reference_data <- FALSE
}


if(load_parsed_reference_data == FALSE){
  clinvar_reference_dat <- read_tsv("/shahlab/archive/misc/sbeatty/reference/variant_summary.txt", col_types=clinvar_column_specs)
  clinvar_reference_dat <- clinvar_reference_dat[c(clinvar_reference_dat[,"Assembly"] == "GRCh37"),] %>% data.frame
 
# elminate records with a start position less than zero. As of April 21
# 2019 this applies to only 2 results

  clinvar_reference_dat <- clinvar_reference_dat[!c(clinvar_reference_dat$Stop < 0),]
  clinvar_reference_dat <- clinvar_reference_dat[!is.na(clinvar_reference_dat$Chromosome),]
  #clinvar_reference_dat[,"inferred_ref"] <- NA
  #clinvar_reference_dat[,"inferred_alt"] <- NA
  clinvar_reference_dat_ranges <-  makeGRangesFromDataFrame(clinvar_reference_dat, seqnames.field ="Chromosome" , start.field="Start", end.field="Stop", keep.extra.columns=TRUE)

    .gene_from_add_info <- function(x){
    output <- NA
    add_info_type <- unlist(strsplit(x, ";")) %>% trimws()
    matches <- str_detect(add_info_type, "gene_name")
    if(length(matches == 1)){
      output <- add_info_type[str_detect(add_info_type, "gene_name")] %>% str_replace("gene_name ","")
    }
    output 
  }

  # five rows are skipped in order to skip the gtf header
  gencode_reference_dat <- read_tsv("/shahlab/archive/misc/sbeatty/reference/gencode.v19.annotation.gtf_withproteinids", skip=5, col_names=FALSE,quote="XXX", col_types=gencode_reference_col_types)
  names(gencode_reference_dat) <- c("chr", "source", "feature_type","start","stop", "score", "strand", "phase", "add_info")
  gencode_reference_dat$chr <- gsub("chr","", gencode_reference_dat$chr)
  gencode_reference_dat$add_info <- gsub(regex('\\"'),"",gencode_reference_dat$add_info)
  gencode_reference_dat[,"gencode_strand"] <- gencode_reference_dat$strand
  gencode_reference_dat[,"genecode_gene_name"] <-mcmapply(gencode_reference_dat$add_info, FUN=.gene_from_add_info, USE.NAMES = FALSE, mc.cores=cores)
  
  gencode_reference_dat_ranges <-    makeGRangesFromDataFrame(data.frame(gencode_reference_dat),seqnames.field ="chr" , start.field="start", end.field="stop", keep.extra.columns=TRUE)

  save(clinvar_reference_dat_ranges,file="/scratch/shahlab_tmp/sbeatty/reference/clinvar_reference_dat_ranges.Rdata")
  save(clinvar_reference_dat, file="/scratch/shahlab_tmp/sbeatty/reference/clinvar_reference_dat.Rdata")
  save(gencode_reference_dat, file="/scratch/shahlab_tmp/sbeatty/reference/genecode_reference.Rdata")
  save(gencode_reference_dat_ranges, file="/scratch/shahlab_tmp/sbeatty/reference/gencode_reference_dat_ranges.Rdata")

} else {
  load("/scratch/shahlab_tmp/sbeatty/reference/clinvar_reference_dat_ranges.Rdata")
  load("/scratch/shahlab_tmp/sbeatty/reference/clinvar_reference_dat.Rdata")
  load("/scratch/shahlab_tmp/sbeatty/reference/genecode_reference.Rdata")
  load("/scratch/shahlab_tmp/sbeatty/reference/gencode_reference_dat_ranges.Rdata")
}

# conversion from vcf object to a data frame even when only 
# converting a single tab (ex INFO) is very slow
# function below is > 1000x fold faster

 source("/scratch/shahlab_tmp/sbeatty/yvr_pipelines/annotate/core_functions.R")

if(target_file_type == "vcf"){
  target_dat <- readVcf(paste(target_file), genome="hg19")
  caller <- ".variant_caller"
  target_dat_ranges <- target_dat
  target_dat_df <- genome_range_to_dataframe(target_dat, extract_reads=TRUE)
  names(target_dat_df) <- gsub("stop", "end", names(target_dat_df))
  target_dat_df[,"genecode_gene_name"] <- NA
  target_dat_df[,"genecode_strand"] <- NA
} else if (target_file_type == "tsv"){
  caller <- ".variant_caller"
  target_dat_df <- read_tsv(paste(target_file), col_types=tsv_column_spec)
  target_dat_df <- target_dat_df %>% data.frame
    names(target_dat_df) <- gsub("chromosome", "chr", names(target_dat_df))
    names(target_dat_df) <- gsub("stop", "end", names(target_dat_df))
    target_dat_df[,"genecode_gene_name"] <- NA
  target_dat_df[,"genecode_strand"] <- NA
  target_dat_df <- target_dat_df[!is.na(target_dat_df$chr),]
  target_dat_ranges <- makeGRangesFromDataFrame(target_dat_df, 
       seqnames.field ="chr" , start.field="start", end.field="end", keep.extra.columns=TRUE)   
} else if (target_file_type == "csv"){
  caller <- ".variant_caller"
  target_dat_df <- read_csv(paste(target_file), col_types=tsv_column_spec)
  target_dat_df <- target_dat_df %>% data.frame
    names(target_dat_df) <- gsub("chromosome", "chr", names(target_dat_df))
    names(target_dat_df) <- gsub("stop", "end", names(target_dat_df))
    target_dat_df[,"genecode_gene_name"] <- NA
  target_dat_df[,"genecode_strand"] <- NA
  target_dat_df <- target_dat_df[!is.na(target_dat_df$chr),]
  target_dat_ranges <- makeGRangesFromDataFrame(target_dat_df, 
       seqnames.field ="chr" , start.field="start", end.field="end", keep.extra.columns=TRUE)   
} 
gencode_matches <- findOverlaps(target_dat_ranges, gencode_reference_dat_ranges, type="any")
target_dat_df[,"genecode_gene_name"][gencode_matches@from] <- gencode_reference_dat_ranges$genecode_gene_name[gencode_matches@to]
target_dat_df[,"genecode_strand"][gencode_matches@from] <- gencode_reference_dat_ranges$gencode_strand[gencode_matches@to]

clin_var_matches <- findOverlaps(target_dat_ranges, clinvar_reference_dat_ranges, type="equal") %>% data.frame

target_ref_column <- names(target_dat_df)[str_detect(names(target_dat_df) %>% toupper,"REF")][1]
target_alt_column <- names(target_dat_df)[str_detect(names(target_dat_df) %>% toupper,"ALT")][1]

clin_var_matches <- extract_exact_matches(matches_df =clin_var_matches, query_df=target_dat_df, query_ref_allele_column=target_ref_column, query_alt_allele_column=target_alt_column, subject_df=clinvar_reference_dat,subject_ref_allele_column="ReferenceAllele",  subject_alt_allele_column="AlternateAllele")


if(nrow(clin_var_matches) > 0){

  target_dat_df[,"unique_id"] <- paste0(target_dat_df$chr,"@",target_dat_df$end)
  clinvar_reference_dat$unique_id <- NA
  clinvar_reference_dat$unique_id[clin_var_matches$subjectHits] <- target_dat_df$unique_id[clin_var_matches$queryHits]
  names(clinvar_reference_dat) <- substitute_column_name("Chromosome", "chr.clinvar", clinvar_reference_dat)
  names(clinvar_reference_dat) <- substitute_column_name("Start", "start.clinvar", clinvar_reference_dat)
  names(clinvar_reference_dat) <- substitute_column_name("Stop", "stop.clinvar", clinvar_reference_dat)


  names(target_dat_df) <- substitute_column_name("start", "start.variant_caller", target_dat_df)
  names(target_dat_df) <- substitute_column_name("end", "end.variant_caller", target_dat_df)
  names(target_dat_df) <- substitute_column_name("chr", "chr.variant_caller", target_dat_df)

  target_dat_clinvar_df <- left_join(data.frame(target_dat_df), data.frame(clinvar_reference_dat), by="unique_id", suffix=c(paste(caller),"clinvar"), KEEP=TRUE)
  } else {
    target_dat_clinvar_df <- target_dat_df
  }


#target_dat_clinvar_df <- target_dat_clinvar_df[,!c(names(target_dat_clinvar_df) %in% c("Chromosome", "Start", "Stop"))] %>% head
names(target_dat_clinvar_df) <- substitute_column_name("ALT","alt", target_dat_clinvar_df)
names(target_dat_clinvar_df) <- substitute_column_name("REF","ref", target_dat_clinvar_df)
names(target_dat_clinvar_df) <- substitute_column_name("start.variant_caller","start", target_dat_clinvar_df)
names(target_dat_clinvar_df) <- substitute_column_name("end.variant_caller","end", target_dat_clinvar_df)
names(target_dat_clinvar_df) <- substitute_column_name("chr.variant_caller","chr", target_dat_clinvar_df)
reordered_columns  <- c("chr", "start", "end",names(target_dat_clinvar_df)[!c(names(target_dat_clinvar_df) %in% c("chr", "start", "end"))])

target_dat_clinvar_df <- target_dat_clinvar_df[,reordered_columns]

if(rangeordf == "range"){
  joined_ranges <- makeGRangesFromDataFrame(target_dat_clinvar_df, seqnames.field ="chr" , start.field='start', end.field='end', keep.extra.columns=TRUE)
  save(joined_ranges,file=output_file)  
} else {
  fwrite(target_dat_clinvar_df, file=output_file)
}
