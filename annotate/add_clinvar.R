# script is setup for GRCh37
# caveats will match clinvar indels that are the same length as the sample indels. SNV matches are exact. 
# in cases where clinvar has a different reference variant than the sample 
# reference, and the variants are 1 bp, the clinvar match will be ignored
.libPaths("/ssd/sda1/sbeatty/software/miniconda3/lib/R/library")

require("stringr")
require("readr")
require("dplyr")
require("parallel")
require("VariantAnnotation")
options(echo=TRUE)
options(verbose=TRUE)
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

# Example command line inputs
#args <- c("--loadreference=FALSE", "--targetfile=clinvar_in/GERM_STRE_INDEL_SA1267NC.vcf.gz", "--outputfile=testfile.csv", "--cores=30", "--rangeordf=range")

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
target_file_type <-gsub( "vcf.gz","vcf",basename(target_file)) %>% str_split("\\.") %>% unlist() %>% last()
load_parsed_reference_data <- args[str_detect(args,"--loadreference")] %>% str_replace("--loadreference=","")

load_parsed_reference_data <- FALSE
if(load_parsed_reference_data == FALSE){



  clinvar_reference_dat <- read_tsv("/shahlab/archive/misc/sbeatty/reference/variant_summary.txt", col_types=clinvar_column_specs)
  clinvar_reference_dat <- clinvar_reference_dat[c(clinvar_reference_dat[,"Assembly"] == "GRCh37"),] %>% data.frame
  clinvar_reference_dat[,"inferred_strand"] <- NA

# elminate records with a start position less than zero. As of April 21
# 2019 this applies to only 2 results

  clinvar_reference_dat <- clinvar_reference_dat[!c(clinvar_reference_dat$Stop < 0),]
  clinvar_reference_dat <- clinvar_reference_dat[!is.na(clinvar_reference_dat$Chromosome),]
  clinvar_reference_dat[,"inferred_ref"] <- NA
  clinvar_reference_dat[,"inferred_alt"] <- NA
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
  save(gencode_reference_dat_ranges, file="gencode_reference_dat_ranges.Rdata")

} else {
  load("/scratch/shahlab_tmp/sbeatty/reference/clinvar_reference_dat_ranges.Rdata")
  load("/scratch/shahlab_tmp/sbeatty/reference/clinvar_reference_dat.Rdata")
  load("/scratch/shahlab_tmp/sbeatty/reference/genecode_reference.Rdata")
  load("/scratch/shahlab_tmp/sbeatty/reference/gencode_reference_dat_ranges.Rdata")
}

# conversion from vcf object to a data frame even when only 
# converting a single tab (ex INFO) is very slow
# function below is > 1000x fold faster

genome_range_to_dataframe <- function(genome_range){
  df <- genome_range
  df_rowRanges <- rowRanges(df)
  range_start <- ranges(df)@start
  range_end <- ranges(df)@start
  range_width <- ranges(df)@width
  df_seqnames <-  as.character(seqnames(df))
  df_elementMetadata <- elementMetadata(rowRanges(df))
  strand_info <- as.character(df_rowRanges@strand)
  df_paramRangeID <- as.character(df_elementMetadata[,"paramRangeID"])
  df_REF <- as.character(df_elementMetadata[,"REF"])
  df_ALT <- rownames(df_elementMetadata) %>% str_sub(start=-1) %>% as.character
  df_QUAL <- as.character(df_elementMetadata[,"QUAL"])
  df_FILTER <- as.character(df_elementMetadata[,"FILTER"])
  df_out <- data.frame(chr=df_seqnames, pos=range_start, start = range_start,
    end=range_end, width=range_width, strand=strand_info,
    paramRangeID=df_paramRangeID, REF=df_REF, ALT=df_ALT, QUAL=df_QUAL, Filter=df_FILTER)
  df_out
}


if(target_file_type == "vcf"){
  target_dat <- readVcf(paste(target_file), genome="hg19")
  caller <- ".variant_caller"
  target_dat_ranges <- target_dat
  target_dat_df <- genome_range_to_dataframe(target_dat)
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
    
} 

gencode_matches <- findOverlaps(target_dat_ranges, gencode_reference_dat_ranges, type="any")
target_dat_df[,"genecode_gene_name"][gencode_matches@from] <- gencode_reference_dat_ranges$genecode_gene_name[gencode_matches@to]
target_dat_df[,"genecode_strand"][gencode_matches@from] <- gencode_reference_dat_ranges$gencode_strand[gencode_matches@to]

clin_var_matches <- findOverlaps(target_dat_ranges, clinvar_reference_dat_ranges, type="equal") %>% data.frame

clin_var_matches[,"gencode_strand"] <- target_dat_df[clin_var_matches[,"queryHits"],"genecode_strand"] 
clin_var_matches[,"query_ref"] <- target_dat_df[clin_var_matches[,"queryHits"],"ref"] 
clin_var_matches[,"query_alt"] <- target_dat_df[clin_var_matches[,"queryHits"],"alt"] 

clin_var_matches[,"subject_ref"] <- clinvar_reference_dat[clin_var_matches[,"subjectHits"],"ReferenceAllele"] 

clin_var_matches[,"subject_alt"] <- clinvar_reference_dat[clin_var_matches[,"subjectHits"],"AlternateAllele"] 


# eliminate results where the snv reference for clinvar and museq vary
clin_var_matches <- clin_var_matches[!c(clin_var_matches$query_ref != clin_var_matches$subject_ref & nchar(clin_var_matches$query_ref) == nchar(clin_var_matches$subject_ref) & nchar(clin_var_matches$query_ref) ==1 ),]




if(nrow(clin_var_matches) > 0){

# include exact matches
  exact_matches <- c(clin_var_matches[,"query_alt"] == clin_var_matches[,"subject_alt"])

# also include multi basepair variants of the same length
  indel_length_matches <- c(nchar(clin_var_matches[,"query_alt"]) == nchar(clin_var_matches[,"subject_alt"]) & nchar(clin_var_matches[,"query_alt"]) > 1)

  clin_var_matches <-  clin_var_matches[exact_matches | indel_length_matches,]

  target_dat_df[,"unique_id"] <- paste0(target_dat_df$chr,"@",target_dat_df$end)
  clinvar_reference_dat$unique_id <- NA
  clinvar_reference_dat$unique_id[clin_var_matches$subjectHits] <- target_dat_df$unique_id[clin_var_matches$queryHits]
  names(target_dat_df) <- gsub(tolower("start"), "start.variant_caller", names(target_dat_df))
  names(target_dat_df) <- gsub(tolower("end"), "end.variant_caller", names(target_dat_df))
  names(target_dat_df) <- gsub(tolower("chr"), "chr.variant_caller", names(target_dat_df))
  target_dat_clinvar_df <- left_join(data.frame(target_dat_df), data.frame(clinvar_reference_dat), by="unique_id", suffix=c(paste(caller),"clinvar"), KEEP=TRUE)
  } else {
    target_dat_clinvar_df <- target_dat_df
  }

names(target_dat_clinvar_df) <- gsub("ALT","alt", names(target_dat_clinvar_df))

target_dat_clinvar_df <- target_dat_clinvar_df[,!c(names(target_dat_clinvar_df) %in% c("Chromosome", "Start", "Stop"))] %>% head

names(target_dat_clinvar_df) <- gsub("start.variant_caller","start", names(target_dat_clinvar_df))
names(target_dat_clinvar_df) <- gsub("end.variant_caller","end", names(target_dat_clinvar_df))
names(target_dat_clinvar_df) <- gsub("chr.variant_caller","chr", names(target_dat_clinvar_df))


joined_ranges <- makeGRangesFromDataFrame(data.frame(target_dat_clinvar_df), seqnames.field ="chr" , start.field='start', end.field='end', keep.extra.columns=TRUE)
if(rangeordf == "range"){
  save(joined_ranges,file=output_file)  
} else {
  fwrite(target_dat_clinvar_df, file=output_file)
}
