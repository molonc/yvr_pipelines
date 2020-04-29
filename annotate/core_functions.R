
genome_range_to_dataframe <- function(genome_range){
  df <- genome_range
  df <- target_dat
  df_rowRanges <- rowRanges(df)
  range_start <- ranges(df)@start
  range_width <- ranges(df)@width
  range_end <- c(range_start + range_width) - 1
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

invert_snv <- function(x, nucleotides_positive_strand=c("a", "g", "c", "t", "A", "G", "C", "T"), nucleotides_positive_negative=c("t","c", "g", "ad", "T", "C", "G", "A")){
  if(x %in% nucleotides_positive_strand){
    output <- nucleotides_out[nucleotides_in %in% x]
    } else {
      output <- "nucleotide inversion error"
      }
      output 
    }




extract_exact_matches <- function(matches_df, query_df,query_ref_allele_column,query_alt_allele_column,  subject_df, subject_ref_allele_column, subject_alt_allele_column, query_strand_column=NA, subject_strand_column=NA){
.invert_snv <- function(x, nucleotides_positive_strand=c("a", "g", "c", "t", "A", "G", "C", "T"), nucleotides_positive_negative=c("t","c", "g", "ad", "T", "C", "G", "A")){
  if(x %in% nucleotides_positive_strand){
    output <- nucleotides_out[nucleotides_in %in% x]
    } else {
      output <- "nucleotide inversion error"
      }
      output 
    }
  .invert_nucleotide <- function(x){
  x <- "AGCT"
  x_split <-  str_split("act","") %>% unlist
  paste(mapply(x_split, FUN=invert_snv, USE.NAMES=FALSE), collapse="")
}

matches_df[,"gencode_strand"] <- query_df[matches_df[,"queryHits"],"genecode_strand"] 
matches_df[,"query_ref"] <- query_df[matches_df[,"queryHits"],paste(query_ref_allele_column)] 
matches_df[,"query_alt"] <- query_df[matches_df[,"queryHits"],paste(query_alt_allele_column)] 
matches_df[,"subject_ref"] <- subject_df[matches_df[,"subjectHits"],paste(subject_ref_allele_column)] 
matches_df[,"subject_alt"] <- subject_df[matches_df[,"subjectHits"],paste(subject_alt_allele_column)] 


if(!is.na(query_strand_column) & !is.na(subject_strand_column){
  matches_df[, "subject_alt_inverted"] <- mapply(matches_df$subject_alt, FUN=.invert_nucleotide)
}


# match reference alleles 
reference_variant_matches <- c(matches_df$query_ref ==  matches_df$subject_ref)
# eliminate matches where the reference alleles vary
matches_df <- matches_df[reference_variant_matches,]

# exact match snv
snv_exact_matches <- c(matches_df$query_alt == matches_df$subject_alt)

# length match indel
indel_length_matches <- c(nchar(matches_df$query_alt) == nchar(matches_df$subject_alt)) & c(nchar(matches_df$query_alt) > 1)

matches_df <- matches_df[c(snv_exact_matches | indel_length_matches), ]
matches_df
}

substitute_column_name <- function(in_name,out_name, input_data){
  target_data_column_names <- names(input_data)
  matches <- which(target_data_column_names == in_name)
  if(length(matches) == 1){
    target_data_column_names[matches] <- out_name
  }
  target_data_column_names
}






