
genome_range_to_dataframe <- function(genome_range, extract_reads=FALSE){
  print(extract_reads)
  df <- genome_range
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
  if(extract_reads == FALSE){
  df_out <- data.frame(chr=df_seqnames, pos=range_start, start = range_start,
    end=range_end, width=range_width, strand=strand_info,
    paramRangeID=df_paramRangeID, REF=df_REF, ALT=df_ALT, QUAL=df_QUAL, Filter=df_FILTER)
  } else {
    geno_AD <- geno(df)$AD
    reference_reads <-  mapply(geno_AD,FUN=function(x){unlist(x)[1]})
    alternative_reads <-  mapply(geno_AD,FUN=function(x){unlist(x)[2]})
    df_out <- data.frame(chr=df_seqnames, pos=range_start, start = range_start,
    end=range_end, width=range_width, strand=strand_info,
    paramRangeID=df_paramRangeID, REF=df_REF, ALT=df_ALT, reference_reads=reference_reads, alternative_reads=alternative_reads,QUAL=df_QUAL, Filter=df_FILTER)
  }
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
#matches_df <- clin_var_matches
#query_df <- target_dat_df
#query_ref_allele_column <- target_ref_column
#query_alt_allele_column <- target_alt_column 
#subject_df<- clinvar_reference_dat
#subject_ref_allele_column <- "ReferenceAllele"
#subject_alt_allele_column <- "AlternateAllele"
#subject_strand_column <-NA
#query_strand_column <- NA



.invert_snv <- function(x, nucleotides_positive_strand=c("a", "g", "c", "t", "A", "G", "C", "T"), nucleotides_negative_strand=c("t","c", "g", "ad", "T", "C", "G", "A")){
  if(x %in% nucleotides_positive_strand){
    output <- nucleotides_negative_strand[nucleotides_positive_strand %in% x]
  } else {
    output <- "nucleotide inversion error"
    }
      output 
    }

.invert_nucleotide <- function(x){
  x_split <-  str_split(paste(x),"") %>% unlist
  paste(mapply(x_split, FUN=.invert_snv, USE.NAMES=FALSE), collapse="")
}

matches_df[,"query_ref"] <- query_df[matches_df[,"queryHits"],paste(query_ref_allele_column)] 
matches_df[,"query_alt"] <- query_df[matches_df[,"queryHits"],paste(query_alt_allele_column)] 
matches_df[,"subject_ref"] <- subject_df[matches_df[,"subjectHits"],paste(subject_ref_allele_column)] 
matches_df[,"subject_alt"] <- subject_df[matches_df[,"subjectHits"],paste(subject_alt_allele_column)] 
matches_df[, "subject_alt_inverted"] <- mapply(matches_df$subject_alt, FUN=.invert_nucleotide)
matches_df[, "subject_ref_inverted"] <- mapply(matches_df$subject_ref, FUN=.invert_nucleotide)

if(!is.na(query_strand_column) & !is.na(subject_strand_column)){
  matches_df[,"query_stand"] <- query_df[matches_df[,"queryHits"],paste(query_strand_column)] 
  matches_df[,"subject_stand"] <- subject_df[matches_df[,"subjectHits"],paste(subject_strand_column)] 

# reference allele match, alternate allele match, strandmatch
  match_all_of_ref_stand_alt_allele <- c(matches_df$query_ref ==  matches_df$subject_ref) & c(matches_df$query_stand ==  matches_df$subject_stand) & c(matches_df$query_alt ==  matches_df$subject_alt) 

# all matches when strands are corrected 
  match_all_of_ref_stand_alt_allele_after_inversion <- c(matches_df$query_stand !=  matches_df$subject_stand) & c(matches_df$query_ref ==  matches_df$subject_ref_inverted) & c(matches_df$query_alt ==  matches_df$subject_alt_inverted)

  match_vector <- c(match_all_of_ref_stand_alt_allele | match_all_of_ref_stand_alt_allele_after_inversion)

} else if(!is.na(subject_ref_allele_column) & !is.na(subject_alt_allele_column)){

# reference allele match, alternate allele match,
  match_vector <- c(matches_df$query_ref ==  matches_df$subject_ref) & c(matches_df$query_alt ==  matches_df$subject_alt)
}

# indel matches 
#indel_match_vector <- c(c(nchar(matches_df$query_alt) > 1) | c(nchar(matches_df$query_ref) > 1) | c(nchar(matches_df$subject_alt) > 1) | c(nchar(matches_df$subject_ref) > 1))


#matches_df <- matches_df[match_vector | indel_match_vector,]
matches_df <- matches_df[match_vector,]
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






