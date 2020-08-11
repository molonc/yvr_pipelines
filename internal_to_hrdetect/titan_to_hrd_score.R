#!/usr/bin/env Rscript
.libPaths("/ssd/sda1/sbeatty/software/miniconda3/lib/R/library")
require("readr", quietly=TRUE)
require("stringr", quietly=TRUE)
require("data.table", quietly=TRUE)
require("dplyr", quietly=TRUE)
require("signature.tools.lib")

args = commandArgs(trailingOnly=TRUE)



if(length(args) == 0){
args <- "Rscript /scratch/shahlab_tmp/sbeatty/yvr_pipelines/internal_to_hrdetect/titan_to_hrd_score.R -target_file=/scratch/shahlab_tmp/danlai/APARICIO-590/PAIR_TITA/titan_pipeline/OUTPUT/RUN/SA1259T_titan/outputs/results/TASK_11_CALC_OPTIMAL_CLUSTERS_SA1259T_titan_optimal_clusters.txt -outputfile=input_files/titan/SA1259T.hrd_score.csv"
args <- str_split(args," ") %>% unlist
}

titan_segs_path <- args[str_detect(args,"-target_file")] %>% str_replace("-target_file=","")
output_file <- args[str_detect(args,"-outputfile")] %>% str_replace("-outputfile=","")
#sample_sheet <- fread("/scratch/shahlab_tmp/sbeatty/ind231/hrdetect/metadata/HRDetect_sample_sheet.csv")
#titan_segs_path <- sample_sheet$titan[1]
#segs_to_HRD_LOD <- function(titan_segs_path){

	get_optimal_clusters <- function(path_i){
		optimal_file <- suppressWarnings(suppressMessages(read_delim(path_i, ": ", col_names=FALSE))) %>% data.frame
		optimal_results <- optimal_file[str_detect(optimal_file[,1], "optimal Clusters"),2] %>% trimws()
		optimal_results <- str_replace_all(optimal_results, pattern="[^.0-9]", replacement=" ") %>% trimws() %>% str_split("  ") %>% unlist %>% trimws
		optimal_ploidy <- optimal_results[2]  %>% as.numeric
		optimal_clusters <- optimal_results[1] %>% as.numeric
		titan_files_sample_i <- list.files(dirname(path_i))
		optimal_seg_file_path <- paste(
			dirname(path_i),
			titan_files_sample_i[
			str_detect(titan_files_sample_i, "TASK_10_ANNOT_PYGENE") & 
			str_detect(titan_files_sample_i, paste0("clusters_", floor(optimal_clusters)))
			][1], sep="/"
			)
		output_df <- fread(optimal_seg_file_path, data.table=FALSE)
		output_df
		}


df <- tryCatch({get_optimal_clusters(titan_segs_path)}, error=function(e){return(NA)})

if(length(df) == 1){
	fwrite(data.frame(df), file=paste(output_file))
	q()
}
ascat.data2 <- data.frame(SampleID=df$Sample, 
	Chromosome=df$Chromosome, 
	Start=df[,"Start_Position(bp)"], 
	End=df[,"End_Position(bp)"],
	nProbes=rep(NA,nrow(df)),
	totalCN=df$Copy_Number,
	nA=df$Copy_Number - df$MinorCN,
	nB=df$MinorCN,
	Ploidy=rep(NA,nrow(df)),
	AberrantCellFraction=rep(NA ,nrow(df)) 
	)
ll<-match(c("SampleID","Chromosome","Start","End","nProbes","totalCN","nA","nB","Ploidy" ,"AberrantCellFraction"),colnames(ascat.data2))
ascat.data2 <-ascat.data2 [,ll]
rm(ll)

ascat.data2[,"Chromosome"] <- as.character(ascat.data2[,"Chromosome"])
ascat.data2[ascat.data2[,"Chromosome"]=="X","Chromosome"] <- "23"
ascat.data2[ascat.data2[,"Chromosome"]=="Y","Chromosome"] <- "24"
ascat.data2[,"Chromosome"] <- as.numeric(ascat.data2[,"Chromosome"])
HRD_LOH <- calc.hrd(ascat.data2, nA=7,check.names=FALSE, return.loc=FALSE)
fwrite(data.frame(HRD_LOH), file=paste(output_file))

