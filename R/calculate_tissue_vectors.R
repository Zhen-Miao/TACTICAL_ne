#' Calculate a tissue score vector for each SNP within a genetic credible set or for each index SNP
#'
#' This function takes the SNP annotations outputted from the annotate_snps function and the
#' provided annotation weigths to yield a set of tissue scores that indicate how much of the
#' fine-mapping/association signal can be apportioned to each tissue/cell type
#'
#' @param snp.annotated.df Dataframe of annotated SNP info, this is output from annotate_snps function
#' @param tissue_annotation_file Path to the input file of tissue annotation names and weights
#' @param genomic_annotation_file Path to the input file of genomic annotation names and weights]
#' @param ess.annot Optional: Name of genomic annotation requiring specificity scores (e.g. "coding")
#' @param ess.file Optional: Path to the specificity score file for ess.annot, (e.g. expression specificity scores for coding annotations)
#' @export
calculate_tissue_vectors <- function(snp.annotated.df,tissue_annotation_file,genomic_annotation_file,
                                     ess.annot=NULL,ess.file=NULL){

  ## require some functions 
  require('data.table')
  
  ## read files 
  tiss.annot.df <- read.table(tissue_annotation_file)
  # gen.annot.df <- data.table::fread(genomic_annotation_file)
  
  ## define some functions
  "%&%" <- function(a,b) paste0(a,b) # just a shortcut for the paste function
  '%>%' <- magrittr::'%>%'
  
  ## pre-order data frames
  rr = 1 
  full.row.df <- snp.annotated.df[rr,]
  row.df <- full.row.df %>%
    dplyr::select(.,-one_of("SIGNAL","SNPID","CHR","POS","VALUE"))
  
  ## get all the tissue types
  name.vec <- names(row.df)
  tiss.vec <- purrr::map(name.vec,function(s){
    vec <- strsplit(x=s,split=".",fixed=T)[[1]]
    ifelse(length(vec)==2,vec[1],NA)
  }) %>% unlist(.) %>% na.omit(.) %>% as.character(.) %>% unique(.)  %>% sort(.)
  
  ## get all the annotation types
  tiss.annot.vec <- purrr::map(name.vec,function(s){
    vec <- strsplit(x=s,split=".",fixed=T)[[1]]
    ifelse(length(vec)==2,vec[2],NA)
  }) %>% unlist(.) %>% na.omit(.) %>% as.character(.) %>% unique(.) %>% sort(.)
  
  ## get all the genome annotation types (only one type, coding)
  gen.annot.vec <- purrr::map(name.vec,function(s){
    vec <- strsplit(x=s,split=".",fixed=T)[[1]]
    ifelse(length(vec)==1,vec[1],NA)
  }) %>% unlist(.) %>% na.omit(.) %>% as.character(.) %>% unique(.) %>% sort(.)
  
  ## step 0. remove the coding part from the data frame 
  snp.annotated.df <- snp.annotated.df[,coding:=NULL]
  
  ## step 1. reorder the annotations 
  ta = vector(length = length(tiss.annot.vec) * length(tiss.vec)  )
  for (i in 1:length(tiss.annot.vec)){
    annot_wo_e <- tiss.annot.vec[i]
    print.warning <- TRUE
    for (e in 1:length(tiss.vec)){
      tiss <- tiss.vec[e]
      ta[(i-1)*length(tiss.vec) + e]  = tiss%&%"."%&%annot_wo_e
    }
  }
  
  ta_snp <-  append( c("SIGNAL","SNPID","CHR","POS","VALUE"),ta)
  setcolorder(snp.annotated.df, ta_snp)
  
  rownames(tiss.annot.df) = tiss.annot.df$V1%&%"."%&%tiss.annot.df$V2
  tiss.annot.df = tiss.annot.df[ta,]
  
  ## step 1.5 calculate some values
  n_tissues <- length(tiss.vec)
  n_snp <- dim(snp.annotated.df)[1]
  n_annos <- dim(snp.annotated.df)[2] ## c("SIGNAL","SNPID","CHR","POS","VALUE")
  
  ### step 2. restructure the matrices by snp 
  
  ## ---- this matrix is only constructed once -------
  annot.matrix_all <- matrix(as.numeric(tiss.annot.df$V3), nrow = n_tissues)
  # rownames(annot.matrix_all) <- tiss.vec
  # colnames(annot.matrix_all) <- tiss.annot.vec
  
  ## use tiss.wgts matrix to store all the tissue weights 
  tiss.wgts <- matrix(nrow = n_snp,ncol = n_tissues)
  colnames(tiss.wgts) <- tiss.vec
  snp_trunc_df <- as.matrix(as.data.frame(snp.annotated.df)[,6:n_annos])
  # snp_trunc_df <- as.numeric(snp_trunc_df)
  
  ## loop for each SNP, note, no nested loop
    for (r in 1:n_snp){
    wgt <- matrix(snp_trunc_df[r,], nrow = n_tissues)
    annot.matrix <- annot.matrix_all * wgt
    
      if (sum(annot.matrix) > 0){
        tiss.wgts[r,] <- rowSums(annot.matrix) / sum(annot.matrix)
      } else{
        tiss.wgts[r,] <- rep(0,n_tissues)
      }
    }
  
  ## to be consistent to the original output 
  out.df <- cbind(as.data.frame(snp.annotated.df)[,1:5], as.data.frame(tiss.wgts))
  return(out.df) 

}

