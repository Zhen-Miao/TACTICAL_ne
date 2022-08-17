## analysis codes for snATAC-seq data
## version 3
## author: Zhen Miao

#' ---
#' Input snATAC-seq data: binary 
#' Input GARFIELD: no
#' Input posterior probability: yes
#' Filter cell types: no
#' Combine cell types with Cell ontology: yes
#' With orthogonal annotations: no
#' Peak filtering: no
#' Forced Sparsity: by rank of column sum
# ---
  
  
  
  


library('GenomicRanges')

## step 0. load related files
# 0.1 snp
snp <- readRDS('step1.snp.df.rds')
snp <- snp[,1:5] ## the VALUE column is the posterior probability
snp.gr <- GenomicRanges::GRanges(seqnames=snp$CHR,
                                 IRanges::IRanges(start=snp$POS,end=snp$POS))

# 0.2 snATAC-seq data
mat <- Matrix::readMM('cCRE_by_cell_type_matrix.tsv') ## region by cell type 1154611
mat <- as.matrix(mat)
pks <- read.table('cCRE_by_cell_type_cCREs.bed')
pks.gr <- GenomicRanges::GRanges(seqnames=pks$V1,
                                 IRanges::IRanges(start=pks$V2,end=pks$V3))

## 0.3 Cell ontology
ct <- read.table('cCRE_by_cell_type_celltypes.txt',sep = ',')
ct = ct$V1
co <- read.csv('Cell_ontology.tsv',sep = '\t')
length(ct)
sum(ct %in% co$Cell.type)  ## 69
# ct[!(ct %in% co$Cell.type)]

co2 <- gsub('\\(|\\)','',co$Cell.type)
ct2 <- gsub('\\(|\\)','',ct)
co2 <- gsub('\\+|\\-|\\,',' ',co2)
ct2 <- gsub('\\+|\\-|\\,',' ',ct2)
co2 <- sub('Fibroblast','Fibro',co2) 
ct2 <- sub('Fibroblast','Fibro',ct2) 
co2 <- sub(' Cell','',co2) 
ct2 <- sub(' Cell','', ct2)
ct2 <- sub(' cell','', ct2)
co2 <- sub(' cell','', co2)
co2 <- sub('Smooth','Sm',co2)
ct2 <- sub('Smooth','Sm',ct2)
co2 <- sub('Muscle','Ms',co2)
ct2 <- sub('Muscle','Ms',ct2)
co2 <- sub('Pancreatic ','',co2)
co2 <- sub('Alverolar','Alveolar', co2)
co2 <- gsub(' AT1','',co2)
co2 <- gsub(' AT2','',co2)
co2 <- sub('GABAergic Neuron','GABA',co2)
co2 <- sub('Glutamatergic Neuron','Glutamatergic',co2)
co2 <- sub('CNS Enteric','Enteric',co2)
co2 <- sub('Blood Brain Barrier ','BBB ',co2)

co2 <- sub('Epithelial', 'Epi',co2)
ct2 <- sub('Epithelial', 'Epi',ct2)

sum(ct2 %in% co2)  ## 222


# 0.4 cell type annotations
#### ------------------ Naively, remove beta_2, only keep beta_1 -----------
colnames(mat) <- ct2
co$Cell.type <- co2

saveRDS(co,'cell_ontology_matched.rds')
saveRDS(ct2, 'cell_type_matched.rds')

co_new <- unique(co$closest.Cell.Ontology.term.s.)
mat_new <- matrix(nrow = nrow(mat), ncol = length(co_new))
colnames(mat_new) <- co_new

for(coi in co_new){
  ct_sel <- co[co$closest.Cell.Ontology.term.s. == coi,'Cell.type']
  if(length(ct_sel) > 1){
    mat_new[,coi] <- as.numeric(rowSums(mat[,ct_sel]) > 0)
  }else{
    mat_new[,coi] <- mat[,ct_sel]
  }
}


# quantile(colSums(mat))
# # 0%    25%    50%    75%   100% 
# # 7663  52839  68388  91420 162264 
quantile(colSums(mat_new))
# 0%    25%    50%    75%   100% 
# 33493  58988  92437 116465 262791

#### replace mat by mat_new
mat <- mat_new


#### ------------------ try normalize the matrix by total read count
mat <- apply(mat, 2, FUN = function(x) x/sum(x) )
mat <- mat * 92437


## step 1. subset to get overlaps between snp and peaks
snp_ol <- findOverlaps(snp.gr, pks.gr)
snp_sub <- snp[as.data.frame(snp_ol)$queryHits,]
mat_sub <- mat[as.data.frame(snp_ol)$subjectHits,]
head(rowSums(mat_sub))
pks_sub <- pks[as.data.frame(snp_ol)$subjectHits,]

## step 2. get scores -- the value column is posterior probability

# 2.1 get multiply by ppi
mat_ppi <- mat_sub * snp_sub$VALUE


# 2.2 upweight by column sum -- decided to take out this step
# ##### -- define enrichment by relative score compared with the median
# csum <- colSums(mat_ppi)
# rel_enrich <- csum / median(csum)
# quantile(rel_enrich) ## -- make sure the relative enrichment fall into reasonable range
# # 0%        25%        50%        75%       100% 
# # 0.07556083 0.77694464 1.00000000 1.30625885 2.11875622
# mat_ppi <- t(t(mat_ppi) * rel_enrich)

# 2.3 row normalization
mat_ppi <- t(apply(mat_ppi, 1, FUN = function(x) x/ (sum(x)+0.01) ))
# head(rowSums(mat_ppi))

# 2.4 sum over all snps 
tiss_sums <- colSums(mat_ppi)
names(tiss_sums) <- (colnames(mat))
tiss_sums <- sort(tiss_sums, decreasing = T)

## display top hits
tiss_sums[1:20]

## results based on cell ontologies 

excitatory neuron 
43.00364 
macrophage 
36.29790 
glutamatergic neuron 
35.04337 
trophoblast giant cell 
34.67637 
enterocyte 
34.21110 
non keratinizing barrier epithelial cell 
31.42845 
extravillous trophoblast 
30.60441 
fibroblast 
29.40888 
endothelial cell 
29.23571 
erythroblast 
29.07306 
fetal thymocyte 
28.90085 
skeletal muscle fiber 
27.66164 
alveolar macrophage 
27.60582 
inhibitory neuron 
27.17573 
hepatoblast 
26.04913 
cortical cell of adrenal gland 
25.98347 
B cell 
25.56557 
colon epithelial cell 
24.69133 
fetal cardiomyocyte 
24.21741 
ventricular cardiac muscle cell 
24.10320 
