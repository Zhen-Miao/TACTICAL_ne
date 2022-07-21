## analysis codes for snATAC-seq data
## version 1
## author: Zhen Miao

#' ---
#' Input snATAC-seq data: binary 
#' Input GARFIELD: no
#' Input posterior probability: yes
#' Filter cell types: no
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
ct <- read.table('cCRE_by_cell_type_celltypes.txt',sep = ',')
ct = ct$V1
pks <- read.table('cCRE_by_cell_type_cCREs.bed')
pks.gr <- GenomicRanges::GRanges(seqnames=pks$V1,
                                 IRanges::IRanges(start=pks$V2,end=pks$V3))

# 0.3 cell type annotations
#### ------------------ Naively, remove beta_2, only keep beta_1 -----------
colnames(mat) <- ct
mat <- mat[,!(colnames(mat) %in% 'Beta 1')]

quantile(colSums(mat))
# 0%    25%    50%    75%   100% 
# 7663  52839  68388  91420 162264 

#### ------------------ A better way is to filter by cell ontology


#### ------------------ try normalize the matrix by total read count
mat <- apply(mat, 2, FUN = function(x) x/sum(x) )
mat <- mat * 68388


## step 1. subset to get overlaps between snp and peaks
snp_ol <- findOverlaps(snp.gr, pks.gr)
snp_sub <- snp[as.data.frame(snp_ol)$queryHits,]
mat_sub <- mat[as.data.frame(snp_ol)$subjectHits,]
head(rowSums(mat_sub))
pks_sub <- pks[as.data.frame(snp_ol)$subjectHits,]

## step 2. get scores -- the value column is posterior probability

# 2.1 get multiply by ppi
mat_ppi <- mat_sub * snp_sub$VALUE


# 2.2 upweight by column sum
##### -- define enrichment by relative score compared with the median
csum <- colSums(mat_ppi)
rel_enrich <- csum / median(csum)
quantile(rel_enrich) ## -- make sure the relative enrichment fall into reasonable range
# 0%        25%        50%        75%       100% 
# 0.07556083 0.77694464 1.00000000 1.30625885 2.11875622
mat_ppi <- t(t(mat_ppi) * rel_enrich)

# 2.3 row normalization
mat_ppi <- t(apply(mat_ppi, 1, FUN = function(x) x/ (sum(x)+0.01) ))
# head(rowSums(mat_ppi))

# 2.4 sum over all snps 
tiss_sums <- colSums(mat_ppi)
names(tiss_sums) <- (colnames(mat))
tiss_sums <- sort(tiss_sums, decreasing = T)

## display top hits
tiss_sums[1:20]


# Esophageal Epithelial                     Follicular 
# 217.06980                      187.56700 
# Glutamatergic 1 Fetal Syncitio+Cytotrophoblast 
# 177.88691                      161.27068 
# Fetal Adrenal Cortical          Fetal V Cardiomyocyte 
# 154.27639                      140.30255 
# Fetal Enterocyte 1                         Beta 2 
# 139.01065                      132.78088 
# V Cardiomyocyte     Fetal Retinal Progenitor 1 
# 125.09076                      124.68867 
# Hepatocyte    Fetal Alveolar Epithelial 1 
# 118.25818                      116.13255 
# Fetal Hepatoblast                Glutamatergic 2 
# 112.88336                      108.38275 
# Alpha 1                         Beta 1 
# 106.98972                      106.23982 
# Fetal Adrenal Neuron                  Fibro General 
# 101.77440                      100.07019 
# Macrophage Gen or Alv              Fetal Metanephric 
# 98.12939                       93.48891



# Glutamatergic 1          Esophageal Epithelial 
# 146.58992                      143.65913 
# Fetal Syncitio+Cytotrophoblast                     Follicular 
# 130.38892                      129.60176 
# Fetal Adrenal Cortical                Glutamatergic 2 
# 119.50371                      113.76717 
# Hepatocyte                V Cardiomyocyte 
# 105.12881                      104.38472 
# Beta 2             Fetal Enterocyte 1 
# 99.06325                       96.67803 
# Fetal Retinal Progenitor 1          Fetal V Cardiomyocyte 
# 96.04882                       93.33094 
# Fetal Extravillous Trophoblast                        Alpha 1 
# 86.62172                       86.05304 
# Fetal Adrenal Neuron               Fetal Chromaffin 
# 85.39431                       84.64344 
# Fetal Hepatoblast          Macrophage Gen or Alv 
# 84.18361                       83.85708 
# A Cardiomyocyte                         Beta 1 
# 83.56684                       82.44296


### after remove beta 1
# Esophageal Epithelial                     Follicular 
# 217.86776                      188.42937 
# Glutamatergic 1 Fetal Syncitio+Cytotrophoblast 
# 178.41971                      161.57264 
# Fetal Adrenal Cortical          Fetal V Cardiomyocyte 
# 154.71120                      140.86673 
# Beta 2             Fetal Enterocyte 1 
# 139.45862                      139.43108 
# Fetal Retinal Progenitor 1                V Cardiomyocyte 
# 126.15310                      125.69496 
# Hepatocyte    Fetal Alveolar Epithelial 1 
# 118.91934                      116.73014 
# Fetal Hepatoblast                        Alpha 1 
# 113.41121                      110.66787 
# Glutamatergic 2           Fetal Adrenal Neuron 
# 108.88542                      102.19937 
# Fibro General          Macrophage Gen or Alv 
# 100.43103                       98.29606 
# Fetal Metanephric                A Cardiomyocyte 
# 93.87285                       93.53321 

#### after adjust for seq depth 
# Glutamatergic 1                     Follicular 
# 81.21125                       78.26656 
# Glutamatergic 2 Fetal Syncitio+Cytotrophoblast 
# 69.51780                       69.26189 
# Fetal Extravillous Trophoblast                         Beta 2 
# 67.88439                       65.57321 
# Fetal Enterocyte 1                     Hepatocyte 
# 65.14149                       65.10245 
# Fetal Adrenal Cortical          Esophageal Epithelial 
# 61.44501                       61.16026 
# Tuft                     Enterocyte 
# 59.88739                       57.65650 
# V Cardiomyocyte      Fetal Inhibitory Neuron 4 
# 56.05960                       54.85920 
# Fetal Hepatoblast                        Alpha 2 
# 54.41039                       52.82449 
# Fetal Photoreceptor          Fetal V Cardiomyocyte 
# 52.80879                       51.67985 
# Fetal Excitatory Neuron 3               Fetal Chromaffin 
# 51.29994                       51.11765 




