## analysis codes for snATAC-seq data

library('GenomicRanges')

## step 0. load related files
# snp
snp <- readRDS('step1.snp.df.rds')
snp <- snp[,1:5]
snp.gr <- GenomicRanges::GRanges(seqnames=snp$CHR,
                                 IRanges::IRanges(start=snp$POS,end=snp$POS))

# snATAC-seq data
mat <- Matrix::readMM('cCRE_by_cell_type_matrix.tsv') ## region by cell type 1154611
ct <- read.table('cCRE_by_cell_type_celltypes.txt',sep = ',')
ct = ct$V1
pks <- read.table('cCRE_by_cell_type_cCREs.bed')
pks.gr <- GenomicRanges::GRanges(seqnames=pks$V1,
                                 IRanges::IRanges(start=pks$V2,end=pks$V3))


## step 1. subset to get overlaps between snp and peaks
snp_ol <- findOverlaps(snp.gr, pks.gr)
snp_sub <- snp[as.data.frame(snp_ol)$queryHits,]
mat_sub <- matrix(as.integer(as.matrix(mat[as.data.frame(snp_ol)$subjectHits,])),
                  ncol = dim(mat)[2]) # 105321
head(rowSums(mat_sub))
pks_sub <- pks[as.data.frame(snp_ol)$subjectHits,]

## step 2. get scores -- where is posterior probability??

# 2.1 row normalization
mat_sub <- t(apply(mat_sub, 1, FUN = function(x) x/sum(x)))
head(rowSums(mat_sub))

# sum over all snps ## -- to be modified
tiss_sums <- colSums(mat_sub)
names(tiss_sums) <- ct
tiss_sums <- sort(tiss_sums, decreasing = T)



# > head(tiss_sums)
# Esophageal Epithelial                Glutamatergic 1 
# 2583.970                       2346.292 
# Fetal Syncitio+Cytotrophoblast                     Follicular 
# 2231.713                       1948.277 
# Glutamatergic 2         Fetal Adrenal Cortical 
# 1788.724                       1614.298 


