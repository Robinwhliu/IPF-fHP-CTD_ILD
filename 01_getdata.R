rm(list = ls())
options(stringsAsFactors = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
# options(BioC_mirror="https://anaconda.org/bioconda/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

getwd()
setwd('C:/Users/LWH/Downloads/bioinfo/IPF&CHD&SSc')

library(sva)
library(limma)
library(dplyr)
library(plyr)
library(tidyr)
library(R.utils)
library(biomaRt)
library(magrittr)
GSE_ID <- c('GSE48149')
destdir<-c("01_getdata")
GSE_file<-paste0(GSE_ID,"_eSet.Rdata") 
# gset <- downGSE(GSE_ID, destdir = destdir)
# save(gset, file = paste0("01_getdata/",GSE_file))
load(paste0("01_getdata/",GSE_file))
func1 = function(l, x) {
  h = c()
  for (i in c(1:l$lengths[x])) {
    p = paste0(l$values[x], i, collapse = "_")
    h = c(h,p)
  }
  return(h)
}

library("GEOquery")
exprSet <- gset[[1]]
assayData <- exprs(exprSet)
phenoData <- pData(exprSet)
colnames(phenoData)
pheno = phenoData[-grep("PAH", phenoData$description),]
exp = assayData[,rownames(pheno)]

pd = pheno[,c("title", "description")]
colnames(pd) = c("title", "Group")
gpl <- exprSet@annotation
featureData <- getGEO(gpl, destdir="01_getdata/")
fd = featureData@dataTable@table
fd = fd[!is.na(fd$Entrez_Gene_ID),]
index = intersect(fd$ID, rownames(assayData))
fd = fd[fd$ID %in% index,]
assayData = assayData[index,]
assayData = as.data.frame(assayData)
assayData$ENTREZID = as.character(fd$Entrez_Gene_ID)
assayData$SD = apply(assayData[,-ncol(assayData)], 1, sd)
assayData = assayData[order(assayData$ENTREZID, assayData$SD, decreasing = T),]
assayData = assayData[!duplicated(assayData$ENTREZID),]
assayData = assayData[,-c(54,55)]
fd = fd[fd$ID %in% rownames(assayData),]
ad_GSE48149 = assayData; pd_GSE48149 = pd; fd_GSE48149 = fd
save(ad_GSE48149, pd_GSE48149, fd_GSE48149, file = "GSE48149.Rdata")
#########################################################################
load(file = "GSE150910.Rdata");table(pd_GSE150910$Group)
phenoData <- read.table(file = "01_getdata/GSE150910_series_matrix.txt", 
                        skip = 30, sep = "\t", header = T)
assayData <- read.csv(file = "01_getdata/GSE150910_de-identified_chp_kallisto_count_summary.csv",
                      header = T)

p <- rownames(assayData) %>% 
  lapply(function(x) unlist(strsplit(x,"\\|"))[2]) %>% 
  lapply(function(x) unlist(strsplit(x,"\\."))[1]) %>% unlist()
assayData$X = p
assayData$SD = apply(X = assayData[,2:289], MARGIN = 1, sd)
assayData = assayData[order(assayData$X, assayData$SD, decreasing = T),]
assayData = assayData[!duplicated(assayData$X),]
rownames(assayData) = assayData$X
mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl")
# searchAttributes(mart, "entrez")
IDs = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
            filters    = "ensembl_gene_id",
            values     = assayData$X, 
            mart       = mart)
fd = IDs[!is.na(IDs$entrezgene_id) & IDs$hgnc_symbol != "",]
fd = fd[!duplicated(fd$ensembl_gene_id),]
ad = assayData[fd$ensembl_gene_id,] %>% .[,-c(1,290)]

l = rle(phenoData[duplicated(phenoData[,1]),]$X.Sample_title)
phenoData[duplicated(phenoData[,1]),]$X.Sample_title = unlist(lapply(c(1:length(l$values)), 
                                         function(x) func1(l, x)))
rownames(phenoData) = phenoData$X.Sample_title 
pd = phenoData %>% t() %>% as.data.frame() %>% 
  .[,c("!Sample_characteristics_ch11", "!Sample_characteristics_ch12","!Sample_characteristics_ch16",
       "!Sample_characteristics_ch113","!Sample_characteristics_ch114", "!Sample_characteristics_ch15",
       "!Sample_characteristics_ch14", "!Sample_characteristics_ch13")]  %>% .[-1,]
colnames(pd) = c("Sex", "Age", "Genotype", "Plate", "Institution", "Race", "Smoker", "Group")
pd1 = apply(pd, 1, function(l) lapply(l, function(x) unlist(strsplit(x, ": "))[2])) %>% unlist() %>% 
  matrix(nrow=288, byrow=T) %>% data.frame(.,stringsAsFactors=FALSE)
rownames(pd1) = rownames(pd)
colnames(pd1) = colnames(pd)                                      
pd_GSE150910 = pd1
ad_GSE150910 = ad[,rownames(pd_GSE150910)]
fd_GSE150910 = fd
rownames(fd_GSE150910) = fd_GSE150910$ensembl_gene_id
save(ad_GSE150910, pd_GSE150910, fd_GSE150910, file = "GSE150910.Rdata")
###############################################################################
load(file = "GSE175457.Rdata");table(pd_GSE175457$Group)
phenoData <- read.table(file = "01_getdata/GSE175457_series_matrix.txt", 
                        skip = 49, sep = "\t", header = T)
assayData <- read.table(file = "01_getdata/GSE175457_geo_sysbio_rnaseq_data.txt",
                      sep = "\t", header = T)
assayData = assayData %>% t() %>% as.data.frame()
IDs = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
            filters    = "ensembl_gene_id",
            values     = rownames(assayData) %>% 
              lapply(function(x) unlist(strsplit(x, ".", fixed = T))[1]) %>% as.character(), 
            mart       = mart)
fd = IDs[!is.na(IDs$entrezgene_id) & IDs$hgnc_symbol != "",]
fd = fd[!duplicated(fd$ensembl_gene_id),]
rownames(fd) = fd$ensembl_gene_id
assayData$ENSG = rownames(assayData) %>% 
  lapply(function(x) unlist(strsplit(x, ".", fixed = T))[1]) %>% as.character()
assayData$SD = apply(assayData, 1, sd)
assayData = assayData[order(assayData$ENSG, assayData$SD, decreasing = T),]
assayData = assayData[!duplicated(assayData$ENSG),]
rownames(assayData) = assayData$ENSG
ad = assayData[rownames(fd),]
l = rle(phenoData[duplicated(phenoData[,1]),]$X.Sample_title)
phenoData[duplicated(phenoData[,1]),]$X.Sample_title = unlist(lapply(c(1:length(l$values)), 
                                                                     function(x) func1(l, x)))
rownames(phenoData) = phenoData$X.Sample_title
pd = phenoData %>% t() %>% as.data.frame() %>% 
  .[,c("!Sample_characteristics_ch16", "!Sample_characteristics_ch15","!Sample_characteristics_ch18",
       "!Sample_characteristics_ch11", "!Sample_characteristics_ch17",
       "!Sample_characteristics_ch19", "!Sample_characteristics_ch12")]  %>% .[-1,]
pd1 = apply(pd, 1, function(l) lapply(l, function(x) unlist(strsplit(x, ": "))[2])) %>% unlist() %>% 
  matrix(nrow=422, byrow=T) %>% data.frame(.,stringsAsFactors=FALSE)
colnames(pd1) = c("Sex", "Age", "Ethnicity", "Institution", "Race", "Smoker", "Type")
rownames(pd1) = rownames(pd) %>% 
  lapply(function(x) unlist(strsplit(x, ".", fixed = T))[1]) %>% as.character()
pd1$Group = rownames(pd1) %>% 
  lapply(function(x) unlist(strsplit(x, "_", fixed = T))[1]) %>% as.character()
ad_GSE175457 = ad[,rownames(pd1)]
pd_GSE175457 = pd1
fd_GSE175457 = fd

save(ad_GSE175457, pd_GSE175457, fd_GSE175457, file = "GSE175457.Rdata")

########################################################################
phenoData <- read.table(file = "01_getdata/GSE199152_series_matrix.txt", 
                        skip = 27, sep = "\t", header = T)
assayData <- read.table(file = "01_getdata/GSE199152_GeneCount.tsv",
                        sep = "\t", header = T)
# CodingLength <- assayData$CodingLength
ad = assayData[,c(2,grep("RPKM", colnames(assayData)))] %>% as.data.frame()
ad$SD <- apply(ad[,2:28], 1, sd)
ad <- ad[order(ad$GeneID, ad$SD, decreasing = T),]
ad <- ad[!duplicated(ad$GeneID),]
rownames(ad) = ad$GeneID
colnames(ad) = lapply(colnames(ad), function(x) gsub("_RPKM", "", x)) %>% unlist()
IDs = getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
            filters    = "hgnc_symbol",
            values     = ad$GeneID, 
            mart       = mart)
IDs <- IDs[order(IDs$hgnc_symbol, IDs$entrezgene_id),]
IDs <- IDs[!duplicated(IDs$hgnc_symbol),]
ad <- ad[IDs$hgnc_symbol,]
rownames(IDs) <- IDs$hgnc_symbol
l = rle(phenoData[duplicated(phenoData[,1]),]$X.Sample_title)
phenoData[duplicated(phenoData[,1]),]$X.Sample_title = unlist(lapply(c(1:length(l$values)), 
                                                                     function(x) func1(l, x)))
rownames(phenoData) = phenoData$X.Sample_title
pd = phenoData %>% t() %>% as.data.frame() %>% 
  .[,c("!Sample_geo_accession", "!Sample_characteristics_ch11")]  %>% .[-1,]
colnames(pd) = c("ID", "Group")
pd$Group = lapply(pd$Group, function(x)unlist(strsplit(x, ": "))[2]) %>% unlist()
pd1 <- pd[-grep("IPF", pd$Group),]
ad_GSE199152 = ad[,rownames(pd1)]
pd_GSE199152 = pd1
fd_GSE199152 = IDs

save(ad_GSE199152, pd_GSE199152, fd_GSE199152, file = "GSE199152.Rdata")
