
rm(list = ls())
options(stringsAsFactors = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
# options(BioC_mirror="https://anaconda.org/bioconda/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

getwd()
setwd('C:/Users/LWH/Downloads/bioinfo/IPF&CHD&SSc')
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(fgsea)
library(customLayout)
library(msigdbr)
library(GOplot)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(dplyr)

load(file = "venndata.Rdata")
load(file = "C:/Users/LWH/Downloads/bioinfo/IPF&CHD&SSc/res.Rdata")
res_IPF3[grep("P53",rownames(res_IPF3),ignore.case = T),]
res_CHP[grep("P53",rownames(res_CHP),ignore.case = T),]
res_CTD[grep("P53",rownames(res_CTD),ignore.case = T),]
res_IPF3[grep("P63",rownames(res_IPF3),ignore.case = T),]
res_CHP[grep("P63",rownames(res_CHP),ignore.case = T),]
res_CTD[grep("P63",rownames(res_CTD),ignore.case = T),]

load(file = "GSE175457.Rdata")

p53 <- ad_GSE175457[fd_GSE175457[grep("^TP53$|^TP63$",fd_GSE175457$hgnc_symbol,
                                      ignore.case = T,perl = T),]$
                      ensembl_gene_id,] %>% mutate_all(as.numeric)
rownames(p53) = c("TP63","TP53")
p53["Group",] <- pd_GSE175457$Group %>% 
  lapply(.,function(x) ifelse(x=="case","IPF","Control"))
data_melted <- reshape2::melt(t(p53) %>% as.data.frame(), 
                              id.vars = "Group") 
data_melted$value <- as.numeric(data_melted$value)
# write.table(data_melted,file = "C:/Users/LWH/Desktop/p53.txt",quote = F,
#             row.names = F, col.names = T, sep = "\t")
theme_set(theme_minimal(base_family = "Arial"))
p <- ggplot(data_melted, aes(x = Group, y = value, color = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5, 
              size = 2) +
  stat_compare_means(aes(group=Group),label="p.signif",label.x=1.5,
                     size=5,hjust=0.8,vjust=.2)+
  scale_color_jco() +
  facet_grid(. ~ variable, scales = "free_x", space = "free_x") +
  xlab("") + ylab("") +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(size = 14, face = "bold.italic"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    legend.position = "none"
  )
print(p)

universe = intersect(rownames(res_CHP), rownames(res_IPF)) %>% 
  intersect(rownames(res_CTD))
universe = fd_GSE48149[fd_GSE48149$Symbol %in% universe,]$Entrez_Gene_ID

ego_up <- enrichGO(gene          = gene_up,
                universe      = as.character(universe),
                ont           = "ALL",
                OrgDb         = org.Hs.eg.db,
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
ego_down <- enrichGO(gene          = gene_down,
                   universe      = as.character(universe),
                   ont           = "ALL",
                   OrgDb         = org.Hs.eg.db,
                   pAdjustMethod = "fdr",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
# save(ego_up,ego_down,file = "egos.Rdata")
load(file = "egos.Rdata")
ego_up@result %>% .[.$ONTOLOGY =="BP",] %>% .[1:5,] %>%  .$Description
ego_down@result %>% .[.$ONTOLOGY =="BP",] %>% .[1:5,] %>%  .$Description
ego_up@result %>% .[.$ONTOLOGY =="CC",] %>% .[1:5,] %>%  .$Description
ego_down@result %>% .[.$ONTOLOGY =="CC",] %>% .[1:5,] %>%  .$Description
ego_up@result %>% .[.$ONTOLOGY =="MF",] %>% .[1:5,] %>%  .$Description
ego_down@result %>% .[.$ONTOLOGY =="MF",] %>% .[1:5,] %>%  .$Description
goplot1 = dotplot(ego_up, x="GeneRatio", showCategory =5,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free') +
  theme(legend.text=element_text(size=8,face = "bold"),
        legend.title=element_text(size=8,face = "bold"),
        axis.text.y=element_text(size=8,face = "bold",lineheight =0.75), 
        axis.text.x=element_text(size=8,face = "bold"),
        axis.title.x=element_text(size=8,face = "bold"),
        strip.text.y = element_text(size = 7,face = "bold")) +
  scale_color_continuous (name="FDR", low = "red", high ="blue", 
                          labels = function(col) format(col, scientific = TRUE))

goplot2 = dotplot(ego_down, x="GeneRatio", showCategory =5,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')  +
  theme(legend.text=element_text(size=8,face = "bold"),
        legend.title=element_text(size=8,face = "bold"),
        axis.text.y=element_text(size=8,face = "bold",lineheight =0.75), 
        axis.text.x=element_text(size=8,face = "bold"),
        axis.title.x=element_text(size=8,face = "bold"),
        strip.text.y = element_text(size = 7,face = "bold")) +
  scale_color_continuous (name="FDR", low = "red", high ="blue", 
                          labels = function(col) format(col, scientific = TRUE))

# library(ReactomePA)
# x <-enrichPathway(gene_up,organism = "human",
#                   pvalueCutoff=0.001, 
#                   pAdjustMethod = "fdr",
#                   qvalueCutoff = 0.01,
#                   readable=T)
# 
# Reactome = barplot(x, color = "pvalue",showCategory=10)+
#   theme(legend.text=element_text(size=8),legend.title=element_text(size=8,face = "bold"),
#         legend.key.size= unit(3, "mm") ,
#         legend.margin= unit(0.5,"mm"),
#         plot.margin = margin(t=0, r=0, b=0, l=1, unit="mm"),
#         axis.text.y=element_text(size=8), axis.text.x=element_text(size=8),
#         axis.title.x=element_text(size=8,face = "bold"))+ 
#   scale_fill_continuous(name="FDR",low = "red", high ="blue", 
#                         labels = function(col) format(col, scientific = TRUE))


kk_up <- enrichPathway(gene         = gene_up,
                       universe     = as.character(universe),
                       organism     = 'human',
                       pvalueCutoff = 0.05, qvalueCutoff = 0.05)
kk_up <- setReadable(kk_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kk_down <- enrichPathway(gene         = gene_down,
                         universe     = as.character(universe),
                         organism     = 'human',
                         pvalueCutoff = 0.05, qvalueCutoff = 0.05)
kk_down <- setReadable(kk_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kk_up@result[1:10,]$Description
kk_down@result[1:10,]$Description

get_KEGGplot<-function(kk){
  kEGGr = kk@result
  kEGGr$p.adjust = kEGGr$pvalue
  kk@result = kEGGr
  kk@result$qvalue = kEGGr$pvalue
  class(kk)
  
  KEGGplot <- barplot(kk, color = "pvalue",showCategory=10)+
    theme(legend.text=element_text(size=8,face = "bold"),
          legend.title=element_text(size=8,face = "bold"),
          axis.text.y=element_text(size=8,face = "bold"), 
          axis.text.x=element_text(size=8,face = "bold"),
          axis.title.x=element_text(size=8,face = "bold"))+
    scale_fill_continuous(name="FDR", low = "red", high ="blue", 
                          labels = function(col) format(col, scientific = TRUE))+
    # theme(legend.position = "none")+
    theme(plot.margin = margin(t=0, r=0, b=0, l=1, unit = "mm"),
          legend.key.size = unit(3, "mm"),
          legend.margin= unit(0.5,"mm"))
  return(KEGGplot)
}
KEGGplot1 <- get_KEGGplot(kk_up)
KEGGplot2 <- get_KEGGplot(kk_down)


library(patchwork)
patch <- ((goplot1/goplot2 + plot_layout(heights = c(3,2))) | (KEGGplot1/KEGGplot2))+
  plot_layout(widths = c(5,4))+
  # plot_annotation(tag_levels = 'A')&
  theme(plot.tag = element_text(size=9,face="bold"))

rownames(fd_GSE48149) = fd_GSE48149$Entrez_Gene_ID
gene_name1 <- fd_GSE48149[gene_up %>% as.character(),]$Symbol
genelist_DO1 = data.frame(ID = gene_up, 
                         logFC = res_CHP[gene_name1,]$FoldChange*
                           res_CTD[gene_name1,]$FoldChange*
                           res_IPF3[gene_name1,]$FoldChange)
gene_name2 <- fd_GSE48149[gene_down %>% as.character(),]$Symbol
genelist_DO2 = data.frame(ID = gene_down, 
                          logFC = res_CHP[gene_name2,]$FoldChange*
                            res_CTD[gene_name2,]$FoldChange*
                            res_IPF3[gene_name2,]$FoldChange)
get_DOplot <- function(genelist_DO,fd){
  gene<-genelist_DO$ID
  DO <- DOSE::enrichDO(
    gene = gene,
    pvalueCutoff = 0.05,
    universe = as.character(universe)
    # pAdjustMethod = "none"
  )
  # DO <- setReadable(DO, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  # t(DO@result[1:10,])
  DO_result = DO@result
  DO_result<-DO_result[DO_result$pvalue<0.05,]
  DO_result = DO_result[-grep("cancer|cinoma|neoplasm", DO_result$Description, 
                              ignore.case = T, perl = T),]
  DO@result = DO_result
  do=data.frame(Category = "DO",ID = DO_result$ID,
                Term = DO_result$Description %>% str_to_title(), 
                Genes = gsub("/", ", ", DO_result$geneID), adj_pval = DO_result$pvalue)
  circ <- circle_dat(do, genelist_DO)
  circ$genes = fd[as.character(circ$genes),]$Symbol
  rownames(genelist_DO) = genelist_DO$ID 
  genelist_DO$ID = fd[rownames(genelist_DO),]$Symbol
  termNum = 8	#闄愬畾GO鏁扮???
  termNum=ifelse(nrow(do)<termNum,nrow(do),termNum)
  geneNum = nrow(genelist_DO)	#闄愬畾鍩哄洜鏁扮???
  chord <- chord_dat(circ, genelist_DO[1:geneNum,], do$Term[1:termNum])
  chord[,"logFC"] = (chord[,"logFC"] %>% as.numeric() %>%  sign(.))*
    (chord[,"logFC"] %>% as.numeric() %>%  abs(.) %>% .^(1/2))
  if(chord[,"logFC"] %>% as.numeric() %>%  sign(.)>0) {
    DOplot= GOChord(chord, lfc.max = 5,lfc.min = 1,process.label = 10,
                    space = 0.001, 
                    gene.order = 'logFC', 
                    gene.space = 0.25,      
                    gene.size = 5,
                    lfc.col=c('firebrick3', 'white','royalblue3')) + 
      theme(legend.text = element_text(size = 13,face = "bold"),
            legend.position = "bottom",legend.direction = "vertical",
            legend.text.align = 0, legend.box.just = "top", 
            legend.title = element_text(size = 14, face = "bold")) + 
      scale_fill_continuous(breaks = seq(1,5,by=0.5), limits=c(1,5), 
                            low = "blue", high ="red")
  }
  else{
    DOplot= GOChord(chord, lfc.max = -1,lfc.min = -5,process.label = 10,
                    space = 0.001, 
                    gene.order = 'logFC', 
                    gene.space = 0.25,      
                    gene.size = 5,
                    lfc.col=c('firebrick3', 'white','royalblue3')) +  
      theme(legend.text = element_text(size = 13,face = "bold"),
            legend.position = "bottom",legend.direction = "vertical",
            legend.text.align = 0, legend.box.just = "top",
            legend.title = element_text(size = 14, face = "bold")) + 
      scale_fill_continuous(breaks = seq(-4,-1,by=0.5), limits=c(-4,-1),
                            low = "blue", high ="red")+
      scale_color_manual(guide = guide_legend(nrow = 3))
  }
  print(do$Term[1:termNum])
  return(DOplot)
}
get_DOplot(genelist_DO1,fd_GSE48149)
get_DOplot(genelist_DO2,fd_GSE48149)

##########################################################################
load(file = "egos.Rdata")
degs1 = res_IPF3[abs(res_IPF3$FoldChange) > 2.2 & res_IPF3$adjusted_P.value <.001,]
degs2 = res_CHP[abs(res_CHP$FoldChange) > 1 & res_CHP$adjusted_P.value <.001, ]
degs3 = res_CTD[abs(res_CTD$FoldChange) > 1 & res_CTD$adjusted_P.value <.05, ]

IPF_up = fd_GSE175457[fd_GSE175457$hgnc_symbol %in% 
                       rownames(degs1[degs1$FoldChange > 0,]),]$entrezgene_id
IPF_down = fd_GSE175457[fd_GSE175457$hgnc_symbol %in% 
                       rownames(degs1[degs1$FoldChange < 0,]),]$entrezgene_id
CHP_up = fd_GSE150910[fd_GSE150910$hgnc_symbol %in% 
                       rownames(degs2[degs2$FoldChange > 0,]),]$entrezgene_id
CHP_down = fd_GSE150910[fd_GSE150910$hgnc_symbol %in% 
                         rownames(degs2[degs2$FoldChange < 0,]),]$entrezgene_id
CTD_up = fd_GSE48149[fd_GSE48149$Symbol %in% 
                       rownames(degs3[degs3$FoldChange > 0,]),]$Entrez_Gene_ID
CTD_down = fd_GSE48149[fd_GSE48149$Symbol %in% 
                         rownames(degs3[degs3$FoldChange < 0,]),]$Entrez_Gene_ID

get_ego <- function(gene, universe, description) {
  ego <- enrichGO(gene          = gene %>% as.character(),
                  universe      = universe %>% as.character(),
                  ont           = "BP",
                  OrgDb         = org.Hs.eg.db,
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  ego@result <- ego@result %>% .[!.$Description %in% description,]
  return(ego)
}
ego_up_IPF <- get_ego(IPF_up, fd_GSE175457$entrezgene_id, ego_up@result$Description)
ego_down_IPF <- get_ego(IPF_down, fd_GSE175457$entrezgene_id, ego_down@result$Description)
ego_up_CHP <- get_ego(CHP_up, fd_GSE150910$entrezgene_id, ego_up@result$Description)
ego_down_CHP<- get_ego(CHP_down, fd_GSE150910$entrezgene_id, ego_down@result$Description)
ego_up_CTD <- get_ego(CTD_up, fd_GSE48149$Entrez_Gene_ID, ego_up@result$Description)
ego_down_CTD<- get_ego(CTD_down, fd_GSE48149$Entrez_Gene_ID, ego_down@result$Description)


# get_goplot <- function(ego) {
#   goplot = dotplot(ego, x="GeneRatio", showCategory =5,split="ONTOLOGY") + 
#     facet_grid(ONTOLOGY~., scale='free') +
#     theme(legend.text=element_text(size=12,face = "bold"),
#           legend.title=element_text(size=12,face = "bold"),
#           axis.text.y=element_text(size=14,lineheight = 0.45,angle=0), 
#           axis.text.x=element_text(size=14),
#           axis.title.x=element_text(size=14,face = "bold"),
#           strip.text.y = element_text(size = 14,face = "bold")) +
#     scale_color_continuous (name="FDR", low = "red", high ="blue", 
#                             labels = function(col) format(col, scientific = TRUE)) + 
#     scale_y_discrete(labels=function(x) str_wrap(x, width=40),expand = c(0,0.6))
#   return(goplot)
# }

get_goplot <- function(ego) {
  goplot = dotplot(ego, x="GeneRatio", showCategory =10) + 
    theme(legend.text=element_text(size=18,face = "bold",family="serif"),
          legend.title=element_text(size=20,face = "bold",family="serif"),
          axis.text.y=element_text(size=20,lineheight = 0.6,face = "bold",family="serif"), 
          axis.text.x=element_text(size=20,face = "bold",family="serif"),
          axis.title.x=element_text(size=20,face = "bold",family="serif")) +
    scale_color_continuous (name="FDR", low = "red", high ="blue", 
                            labels = function(col) format(col, scientific = TRUE)) + 
    scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=35) %>% 
                       str_to_title(), expand = c(0,0.6))
  return(goplot)
}

# save(ego_up_IPF,ego_down_IPF,ego_up_CHP,ego_down_CHP,ego_up_CTD,ego_down_CTD,
#      file = "egos2.Rdata")
load(file = "egos2.Rdata")
# ego_up_CTD@result %>% .[.$ONTOLOGY == "BP",] %>% .$Description
ego_up_CTD@result$Description <- ego_up_CTD@result$Description %>% 
  replace(.,.=="chemokine-mediated signaling pathway","chemokine signaling pathway")
library(stringr)
goplot_up_IPF <- get_goplot(ego_up_IPF)
goplot_down_IPF <- get_goplot(ego_down_IPF)
goplot_up_CHP <- get_goplot(ego_up_CHP)
goplot_down_CHP <- get_goplot(ego_down_CHP)
goplot_up_CTD <- get_goplot(ego_up_CTD)
goplot_down_CTD <- get_goplot(ego_down_CTD)
goplot_up_CTD/goplot_down_CTD
# vars_go <- ls(pattern = "^goplot_")
save(goplot_up_IPF,goplot_down_IPF,goplot_up_CHP,
     goplot_down_CHP,goplot_up_CTD,goplot_down_CTD,file = "vars_go.Rdata")
