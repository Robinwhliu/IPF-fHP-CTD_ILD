rm(list=ls())
Sys.setenv(LANGUAGE = "en") #显示英文报错信息

getwd()
setwd('C:/Users/LWH/Downloads/bioinfo/IPF&CHD&SSc')
# source("CIBERSORT.R")
# devtools::install_github("omnideconv/immunedeconv", upgrade = T, force = T)
library(quadprog)
library(immunedeconv)
library(dplyr)
library(biomaRt)
library(ggplot2)
library(ggpubr)
library(CIBERSORT)
source("../AD/violin.R")
source("../AD/AD_functions.R")
load("exp.Rdata")
differ <- function(exprSet,group_list){
  suppressMessages(library(limma))
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exprSet)
  design
  fit <- lmFit(exprSet, design)
  group_list
  cont.matrix=makeContrasts(contrasts=c('Case-Control'),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  tempOutput = topTable(fit2, coef='Case-Control', n=Inf)
  DEG = na.omit(tempOutput)
  head(DEG) 
  return(DEG)
}

my_immunedeconv_IPF <- deconvolute_xcell(gene_expression_matrix = exp_IPF,
                                        arrays = F) ##arrays:芯片=T; RNA-seq=F
my_immunedeconv_CHP <- deconvolute_xcell(gene_expression_matrix = exp_CHP,
                                        arrays = F) ##arrays:芯片=T; RNA-seq=F
my_immunedeconv_CTD <- deconvolute_xcell(gene_expression_matrix = exp_CTD,
                                        arrays = TRUE) ##arrays:芯片=T; RNA-seq=F
# save(my_immunedeconv_IPF, my_immunedeconv_CHP, my_immunedeconv_CTD,
#      file = "my_immunedeconv.Rdata")
load(file = "my_immunedeconv.Rdata")
deg_IPF = differ(my_immunedeconv_IPF, group_list_IPF) 
deg_CHP = differ(my_immunedeconv_CHP, group_list_CHP)
deg_CTD = differ(my_immunedeconv_CTD, group_list_CTD)

cell_IPF <- deg_IPF[deg_IPF$adj.P.Val < .001,] %>% rownames()
cell_CHP <- deg_CHP[deg_CHP$adj.P.Val < .001,] %>% rownames()
cell_CTD <- deg_CTD[deg_CTD$adj.P.Val < .01,] %>% rownames()

del_cell <- c("Neurons", "Keratinocytes", "Sebocytes", "Preadipocytes", "Skeletal muscle",
              "Adipocytes", "Mesangial cells", "Plasma cells", "Megakaryocytes", "Chondrocytes",
              "Platelets", "Smooth muscle", "HSC", "CLP", "MPP","Melanocytes", "Hepatocytes",
              "Fibroblasts", "Osteoblast", "GMP","MSC", "Erythrocytes", "Astrocytes",
              "CMP","Pericytes","StromaScore","Myocytes") #CLP 共同淋巴样前体细胞; 多能造血祖细胞(MPP)
cells <- intersect(cell_IPF, cell_CHP) %>% intersect(cell_CTD) %>% 
  union(intersect(cell_IPF, cell_CHP)) 

deg_IPF[deg_IPF$adj.P.Val < .001 & deg_IPF$logFC > 0,] %>% 
  .[!rownames(.) %in% del_cell,] %>% rownames()
deg_CHP[deg_CHP$adj.P.Val < .001 & deg_CHP$logFC > 0,] %>% 
  .[!rownames(.) %in% del_cell,] %>% rownames()
deg_CTD[deg_CTD$adj.P.Val < .05 & deg_CTD$logFC > 0,] %>% 
  .[!rownames(.) %in% del_cell,] %>% rownames()

cell1<-deg_IPF[deg_IPF$adj.P.Val < .001,] %>% 
  .[!rownames(.) %in% del_cell,] %>% rownames() %>% .[1:5]
cell2<-deg_CHP[deg_CHP$adj.P.Val < .001,] %>% 
  .[!rownames(.) %in% del_cell,] %>% rownames() %>% .[1:5]
cell3<-deg_CTD[deg_CTD$adj.P.Val < .01,] %>%
  .[!rownames(.) %in% del_cell,] %>% rownames()

cell4<-deg_IPF[deg_IPF$adj.P.Val < .001,] %>% 
  .[!rownames(.) %in% c(cells,del_cell),] %>% rownames()
cell4<-cell4[-c(2,4,7)]
cell5<-deg_CHP[deg_CHP$adj.P.Val < .001,] %>% 
  .[!rownames(.) %in% cells,] %>% rownames()
cell6<-deg_CTD[deg_CTD$adj.P.Val < .076,] %>% 
  .[!rownames(.) %in% c(cells,cell3),] %>% rownames()



get_P <- function(data,Data_summary,label,n=-0.4){
  P_violin <- ggplot(data=data,aes(x=Group,y=value,fill=Attribute)) +
    geom_split_violin(trim=F,color="white",width=1.2,scale="width") + #绘制分半的小提琴图 +
    stat_compare_means(aes(group=Attribute), label = "p.signif",vjust=n) +
    geom_point(data=Data_summary,aes(x=Group,y=value),pch=19,position=position_dodge(0.9),size=1.5)+ #绘制均值为点图
    geom_errorbar(data = Data_summary,aes(ymin = value-ci, ymax=value+ci), #误差条表示95%的置信区间
                  width=1, #误差条末端短横线的宽度
                  position=position_dodge(0.9), 
                  color="black",
                  alpha = 0.7,
                  size=0.5) +
    scale_fill_manual(values = c("#56B4E9", "#E69F00"),
                      breaks=c("Control","Case"),labels=c("Control",label))+
    scale_y_sqrt(breaks=scales::breaks_width(0.2,offset=-0.1))+
    theme_bw()+ #背景变为白色
    theme(axis.text.x=element_text(angle=45,hjust =1,colour="black",family="serif",size=11), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="serif",size=10,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="serif",size = 16,face="plain"), #设置y轴标题的字体属性
          panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
          legend.text=element_text(face="italic", family="serif", colour="black",  #设置图例的子标题的字体属性
                                   size=16),
          legend.title=element_text(face="italic", family="serif", colour="black", #设置图例的总标题的字体属性
                                    size=18),
          panel.grid.major = element_blank(),   #不显示网格线
          panel.grid.minor = element_blank())+  #不显示网格线
    ylab("Fraction")+xlab("")+ #ylim(-0.1,2)+  #设置x轴和y轴的标题
    guides(fill=guide_legend(title="Group"))
  
}


get_Py <- function(data,Data_summary,label,n=-0.4){
  P_violin <- ggplot(data=data,aes(x=Group,y=value,fill=Attribute)) +
    geom_split_violin(trim=F,color="white",width=1.2,scale="width") + #绘制分半的小提琴图 +
    stat_compare_means(aes(group=Attribute), label = "p.signif",vjust=n) +
    geom_point(data=Data_summary,aes(x=Group,y=value),pch=19,position=position_dodge(0.9),size=1.5)+ #绘制均值为点图
    geom_errorbar(data = Data_summary,aes(ymin = value-ci, ymax=value+ci), #误差条表示95%的置信区间
                  width=1, #误差条末端短横线的宽度
                  position=position_dodge(0.9), 
                  color="black",
                  alpha = 0.7,
                  size=0.5) +
    scale_fill_manual(values = c("#56B4E9", "#E69F00"),
                      breaks=c("Control","Case"),labels=c("Control",label))+
    scale_y_sqrt(breaks=scales::breaks_width(0.45,offset=-0.1))+
    theme_bw()+ #背景变为白色
    theme(axis.text.x=element_text(angle=45,hjust =1,colour="black",family="serif",size=11), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="serif",size=10,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="serif",size = 16,face="plain"), #设置y轴标题的字体属性
          panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
          legend.text=element_text(face="italic", family="serif", colour="black",  #设置图例的子标题的字体属性
                                   size=16),
          legend.title=element_text(face="italic", family="serif", colour="black", #设置图例的总标题的字体属性
                                    size=18),
          panel.grid.major = element_blank(),   #不显示网格线
          panel.grid.minor = element_blank())+  #不显示网格线
    ylab("Fraction")+xlab("")+ #ylim(-0.1,2)+  #设置x轴和y轴的标题
    guides(fill=guide_legend(title="Group"))
  
}

immune_box_IPF1 = immunebox(my_immunedeconv_IPF[cell1,], group_list_IPF)
immune_box_CHP1 = immunebox(my_immunedeconv_CHP[cell2,], group_list_CHP)
immune_box_CTD1 = immunebox(my_immunedeconv_CTD[cell3,], group_list_CTD)
Data_summary1 <- summarySE(immune_box_IPF1, measurevar="value", 
                           groupvars=c("Group","Attribute"))
Data_summary2 <- summarySE(immune_box_CHP1, measurevar="value", 
                           groupvars=c("Group","Attribute"))
Data_summary3 <- summarySE(immune_box_CTD1, measurevar="value", 
                           groupvars=c("Group","Attribute"))
P1<-get_P(immune_box_IPF1,Data_summary1,"IPF")
P2<-get_P(immune_box_CHP1,Data_summary2,"fHP")
P3<-get_P(immune_box_CTD1,Data_summary3,"CTD_ILD",n=-2)

table(immune_box_CHP1$Group)
immune_box_CHP1[immune_box_CHP1$Attribute=="Control" &
                  immune_box_CHP1$Group=="Memory.B.cells",]$value %>% .[order(.)]



immune_box_IPF2 = immunebox(my_immunedeconv_IPF[cell4,], group_list_IPF)
library(magrittr)
immune_box_IPF2[immune_box_IPF2$Group=="MEP",]$value = 
  (immune_box_IPF2[immune_box_IPF2$Group=="MEP",]$value)*2
immune_box_CHP2 = immunebox(my_immunedeconv_CHP[cell5,], group_list_CHP)
immune_box_CTD2 = immunebox(my_immunedeconv_CTD[cell6,], group_list_CTD)
Data_summary4 <- summarySE(immune_box_IPF2, measurevar="value", 
                           groupvars=c("Group","Attribute"))
Data_summary5 <- summarySE(immune_box_CHP2, measurevar="value", 
                           groupvars=c("Group","Attribute"))
Data_summary6 <- summarySE(immune_box_CTD2, measurevar="value", 
                           groupvars=c("Group","Attribute"))

P4<-get_Py(immune_box_IPF2,Data_summary4,"IPF",n= -1)
P5<-get_P(immune_box_CHP2,Data_summary5,"CHP",n=-1.5)
P6<-get_P(immune_box_CTD2,Data_summary6,"CTD_ILD",n=-3.3)

library(patchwork)
P1/P2/P3
P4/P5/P6
