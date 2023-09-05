library(Seurat)
library(ggplot2)
library(dplyr)
library(plyr)
load('./immuneCell_pure/recluster/Sdata_recluster.Rda')

constitutePath <- paste0(sampleoutpath_recluster,"constitutePath/")#
if(!dir.exists(constitutePath)){dir.create(constitutePath,recursive = T)}

CellinEachcluster <- table(Sdata$Tcelltype,Sdata$sample)## total cell type
TCellinEachcluster <- as.data.frame(CellinEachcluster)
colnames(TCellinEachcluster) <- c("Cluster","Sample","Freq")
TCellinEachcluster <- TCellinEachcluster[order(TCellinEachcluster$Cluster,TCellinEachcluster$Cluster),]

#### Normalize the number of cell by each sample
normalFUN <- function(df=NULL){
  TC <- as.data.frame(tapply(df$Freq, df$Sample, sum))
  TC <- cbind(TC,rownames(TC))
  colnames(TC) <- c("Totalcell","Sample")
  for(sample in TC$Sample){
    df$Totalcell[df$Sample==sample] <- TC$Totalcell[TC$Sample==sample]
  }
  df$Freq_normalized <- (df$Freq/df$Totalcell)*10000
  return(df)
}
TCellinEachcluster <- normalFUN(df=TCellinEachcluster)

##########画百分比的堆砌图
library(plyr)
#immuneCell_pure 
#Tcelltype-samples
sample="allImmuneCell"
tissueCols = c("#0a8cf7","#febe08","#ed2b2b")
col9 = c( "#e72287","#377eb8",  "#4daf4a",  "#984ea3",  "#ff7f00",  "#ffff33",  "#a65628",  "#f781bf",  "#999999")
ce2 <- ddply(TCellinEachcluster,"Cluster",transform, percent_Freq = Freq_normalized / sum(Freq_normalized) * 100)
pdf(paste0(constitutePath,"Tcelltype-sample_bottomlegend_","Barplot_Percent_",sample,".pdf"),
                                                            height=2+ length(levels(Sdata$Tcelltype))/4, width = 4)
p <- ggplot(ce2,aes(x=Cluster,y=percent_Freq,fill=Sample))+geom_bar(stat="identity") +
  # scale_fill_brewer(palette = "Paired",
  scale_fill_manual(values = col9,
                    labels = c(paste0("CRC1",c("N","P","T")),
                                                paste0("CRC2",c("N","P","T")), paste0("CRC3",c("N","P","T")))) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size=13,color = "black"),
        axis.text.y = element_text(size=15,color = "black"),
        legend.text = element_text(size = 11,color = "black"),
        legend.title = element_text(size = 11,color = "black"),
        legend.position = "bottom",
        panel.border = element_rect(color = "black",size=0.8),
        legend.margin=margin(t = -0.3, unit='cm'),
        axis.title = element_text(size=13,color = "black")) +
  # legend精准调整
  guides(
    fill = guide_legend(
      title.position = "top",
      ncol = 3, 
      bycol = TRUE,
      reverse = F)
  )+
  labs(x="",y="Percentage of frequence(%)")+
  coord_flip()+
  scale_x_discrete(limits= rev(levels(Sdata$Tcelltype))) ##确定X轴的顺序
print(p)
dev.off()
#location
##按每个样本单独Normalization，再计算不同tissue中的各个细胞类型的百分比
TCellinEachcluster$location = paste0(substr(TCellinEachcluster$Sample,1,3),
                                     substr(TCellinEachcluster$Sample,5,5))
ce2 <- ddply(TCellinEachcluster,"Cluster",transform, percent_Freq = Freq_normalized / sum(Freq_normalized) * 100)
ce3=data.frame(
  Cluster = NULL,
  Sample = NULL,
  percent_Freq = NULL
)
for(i in unique(ce2$Cluster)){
  for(j in unique(ce2$location)){
    x = sum(ce2$percent_Freq[ce2$Cluster==i & ce2$location==j])
    df = data.frame(Cluster = i,Sample =j,percent_Freq =x)
    ce3 = rbind(ce3,df)
  }
}
pdf(paste0(constitutePath,"Tcelltype-location_bottomlegend_20210914_","Barplot_Percent_",sample,".pdf"),height=2+ length(levels(Sdata$Tcelltype))/4, width = 3)
p <- ggplot(ce3,aes(x=Cluster,y=percent_Freq,fill=Sample))+geom_bar(stat="identity") +
  scale_fill_manual(values = tissueCols)+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size=13,color = "black"),
        axis.text.y = element_text(size=15,color = "black"),
        legend.text = element_text(size = 11,color = "black"),
        legend.title = element_text(size = 11,color = "black"),
        legend.position = "bottom",
        panel.border = element_rect(color = "black",size=0.8),
        legend.margin=margin(t = -0.3, unit='cm'),
        axis.ticks.y = element_blank(),
        legend.direction ="vertical",#legend纵向排列
        axis.title = element_text(size=14,color = "black")) +
  # legend精准调整
  guides(
    fill = guide_legend(
      title.position = "top",
      nrow = 1, 
      byrow = F,
      reverse = F)
  )+
  labs(x="Celltype",y="")+#Percentage of frequence(%)
  coord_flip()+ guides(fill = guide_legend(title = 'Tissue'))+
  scale_x_discrete(limits= rev(levels(Sdata$Tcelltype))) 
print(p)
dev.off()

#### compare the proportion between different tissues by scDC ####
library(ggplot2)
library(dplyr)
library(tidyverse)
library(scDC)
library(broom.mixed)
library(lme4)
setwd("./scDC/")

## human immune cell
Sys.time()
## 1.calculate cell proportion and confidence intervals
load("./immuneCell_pure/recluster/immuneCell_pure.Rda")
cellTypes = Sdata$Tcelltype
cellTypes=factor(cellTypes,levels = c("Macrophage","Neutrophil","DC","Eosinophil" ,"Mast cell", "T cell","B cell","Plasma cell","ILC","NK"))
subject = paste0(Sdata$location,"_",Sdata$sample)
res_BCa = scDC_noClustering(cellTypes, subject,
                         calCI = TRUE, calCI_method = "BCa",nboot=1000,ncores = 10)
save(res_BCa,file = paste0("scDC_result_BCa_human.Rda"))

## 2.fit GLM/GLMM
source("~/data/tools/githubLib/scDC/R/fitGLM_my.R")#revise the cellTypes in original code to assign the reference cell type
load(paste0("scDC_result_BCa_human.Rda"))
condition=c(rep("CRCN",10*3),rep("CRCP",10*3),rep("CRCT",10*3))
# barplotCI(res_BCa, condition)

res_BCa_Mast=fitGLM_human_my(res_BCa, condition ) #regard Mast cell as reference, need to change fitGLM.R to fitGLM_my.R
save(res_BCa_Mast,file = paste0("scDC_fitGLM_result_BCa_human_Mast.Rda"))
Sys.time()

## 3. plot
sample="human"
condition=c(rep("CRCN",10*3),rep("CRCP",10*3),rep("CRCT",10*3))

## plot barplot or density plot or cell type proportion and confidence interval
load(paste0("./scDC/scDC_result_BCa_",sample,".Rda"))
barplotCI(res_BCa,condition)
densityCI(res_BCa,condition)

## write out glm result
load(paste0("./scDC/scDC_fitGLM_result_BCa_",sample,"_NaiveT.Rda"))
glmresult=summary(res_BCa_Mast$pool_res_fixed)
write.csv(glmresult,file = paste0("./scDC/glmresult_",sample,".csv"),row.names = F)

## plot boxplot of cell type proportion
# data preparation
load(paste0("./scDC/scDC_result_BCa_",sample,".Rda"))
res=res_BCa
df1=as.data.frame(res$thetastar)
df1$celltype=res$results$cellTypes
df1$subject=res$results$subject

library(reshape2)
df=melt(df1,value.name = c("cellProportion"))
df$tissue=substr(df$subject,1,4) #for human data

colnames(df)=c("Cluster","subject","variable","Freq_normalized","location")
head(df)
df$Cluster=factor(df$Cluster,levels = c("Macrophage","Neutrophil","DC","Eosinophil" ,"Mast cell",
                                          "T cell","B cell","Plasma cell","ILC","NK")) #human

tissueCols = c("#0a8cf7","#febe08","#ed2b2b") #human 

# Make box plots. 
plot_list = list() 
for(i in 1:length(levels(df$Cluster))){
  # i=1
  ct=levels(df$Cluster)[i]
  tmp=df[df$Cluster==ct,]
  ymax=max(tmp$Freq_normalized)*100
  compaired <- list(c("CRCP","CRCN"),c("CRCT","CRCP"),c("CRCT", "CRCN"))
  p=ggplot(tmp,aes(x=location,y=Freq_normalized*100,fill=location)) + 
    geom_boxplot(width=0.7,alpha=1,color=tissueCols,notch=F,notchwidth=0.7) + 
    stat_boxplot(geom="errorbar",width=0.5,linetype=1)+ #添加error bar
    geom_boxplot(linetype=1)+
    stat_summary(fun="median",geom="point",shape=23,size=1.5,fill="white")+
    theme_classic() +
    theme(plot.title=element_text(size = 13,hjust = 0.5,colour = "black"),
          axis.text.x=element_text(size=13,angle=45,hjust = 1,colour = "black"),
          axis.text.y=element_text(size=12,colour = "black"),
          axis.title=element_text(size = 13,colour = "black"),
          axis.line.x.bottom = element_line(colour = "black",size=0.6),
          axis.line.y.left = element_line(colour = "black",size=0.6),
          panel.grid = element_blank(),
          legend.position = "none")+
    ylim(c(0,NA))+
    labs(title = ct,x='', y= paste0('Cell type proportion (%)')) +
    scale_fill_manual(values = tissueCols)+
    scale_color_manual(values = tissueCols)+
    geom_signif(data=tmp,comparisons = compaired,
                tip_length = 0.01,#连线的长度
                step_increase = 0.1,map_signif_level = T,vjust = 0.4,
                # map_signif_level=function(s)sprintf("p = %.2g", s),
                test = "t.test") #仅用于画出标注线，具体检验结果需根据glm result来修改
  plot_list[[i]] = p 
}

library(gridExtra)
pdf(paste0("./scDC/boxplot_res_",sample,".pdf"),w=7.2,h=6,useDingbats = F) 
grid.arrange(grobs=plot_list,top = NULL,nrow = 2) 
dev.off()


