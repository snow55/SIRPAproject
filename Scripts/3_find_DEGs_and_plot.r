library(Seurat)
library(ggplot2)
load('./immuneCell_pure/recluster/Sdata_recluster.Rda')

## DEGs of different cell types
Idents(Sdata)="Tcelltype"
DFgenes <- FindAllMarkers(Sdata,only.pos = TRUE)#
DFgenes <- DFgenes[DFgenes$p_val_adj<0.05,]
write.table(DFgenes,file=paste0("./immuneCell_pure/recluster/","Tcelltype_DF.csv"),sep = ",",col.names = T,row.names = T)

## DEGs of tumor VS other
subNeu = subset(Sdata,cells=rownames(Sdata@meta.data)[which(Sdata$Tcelltype=="Macrophage")])
Idents(subMac)="group"
DFgenes=FindMarkers(subMac,ident.1 = "tumor",ident.2 = "non-tumor",logfc.threshold = 0.01,min.pct = 0.01)
write.table(DFgenes,file=paste0(sampleoutpath_subMye_pure,"DFgenes_tumor_VS_other_subMac.csv"),sep = ",",col.names = T,row.names = T)

## volcano plot
DFgene = read.table(file=paste0(sampleoutpath_subMye_pure,"DFgenes_tumor_VS_other_subMac.csv"),sep = ",")
sampleoutpath_volcano = paste0(sampleoutpath_subMye_pure,"volcanoPlot/")
if(!dir.exists(sampleoutpath_volcano)){dir.create(sampleoutpath_volcano,recursive = T)}
genes.to.label = c("SIRPA","PILRA","LILRB2","CLEC4A","LILRB3","LILRA5","SIGLEC9")#selected sig ITIMreporter genes

volcanoFUN = function(dataset=NULL,title=NULL,sampleoutpath=NULL,sample=NULL,cut_off_logFC=NULL,
                      labelUp=NULL,labelDown=NULL,w=4,h=4.6){
  # 设置pvalue和logFC的阈值
  cut_off_pvalue =0.05# 0.001
  cut_off_logFC = cut_off_logFC #0.5
  # 根据阈值分别为上调基因设置‘up’，下调基因设置‘Down’，无差异设置‘Stable’，保存到change列
  # 这里的change列用来设置火山图点的颜色
  # dataset = DFIN
  dataset = dataset[dataset$p_val_adj!=1,]
  dataset$gene = rownames(dataset)
  dataset$change = ifelse(dataset$p_val_adj < cut_off_pvalue & abs(dataset$avg_logFC) >= 0.1, 
                          ifelse(dataset$avg_logFC> 0.1 ,labelUp,labelDown),
                          '')
  ## 将-log10(p-value)>100的值均写为100 
  # dataset[dataset$p_val_adj<1E-20,]$p_val_adj = 1E-20
  # 绘制火山图
  pdf(paste0(sampleoutpath,sample,"Volcanoplot_",title,".pdf"),width = w,height = h)
  p = ggplot(
    #设置数据
    dataset,aes(x = avg_logFC,
                y = -log10(p_val_adj),
                colour=change)) +
    geom_point(alpha=1, size=0.5) +
    scale_color_manual(values=c("grey","#ff4757","#546de5"),
                       breaks=c(labelUp,labelDown),
                       labels=c(labelUp,#expression("IFN"^"hi"*" T special")
                                labelDown)
    )+ 
    
    # 辅助线
    # geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6) +
    geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="grey",lwd=0.6) +
    
    # 坐标轴
    labs(x=expression("Log"["2"]*"(Fold change)"),
         y=expression("-Log"["10"]*"(adjust P vaule)"),
         title = title)+ 
    theme_bw()+
    
    # 图例
    theme(plot.title = element_text(hjust = 0.5), 
          legend.text.align = 0,
          legend.title = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(size = 13),
          axis.title = element_text(size=15),
          legend.position = "top",
          legend.text = element_text(size = 13)) +
    guides(colour = guide_legend(override.aes = list(size=3),reverse = T)) #可使得legend中的圆点变大
  
  library(ggrepel)
  # 将需要标记的基因放置在label列
  genes.to.label = c("SIRPA","PILRA","LILRB2","CLEC4A","LILRB3","LILRA5","SIGLEC9")#selected sig ITIMreporter genes
  dataset$label = ifelse((dataset$p_val_adj < cut_off_pvalue & abs(dataset$avg_logFC) >= 1.5 ) | (dataset$gene %in% genes.to.label), as.character(dataset$gene),"")
  dataset$label[dataset$label=="IGHA2"]=""
  p=p+geom_text_repel(data = dataset, aes(x = avg_logFC, 
                                          y = -log10(p_val_adj), 
                                          label = label),
                      segment.color = "black",
                      show.legend = FALSE)
  print(p)
  dev.off()
}
volcanoFUN(dataset = DFgene[DFgene$p_val_adj!=1,],
           title = "tumor_VS_other_2",
           sampleoutpath = sampleoutpath_volcano,
           cut_off_logFC=0.1,
           sample = "subMac",
           labelUp="tumor high",
           labelDown = "tumor low",w=5,h=5.6)