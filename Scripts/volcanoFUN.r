## basic volcanoplot
volcanoFUN = function(dataset=NULL,title=NULL,
                      sampleoutpath=NULL,sample=NULL,
                      labelUp=NULL,labelDown=NULL,
                      cut_off_pvalue = 0.001,
                      cut_off_logFC = 1,
                      label_logFC=0.8,
                      w=5.3,h=4){
  # 设置pvalue和logFC的阈值
  # cut_off_pvalue = 0.0000001#0.05
  # cut_off_logFC = 0.5 #0.5
  # 根据阈值分别为上调基因设置‘up’，下调基因设置‘Down’，无差异设置‘Stable’，保存到change列
  # 这里的change列用来设置火山图点的颜色
  dataset$gene = rownames(dataset)
  dataset$change = ifelse(dataset$p_val_adj < cut_off_pvalue & abs(dataset$avg_logFC) >= cut_off_logFC, 
                          ifelse(dataset$avg_logFC> cut_off_logFC ,labelUp,labelDown),
                          '')
  dataset$change = factor(dataset$change,levels=c("",labelDown,labelUp))
  ## 将-log10(p-value)>100的值均写为100 
  dataset[dataset$p_val_adj<1E-100,]$p_val_adj = 1E-100
  # 绘制火山图
  pdf(paste0(sampleoutpath,sample,"Volcanoplot_",title,".pdf"),width = w,height = h)
  p = ggplot(
    #设置数据
    dataset,aes(x =avg_logFC,  #dataset$meanWT,
                y = -log10(p_val_adj), #dataset$meanKO,
                colour=change)) +
    geom_point(alpha=1, size=0.5) +
    scale_color_manual(values=c("grey","#546de5","#ff4757"),
                       breaks=c(labelUp,labelDown),
                       labels=c(labelUp,
                                labelDown)
    )+
    
    # 辅助线
    geom_vline(xintercept=c(-0.1,0.1),lty=4,col="black",lwd=0.6) +
    geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.6) +
    
    # 坐标轴
    labs(x=expression("Log"["2"]*"(Fold change)"),#"mean expression in WT",
         y=expression("-Log"["10"]*"(P vaule)"),#"mean expression in KO",
         title = title)+ 
    theme_bw()+
    
    # 图例
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position="right", 
          legend.text.align = 0,
          legend.title = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(size = 13),
          axis.title = element_text(size=15),
          legend.text = element_text(size = 13)) +
    guides(colour = guide_legend(override.aes = list(size=3)))#+
  # scale_x_continuous(breaks = c(-1,-0.5,0,0.5,1))#+ylim(0,100)
  
  library(ggrepel)
  # 将需要标记的基因放置在label列
  # 这里设置logFC值大于5的差异基因来标记
  # ！！！需要注意的是标记的基因不能太多，Rstudio容易卡死
  dataset$label = ifelse(dataset$p_val_adj < cut_off_pvalue & abs(dataset$avg_logFC) >= label_logFC, as.character(dataset$gene),"")
  p=p+geom_text_repel(data = dataset, aes(x = dataset$avg_logFC, #dataset$meanWT,
                                          y = -log10(dataset$p_val_adj), #dataset$meanKO,
                                          label = dataset$label),
                                          size = 3,
                                          box.padding = unit(0.1, "lines"),
                                          # point.padding = unit(0.8, "lines"), 
                                          segment.color = "black", 
                                          show.legend = FALSE)
  print(p)
  dev.off()
}

eg.
volcanoFUN(dataset = DFgenes_TAM[!DFgenes_TAM$p_val_adj==1,],
           title = "KO vs WT",
           sampleoutpath = sampleoutpath_DFgene,
           cut_off_logFC = 0.1,
           label_logFC=0.5,
           sample = "TAM",
           labelUp="KO high",
           labelDown = "WT high")
           
## beautiful volcanoplot
volcanoFUN_beautiful=function(dataset=NULL,title=NULL,
                              sampleoutpath=NULL,sample=NULL,
                              labelUp=NULL,labelDown=NULL,
                              cut_off_pvalue = 0.05,
                              cut_off_logFC = 0.1,
                              label_logFC=0.8,
                              label_pvalue=NULL,
                              genes.to.label=NULL,
                              p_val_min=NULL,
                              w=5.3,h=4){
  library(ggrepel)
  library(ggrastr)
  
  dataset = dataset[dataset$p_val_adj!=1,]
  dataset$gene = rownames(dataset)
  
  dataset$change = ifelse(dataset$p_val_adj < cut_off_pvalue & abs(dataset$avg_logFC) >= 0, 
                          ifelse(dataset$avg_logFC> 0 ,labelUp,labelDown),"") 
  # dataset$change = factor(dataset$change,levels = c("",labelUp,labelDown))
  if(!is.null(genes.to.label) & is.null(label_pvalue)){
    dataset$label = ifelse(dataset$gene %in% genes.to.label | 
                             (dataset$p_val_adj < cut_off_pvalue & abs(dataset$avg_logFC)>label_logFC), 
                           dataset$gene,"")
  } else if(is.null(genes.to.label) & !is.null(label_pvalue)){
    dataset$label = ifelse( dataset$p_val_adj < cut_off_pvalue & abs(dataset$avg_logFC)>label_logFC |
                              dataset$p_val_adj<label_pvalue, dataset$gene,"")
  }else if(!is.null(genes.to.label) & !is.null(label_pvalue)){
    dataset$label = ifelse(dataset$gene %in% genes.to.label | 
                             (dataset$p_val_adj < cut_off_pvalue & abs(dataset$avg_logFC)>label_logFC) |
                             dataset$p_val_adj<label_pvalue, 
                           dataset$gene,"")
  }else {
    dataset$label = ifelse( dataset$p_val_adj < cut_off_pvalue & abs(dataset$avg_logFC)>label_logFC, dataset$gene,"")
  }
  dataset = dataset[order(dataset$change,decreasing = T),]
  dataset$pointSize = ifelse(dataset$change==labelUp | dataset$change==labelDown,0.5,0.1)
  if(!is.null(p_val_min)){
    dataset[which(dataset$p_val_adj < p_val_min),]$p_val_adj=p_val_min
  }else{
    dataset$p_val_adj=dataset$p_val_adj
  }
  
  pdf(paste0(sampleoutpath,sample,"Volcanoplot_",title,".pdf"),width = w,height = h)
  p = ggplot(dataset,aes(x = avg_logFC, 
                         y = -log10(p_val_adj),
                         fill=change,color=change)) +
    # geom_point(size=2,alpha=1,shape=16) +
    geom_point_rast(size=2,alpha=1,shape=16, raster.dpi = getOption("ggrastr.default.dpi", 300),)+ #图片瘦身1：点图不是是矢量，文字是
    # scattermore::geom_scattermore(pixels = c(512, 512), pointsize = pt.size*50) #raster.dpi=c(512, 512) #图片瘦身2
    scale_color_manual(values=c("grey80","#a50f15","#4169e1"),
                       breaks=c(labelUp,labelDown),
                       labels=c(labelUp,
                                labelDown)
    )+ 
    scale_fill_manual(values=c("grey80","#a50f15","#4169e1")
    )+ 
    
    geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.6) +
    geom_vline(xintercept = c(0),lty=4,col="black",lwd=0.6) +
    labs(x=expression("Log"["2"]*"(Fold change)"),
         y=expression("-Log"["10"]*"(adjust P vaule)"),
         title = "")+ 
    theme_classic(base_size = 15) + 
    theme(plot.title = element_text(hjust = 0.5), 
          legend.text.align = 0,
          legend.title = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(size = 13,color = "black"),
          axis.title = element_text(size=15,color = "black"),
          axis.line.x.bottom = element_line(color = "black",size=0.8),
          axis.line.y.left = element_line(color = "black",size=0.8),
          legend.margin = margin(b=-0.3),
          legend.position = "top",
          legend.text = element_text(size = 13)) +
    guides(fill = "none",
           colour = guide_legend(override.aes = list(size=3),reverse = T))+ 
    geom_label_repel(aes(label=label), 
                     color="white", 
                     box.padding=unit(0.15, "lines"), 
                     point.padding=unit(0.5, "lines"), 
                     segment.size=0.2,segment.colour = "black")
  print(p) 
  dev.off()  
}

eg.
volcanoFUN_beautiful(dataset = DFgenes,
           title = "KO vs WT",
           sampleoutpath = sampleoutpath_DFgene,
           cut_off_logFC = 0.1,
           label_logFC=0.4,
           label_pvalue=NULL,
           p_val_min=1E-100,
           sample = ct,
           labelUp="KO high",
           labelDown = "WT high",
           w=5.3,h=6)
