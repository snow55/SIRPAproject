## multiple boxplot using ggplot2
library(ggplot2)
multiBoxplotFUN=function(Sdata=NULL,ct=NULL,Genes=NULL,
                         groupVar=NULL,cols.use=NULL,
                         compaireList=NULL,
                         jitterSize=0.8,
                         sampleoutpath=NULL,w=6,h=3.5,nrow=NULL,ncol=NULL,
                         plotName=NULL){
  # ct="Macrophage"
  cellUse=rownames(Sdata@meta.data)[Sdata$Tcelltype==ct]
  # Make plots. 
  plot_list = list() 
  for(i in 1:length(Genes)){
    gene=Genes[i]
    expression = Sdata[["RNA"]]@data[gene,]
    df = data.frame(
      exp = expression,
      group = Sdata@meta.data[,groupVar]
    )
    df=df[cellUse,]
    p=ggplot(df,aes(x=group,y=exp,fill=group)) + 
      # stat_summary(fun=mean, geom="point", size=1,color="white") +#shape=23, 
      # geom_dotplot(binwidth=0.05,binaxis='y',stackdir='center',dotsize = .5,fill="yellow") +
      # stat_boxplot(geom="errorbar",width=0.3,aes(color=group))+#
      geom_jitter(width=0.23,alpha=0.3,aes(color=group),size=jitterSize) +
      geom_violin(width=0.8,aes(color=group),trim = T) +
      geom_boxplot(width=0.15,outlier.shape = NA,color="black",fill="white",alpha=1) + #outlier.shape = NA,color="white" #,color=tissueCols,notch=TRUE,notchwidth=0.7
      theme_classic() +
      theme(plot.title=element_text(size = 13,hjust = 0.5,colour = "black"),
            axis.text.x=element_text(size=13,angle=45,hjust = 1,colour = "black"),
            axis.text.y=element_text(size=12,colour = "black"),
            axis.title=element_text(size = 13,colour = "black"),
            axis.line.x.bottom = element_line(colour = "black",size=0.6),
            axis.line.y.left = element_line(colour = "black",size=0.6),
            panel.grid = element_blank(),
            legend.position = "none")+
      labs(title = gene,x='', y= 'Expression level') +
      # facet_wrap(~df$celltype) +
      # scale_fill_manual(values = mycol_Mac)#+
      scale_fill_manual(values = cols.use)+
      scale_color_manual(values = cols.use)+
      geom_signif(data=df,comparisons = compaireList,
                  # y_position = c(0.31,0.27,0.18),#横线标记的位置c(0.3,0.28,0.26,0.12,0.1)
                  tip_length = 0.01,#连线的长度
                  step_increase = 0.1,map_signif_level = T,vjust = 0.4,
                  test = "wilcox.test")#"t.test"
    plot_list[[i]] = p 
  }
  library(gridExtra)
  pdf(paste0(sampleoutpath,ct,"_",plotName,"_box_exp_tissues.pdf"),
      w=w,h=h,useDingbats = F)
  grid.arrange(grobs=plot_list,
               top = ct,nrow = nrow, ncol=ncol) 
  dev.off()
}

eg.
multiBoxplotFUN(Sdata=subMye,ct="Macrophage",Genes=Genes,
                groupVar = "location",cols.use=tissueCols,
                compaireList=list(c("CRCP","CRCN"),c("CRCT","CRCP"),c("CRCT", "CRCN")),
                jitterSize=0.1,
                sampleoutpath=sampleoutpath_boxplot,
                plotName = "ITIM",
                w=9,h=4,ncol=6)
