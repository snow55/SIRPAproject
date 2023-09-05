library(ggplot2)
library(Seurat)
library(dplyr)
library(qugsea)

######## calculate immune inhibitory geneset core ########
Genes = readLines(paste0("/data/jinwf/wangxf/projects/CRCprojects/ImmuneInhibitory_selected.txt"))
Genes = intersect(Genes,rownames(Sdata[["RNA"]]@data))
## 1.construct dataframe to plot
cells = rownames(Sdata@meta.data)[Sdata$location == "CRCT"]
subData = as.matrix(t(Sdata[["RNA"]]@data[,cells]))
avg_matrix = data.frame(
  row.names = rownames(subData),
  immuneInhibitory = apply(subData[,Genes],1,mean)
)
avg_matrix$Tcelltype = Sdata$Tcelltype[rownames(avg_matrix)]
write.table(avg_matrix,paste0(sampleoutpath_recluster,"allImmunCell_InhibitScoreMatrix_onlyCRCT.csv"))

sampleoutpath_recluster = "./immuneCell_pure/recluster/"
avg_matrix=read.table(paste0(sampleoutpath_recluster,"allImmunCell_InhibitScoreMatrix_onlyCRCT.csv"))

avg_matrix$Tcelltype = factor(avg_matrix$Tcelltype,levels=c("Macrophage","Neutrophil","DC","T cell","ILC","NK",
                                                            "Basophil","Mast cell","B cell","Plasma cell"))
## 2.plot
pdf(paste0(sampleoutpath_recluster,"Immuneinhibitory_allImmuneCell_onlyCRCT.pdf"),w=5.5,h=4.5)
par(mar=c(8, 5, 1, 1) + 0.1)
boxplot(immuneInhibitory ~ Tcelltype, data = avg_matrix, 
    col=c("#FD6349","#00ab9f","#2fa3e9","#dbbca7","#31CC30","#078309","#be2583","#FF6FB5","#cb6318","#d5c0e0"),
        par(las="2"),
        cex.lab = 1.5,cex.axis=1.5,
        xlab="",
        xaxt="n", #不画横轴的刻度和标签
        ylab="Immune inhibitory score",
        outline = FALSE     ## avoid double-plotting outliers, if any
)
axis(side=1, at=1:10, labels=FALSE) #增加横轴刻度
text(levels(avg_matrix$Tcelltype), x=1:10,cex=1.3, y=-0.08,adj = 0.9, xpd=T, srt=45) #调整横轴文字的大小、方向和位置
dev.off()


####### compare the difference of one geneset between two groups by QuSAGE ########
#http://bioconductor.riken.jp/packages/3.0/bioc/vignettes/qusage/inst/doc/qusage.pdf
install.packages("http://bioconductor.riken.jp/packages/3.0/bioc/src/contrib/qusage_1.6.0.tar.gz", repos=NULL, type="source")
library(qusage)
# transfer p value for base R plotting
transPvalue=function(p0){
  # p0=1.2e-3; p0
  if (p0<2.2e-16) {
    p0=2.2e-16
    p0=formatC(p0, format="e", digits=2);#p0 #"1.2e-52"
    # split by e
    p0=strsplit(p0, "e")[[1]]; #p0 #[1] "1.2" "-52"
    label = bquote(italic(P)~"<"~.(p0[1])~"\u00d7"~10^.(p0[2]) )
  }else if(p0<0.01){
    p0=formatC(p0, format="e", digits=2);#p0 #"1.2e-52"
    # split by e
    p0=strsplit(p0, "e")[[1]]; #p0 #[1] "1.2" "-52"
    # plot
    label = bquote(italic(P)~"="~.(p0[1])~"\u00d7"~10^.(as.numeric(p0[2]) ))
  }else{
    p0=round(p0,3); #p0 #保留3位
    label = bquote(italic(P)~"=" ~ .(p0) )
  }
  return(label)
}

#### for single geneSet -- chemokine activity
chemokine = readLines("chemokine.txt")#genelist
## 1.calculation
#mMDSC, TAM, gMDSC
eset=as.data.frame(Sdata[["RNA"]]@data[,rownames(Sdata@meta.data)[which(Sdata$Scelltype %in% c("mMDSC", "TAM", "gMDSC"))]])
labels=as.character(paste(Sdata$Scelltype[colnames(eset)],Sdata$group[colnames(eset)],sep="."))
sf_mMDSC=qusage(eset,labels,"mMDSC.KO-mMDSC.WT",chemokine)
sf_TAM=qusage(eset,labels,"TAM.KO-TAM.WT",chemokine)
sf_gMDSC=qusage(eset,labels,"gMDSC.KO-gMDSC.WT",chemokine)

# p value adjustment
p.vals_TAM=pdf.pVal(sf_TAM)
q.vals_TAM = p.adjust(p.vals_TAM, method="fdr")
p.vals_gMDSC=pdf.pVal(sf_gMDSC)
q.vals_gMDSC = p.adjust(p.vals_gMDSC, method="fdr")
p.vals_mMDSC=pdf.pVal(sf_mMDSC)
q.vals_mMDSC = p.adjust(p.vals_mMDSC, method="fdr")

## 2.plot
pdf(paste0(sampleoutpath,"genesetScore_","mMDSC-TAM-gMDSC","_QuSAGE.pdf"),w=5.5,h=4.5,useDingbats = F)
par(oma=c(1,1,1,1),mar=c(4,4.5,3,0))
plotDensityCurves(sf_TAM,col=ScelltypeCols[2],lwd=2, cex.axis=1.2,cex.lab=1.6,cex.main=1.6, 
                  main=expression("Chemokine activity: "*italic('Sirpa')^italic('-/-')*' vs WT') ,#
                  xlab=expression(Log[2](Activity)),ylab="Density",xlim = c(-0.005,0.07),ylim=c(0,150))
plotDensityCurves(sf_gMDSC,col=ScelltypeCols[5],lwd=2, xlab=expression(log[2](Activity)),add = T)
plotDensityCurves(sf_mMDSC,col=ScelltypeCols[1],lwd=2, xlab=expression(log[2](Activity)),add = T)
box(lwd=2)
legend("topright", legend=c("TAM","gMDSC","mMDSC"),lty=1,lwd=2,col=ScelltypeCols[c(2,5,1)],cex=1.2)
text(x = c(0.061), y = c(80), labels =transPvalue(q.vals_TAM), col = ScelltypeCols[2],cex = c(0.8))
text(x = c(0.01), y = c(140), labels = transPvalue(q.vals_gMDSC),col = ScelltypeCols[5],cex = c(0.8))
text(x = c(0.005), y = c(50), labels =transPvalue(q.vals_mMDSC), col = ScelltypeCols[1],cex = c(0.8))
dev.off()
