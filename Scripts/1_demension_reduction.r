library(Seurat)
library(ggplot2)
library(dplyr)
library(plyr)

######### 1.basic analysis pipeline #####
reclusterFun = function(Sdata,sample,sampleoutpath,nfeatures=3000,resolution=0.8){
  Sdata<- CreateSeuratObject(counts = Sdata[["RNA"]]@counts, min.cells = 0, min.features = 0)
  dim(Sdata[["RNA"]]@counts)
  Sdata[["percent.mt"]] <- PercentageFeatureSet(Sdata, pattern = "^MT-")
  head(Sdata@meta.data, 5)
  
  Sdata <- NormalizeData(Sdata, normalization.method = "LogNormalize", scale.factor = 10000)
  Sdata <- FindVariableFeatures(Sdata, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(Sdata), 10)
  pdf(paste0(sampleoutpath,"VarFeatureplot_filtered_",sample,".pdf"),width = 12,height = 6)
  plot1 <- VariableFeaturePlot(Sdata)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  print(CombinePlots(plots = list(plot1, plot2)))
  dev.off()
  all.genes <- rownames(Sdata)
  Sdata <- ScaleData(Sdata, features = all.genes)
  RPgenes <- rownames(Sdata[["RNA"]]@data)[grep("^RP",rownames(Sdata[["RNA"]]@data))]
  MTgenes <- rownames(Sdata[["RNA"]]@data)[grep("^MT-",rownames(Sdata[["RNA"]]@data))]
  Sdata <- RunPCA(Sdata, features = VariableFeatures(object = Sdata))
  
  PCAplotFun = function(sample){
    pdf(paste0(sampleoutpath,"PCAplot_filtered_",sample,"_0.pdf"),width = 5.5,height = 5)
    p1 <- VizDimLoadings(Sdata, dims = 1:2, reduction = "pca")
    p2 <- DimPlot(Sdata, reduction = "pca", dims = c(1, 2))
    p3 <- DimPlot(Sdata, reduction = "pca", dims = c(2, 3))
    p4 <- DimHeatmap(Sdata, dims = 1, cells = 500, balanced = TRUE)
    p5 <- ElbowPlot(Sdata,ndims = 50)
    print(list(p1, p2,p3,p4,p5))
    dev.off()
  }
  PCAplotFun(sample)
  
  Sdata <- FindNeighbors(Sdata, dims = 1:50)
  Sdata <- FindClusters(Sdata, resolution = 0.8)
  
  Sdata <- RunTSNE(Sdata, dims = 1:50,check_duplicates = FALSE)
  pdf(paste0(sampleoutpath,"tSNEplot_filtered_",sample,".pdf"),width = 5,height = 5)
  p <- DimPlot(Sdata, reduction = "tsne",label = TRUE,repel = TRUE,pt.size = 0.001)+ NoLegend()
  print(p)
  dev.off()

  Sdata <- RunUMAP(Sdata, reduction = "pca", dims = 1:50, min.dist = 0.1)
  pdf(paste0(sampleoutpath,"UMAPplot_filtered_batch_",sample,".pdf"),width = 5,height = 5)
  p <- DimPlot(Sdata, reduction = "umap",label = TRUE,repel = TRUE,pt.size = 0.001) + NoLegend()
  print(p)
  dev.off()
  
  save(Sdata,file = paste0(sampleoutpath,sample,".Rda"))
}
reclusterFun(Sdata,sample="Sdata_recluster",sampleoutpath="./immuneCell_pure/recluster/",nfeatures=3000,resolution=0.8)

######### 2.further use FIt-SNE #########
load('./immuneCell_pure/recluster/Sdata_recluster.Rda')
PC50 <- Sdata@reductions$pca@cell.embeddings
PCAinit = PC50[,1:2]/sd(PC50[,1]) * 0.0001
Sdata <- RunTSNE(Sdata,tsne.method = "FIt-SNE",  dims = 1:50,
                 initialization=PCAinit, learning_rate=1000,
                 fast_tsne_path = "/home/wangxf/data/tools/FIt-SNE-master/bin/fast_tsne",
                 reduction.key = "FIt_SNE")

######### 3.plot #########
mycol_subMye = c('#fdc086','#f2991a',  '#cc051a', '#cf6eb8', '#9c21b8','#d94c1a','#e6cc4c', '#b2d199','#b7a39f', '#72daf2', '#3387b5',
                 '#beaed4', '#76c6ba', '#7fc97f','#1a734c', '#dbb972',"#3abca5","#b49674","#e2a226","#e26a46","#e58f70","#d7ce96","#e57371","#e76da2","#c5aac3","#be71af")
tissueCols = c("#0a8cf7","#febe08","#ed2b2b")
TcelltypeCols=c("#e11f27", "#359e4a","#1f78b3","#6b3f97", "#f79d9c",  "#dbbca7",  "#a6cce0",  "#f47f21","#b2d589","#d5c0e0")
sampleoutpath_recluster="./immuneCell_pure/recluster/"
splitPath = paste0("./immuneCell_pure/recluster/splitPath/")
#allImmune cell
pdf(paste0(sampleoutpath_recluster,"splitPath/","Batch_Tcelltype_","allImmu_pure","_","tsne",".pdf"),width = 5.8,height = 4.3)#
p <- DimPlot(Sdata, group.by = "Tcelltype",label=T,label.size=4.2,repel = T,pt.size=0.3,reduction = "tsne",cols=TcelltypeCols)  + #
  theme_bw() + theme(panel.grid = element_blank() ,axis.text = element_text(size = 13,color="black"),
                     panel.border = element_rect(color="black",size=0.8),
                     axis.title = element_text(size = 15,color="black"),
                     legend.text = element_text(size = 13,color="black"),
                     legend.key.height=unit(0.7,"cm"),
                     legend.key.width=unit(2,'mm'),
                     legend.margin=margin(b = -0.3, unit='cm'),
                     legend.title = element_blank())+ 
  labs(x="tSNE 1",y="tSNE 2")
print(p)
dev.off()
##location
pdf(paste0(sampleoutpath_recluster,"splitPath/","Batch_location_","allImmu_pure","_","tsne",".pdf"),width = 5.3,height = 4.3)#
p <- DimPlot(Sdata, group.by = "location",label=F,label.size=4.2,repel = T,pt.size=0.3,reduction = "tsne",cols=tissueCols)  + #
  theme_bw() + theme(panel.grid = element_blank() ,axis.text = element_text(size = 13,color="black"),
                     panel.border = element_rect(color="black",size=0.8),
                     axis.title = element_text(size = 15,color="black"),
                     legend.text = element_text(size = 13,color="black"),
                     legend.key.height=unit(0.7,"cm"),
                     legend.key.width=unit(2,'mm'),
                     legend.margin=margin(b = -0.3, unit='cm'),
                     legend.title = element_blank())+ 
  labs(x="tSNE 1",y="tSNE 2")
print(p)
dev.off()

#subMye
pdf(paste0(sampleoutpath_subMye_pure,"splitPath/","Batch_location_","subMye_pure","_","tsne","_20230614.pdf"),width = 8.15,height = 4.3)#
p <- DimPlot(Sdata, group.by = "celltype",label=T,label.size=4.2,repel = T,pt.size=0.3,reduction = "tsne")  + #,cols=mycol_subMye
  theme_bw() + theme(panel.grid = element_blank() ,axis.text = element_text(size = 13,color="black"),
                     panel.border = element_rect(color="black",size=0.8),
                     plot.title = element_blank(),
                     legend.text = element_text(size = 13,color="black"),
                     legend.key.height=unit(0.7,"cm"),
                     legend.key.width=unit(2,'mm'),
                     legend.margin=margin(b = -0.3, unit='cm'),
                     legend.title = element_blank())+ 
  scale_color_manual(values=mycol_subMye,labels=as.character(levels(Sdata$newclusterLabel)))+
  labs(x="tSNE 1",y="tSNE 2",color="Tissue")
print(p)
dev.off()
