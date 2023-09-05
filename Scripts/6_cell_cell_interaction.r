library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(Cairo)
options(stringsAsFactors = FALSE)
mycol_45 = c('#fdc086','#f2991a',  '#cc051a', '#cf6eb8', '#9c21b8','#d94c1a','#e6cc4c', '#b2d199','#b7a39f', '#72daf2', '#3387b5',#Macrophage
             '#beaed4', '#76c6ba', '#7fc97f','#1a734c','#dbb972', #DC
             "#3abca5","#b49674","#e2a226","#e26a46","#e58f70","#d7ce96","#e57371","#e76da2","#c5aac3",#Neutrophil
             "#be71af","#648bc9","#be2583",#Mast cell
             "#9b9a9b", "#42a16c", "#fdcf16", "#c26257", "#368b9d", "#9c425c", "#9452a0", "#bf6067",
             "#d4af29", "#e986b7", "#ea7425", "#1a79b5", "#2fa147", "#fdbf6d", "#b3d78b", #T cell
             "#f6999a", "#cab2d6", "#a6cfe4")#B cell
ScelltypeCols=c( "#b49674", "#eb0875", '#b2d199',"#e57371", "#5298e7",#'#3387b5', #20210917 renew
                 '#d94c1a','#beaed4',"#3abca5", '#f2991a',"#e58f70",
                 "#e76da2",'#7fc97f',"#d7ce96",'#72daf2',"#be71af")
load("./mergedImmuneCell_pure/mergedImmuneCell_recluster.Rda")

######## 1.basic cell-cell communication analysis
Idents(Sdata) = "Scelltype"
type = "Scelltype"
for(groupName in unique(Sdata$group)){
  # for(groupName in "ALL"){
  # groupName = "WT"
  # subData <- subset(Sdata,cells=rownames(Sdata@meta.data)[which(Sdata$group %in% groupName)])
  subData = Sdata
  data.input <- GetAssayData(subData, assay = "RNA", slot = "data") # normalized data matrix
  labels <- Idents(subData)
  meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
  
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
  levels(cellchat@idents) # show factor levels of the cell labels
  
  groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
  
  LRtypes = unique(CellChatDB.mouse$interaction$annotation)
  # [1] "Secreted Signaling" "ECM-Receptor"       "Cell-Cell Contact"
  for(i in 1:3){
    # i=2 
    cellchat@DB <-  subsetDB(CellChatDB.mouse, search = LRtypes[i])
    
    ## Preprocessing the expression data for cell-cell communication analysis
    cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
    future::plan("multiprocess", workers = 4) 
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat,PPI.mouse) # PPI.human

    # Part II: Inference of cell-cell communication network
    cellchat <- computeCommunProb(cellchat)  #time-consuming
    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat <- filterCommunication(cellchat, min.cells = 5)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    # Part III: Visualization and systems analysis of cell-cell communication network
    ## Create a directory to save figures
    subName = paste0(groupName)#"subMye_",
    data.dir <- paste0("/data/jinwf/wangxf/projects/CRCprojects/seuratv3_analysis_rawdata/MDSC/mergedMDSC/mergedImmuneCell_pure/cellchat_results/",subName,"_",type,"_allImmuCell_CellChat_results/",LRtypes[i],"/")
    if(!dir.exists(data.dir)){dir.create(data.dir,recursive = T)}
    setwd(data.dir)
    
    levels(cellchat@idents)
    
    # Hierarchy plot
    vertex.receiver = seq(1,ceiling(length(levels(cellchat@idents))/2))# a numeric vector

    # Compute the network centrality scores
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
    
    # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
    par(mfrow=c(1,2))
    pdf(file = paste0("Allpathways","_netAnalysis_signalingRole_scatter.pdf"),w=10,h=6)
    gg1 <- netAnalysis_signalingRole_scatter(cellchat,color.use = ScelltypeCols)
    #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
    # Signaling role analysis on the cell-cell communication networks of interest
    # gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
    #> Signaling role analysis on the cell-cell communication network from user's input
    print(gg1)
    # print(gg2)
    dev.off()
    
    # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
    par(mfrow=c(1,2))
    pdf(file = paste0("allsigPathways","_netAnalysis_signalingRole_heatmap.pdf"),w=10,h=6)
    ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",color.use = ScelltypeCols)
    ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",color.use = ScelltypeCols)
    # print(list(ht1 ,ht2))
    print(ht1)
    print(ht2)
    dev.off()
    
    par(mfrow=c(1,1))
    ## total ligand-receptor communication
    CairoPDF("netVisual_circle plot_top0.2_all.pdf",w=5,h=5)#
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize,top=0.1, weight.scale = T, label.edge= F, edge.label.cex = 0.8, 
                     vertex.label.cex = 1,color.use = ScelltypeCols)
    dev.off()
    CairoPDF("netVisual_circle plot_labeled_top0.2_all.pdf",w=5,h=5)
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize,top=0.1,weight.scale = T, label.edge= T, edge.label.cex = 0.8, 
                     vertex.label.cex = 1,color.use = ScelltypeCols)
    dev.off()
    
    CairoPDF("netVisual_circle plot_all_number.pdf",w=5,h=5)#
    netVisual_circle(cellchat@net$count,  weight.scale = T, label.edge= F, color.use = ScelltypeCols,vertex.weight = groupSize,
                     edge.label.cex = 0.8, vertex.label.cex = 1, title.name = "Number of interactions")
    dev.off()
    CairoPDF("netVisual_circle plot_all_strength.pdf",w=5,h=5)#
    netVisual_circle(cellchat@net$weight,  weight.scale = T, label.edge= F, color.use = ScelltypeCols,vertex.weight = groupSize,
                     edge.label.cex = 0.8, vertex.label.cex = 1, title.name = "Interaction weights/strength")
    dev.off()
    
    mat <- cellchat@net$weight
    par(mfrow = c(ceiling(length(unique(Idents(Sdata)))/4),4), xpd=TRUE)
    shellPath="./shellPath/"
    if(!dir.exists(shellPath)){dir.create(shellPath)}
    for (j in 1:nrow(mat)) {
      mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      mat2[j, ] <- mat[j, ]
      CairoPDF(paste0(shellPath,"netVisual_circle plot_", rownames(mat)[j],"_strength.pdf"),w=5,h=5)#
      netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[j])
      dev.off()
    }
    
    ## Part IV: Save the CellChat object
    cellchat <- rankNetPairwise(cellchat)
    saveRDS(cellchat, file = paste0("cellchat_",subName,"_",LRtypes[i],"_",type,"_allImmuCell.rds"))
    # cellchat <- readRDS(file=paste0("cellchat_",subName,".rds"))
  }
}

######## 2.visualize the cell-cell communication network
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(Cairo)
library(patchwork)
options(stringsAsFactors = FALSE)
cellchat <- readRDS(file=paste0("cellchat_",subName,".rds"))

mycol_subMye = c("#7fc67f","#e58f70",  "#ba7fb7","#5b86ab","#eb0875",'#ffe135',#Macrophage
                 '#e76da2',  "#b49674","#d531ae",'#72daf2', #DC
                 '#975fe3',"#a9a9a9","#ffc1cc","#b2aa0d",'#5298e7','#e32636','#f2991a',"#3abca5",#Neutrophil
                 '#dbb972',"#c5aac3","#e26a46")#Mast cell ##new colors
ScelltypeCols=c( "#b49674", "#eb0875", '#b2d199',"#e57371", "#5298e7",#'#3387b5', #20210917 renew
                 '#d94c1a','#beaed4',"#3abca5", '#f2991a',"#e58f70",
                 "#e76da2",'#7fc97f',"#d7ce96",'#72daf2',"#be71af")

####### compare two samples ########
##分不同种类LRpairs -- Scelltype
LRtypes = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact")
m=3
type="Scelltype"
# cellchatDir = "/data/jinwf/wangxf/projects/CRCprojects/seuratv3_analysis_rawdata/MDSC/mergedMDSC/mergedAllCell_pure/cellchat_results/"
cellchatDir = "/data/jinwf/wangxf/projects/CRCprojects/seuratv3_analysis_rawdata/MDSC/mergedMDSC/mergedImmuneCell_pure/cellchat_results/"
for(i in 1:2){
  if(i==1){
    subName = "KO"
    Scelltype_KO = readRDS(paste0(cellchatDir,subName,"_",type,"_allImmuCell_CellChat_results/",LRtypes[m],"/cellchat_",subName,"_",LRtypes[m],"_",type,"_allImmuCell.rds")) #
    Scelltype_KO <- rankNetPairwise(Scelltype_KO)
    # Scelltype_KO = updateCellChat(Scelltype_KO)
  }else{
    subName = "WT"
    Scelltype_WT = readRDS(paste0(cellchatDir,subName,"_",type,"_allImmuCell_CellChat_results/",LRtypes[m],"/cellchat_",subName,"_",LRtypes[m],"_",type,"_allImmuCell.rds"))  
    Scelltype_WT <- rankNetPairwise(Scelltype_WT)
    # Scelltype_WT = updateCellChat(Scelltype_WT)
  }
}

sampleoutpath=paste0(cellchatDir,"plot_",LRtypes[m],"_",type,"_allImmu/")#secreting ECM
if(!dir.exists(sampleoutpath)){dir.create(sampleoutpath)}
object.list <- list(WT = Scelltype_WT, KO = Scelltype_KO)

group.cellType = levels(Scelltype_WT@idents)
group.cellType <- factor(group.cellType, levels = unique(group.cellType))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
## Macrophage and T cell are the hub of immune regulation -- top 0.2
pdf(paste0(sampleoutpath,type,"_network_ImmuCell_top0.1_labeled",".pdf"),w=9,h=5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], #vertex.size = 10,
                   vertex.label.cex = 1.2,  edge.label.cex = 1, edge.curved = 0.2,arrow.size = 0.4,color.use = ScelltypeCols,
                   edge.width.max = 12, top=0.1,title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

## upregulated/downregulated LRpairs in KO and WT
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "KO"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, 
                                       only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 0.05)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "KO",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "WT",ligand.logFC = -0.2, receptor.logFC = NULL)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]#c(1,2,3,4)c(9:11)c(1,2,8,9,12)
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(1,2,5), targets.use = c(8,9), comparison = c(1, 2),  
                        angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 2, targets.use = c(8,9), comparison = c(1, 2),  
                        angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

#Visualize the upgulated and down-regulated signaling ligand-receptor pairs using Chord diagram
# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
pdf(paste0(sampleoutpath,"mMDSC-TAM-gMDSC_CD4-8T_upgene_chordPlot.pdf"),h=4,w=5) #allMye
netVisual_chord_gene(object.list[[1]], sources.use = c(1,2,5), targets.use = c(8,9), #sources.use = c(5), targets.use = c(seq(22,33)), #signaling =c("CCL","CXCL","XCR"), 
                     slot.name = 'net', net = net.up, color.use = ScelltypeCols,
                     lab.cex = 1, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
dev.off()
pdf(paste0(sampleoutpath,"mMDSC-TAM-gMDSC_CD4-8T_downgene_chordPlot.pdf"),h=4,w=5)
netVisual_chord_gene(object.list[[1]], sources.use = c(1,2,5), targets.use = c(8,9),
                     slot.name = 'net', net = net.down,  color.use = ScelltypeCols,
                     lab.cex = 1.2, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()

#### cell-cell contact for CD4/CD8 T cell
LR_CD4T = data.frame(pairLR.use.up[!grepl("CD8",pairLR.use.up$interaction_name),])
colnames(LR_CD4T)=colnames(pairLR.use.up)
pdf(paste0(sampleoutpath,"allMye_CD4T_upgene_chordPlot.pdf"),h=5,w=6) #allMye
netVisual_chord_gene(object.list[[1]], sources.use = c(1,2,3,4,5), targets.use = c(8), pairLR.use = LR_CD4T, 
                     slot.name = 'net', net = net.up, color.use = ScelltypeCols,
                     lab.cex = 1, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
dev.off()

LR_CD8T = data.frame(pairLR.use.up[!grepl("CD4",pairLR.use.up$interaction_name),])
colnames(LR_CD8T)=colnames(pairLR.use.up)
pdf(paste0(sampleoutpath,"allMye_CD8T_upgene_chordPlot.pdf"),h=5,w=6) #allMye
netVisual_chord_gene(object.list[[1]], sources.use = c(1,2,3,4,5), targets.use = c(9), pairLR.use = LR_CD8T, 
                     slot.name = 'net', net = net.up, color.use = ScelltypeCols,
                     lab.cex = 1, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
dev.off()


###### 单个object ############
LRtypes = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact")
m=3
type="Scelltype"
subName="ALL"
cellchatDir = "/data/jinwf/wangxf/projects/CRCprojects/seuratv3_analysis_rawdata/MDSC/mergedMDSC/mergedImmuneCell_pure/cellchat_results/"
cellchat = readRDS(paste0(cellchatDir,subName,"_",type,"_allImmuCell_CellChat_results/",LRtypes[m],"/cellchat_",subName,"_",LRtypes[m],"_",type,"_allImmuCell.rds"))

sampleoutpath=paste0(cellchatDir,"plot_",LRtypes[m],"_",type,"_allImmu/")#secreting ECM
if(!dir.exists(sampleoutpath)){dir.create(sampleoutpath)}
library(CellChat)
library(patchwork)
pdf(paste0(sampleoutpath,type,"_network_ImmuCell_top0.1_all",".pdf"),w=9,h=5)
netVisual_circle(cellchat@net$count, weight.scale = T, label.edge= F, edge.label.cex = 1, edge.curved = 0.2,arrow.size = 0.6,
                 edge.width.max = 12, top=0.1,
                 vertex.label.cex = 1.2,color.use = ScelltypeCols)
dev.off()

pdf(paste0(sampleoutpath,type,"_network_ImmuCell_top0.7_all",".pdf"),w=7,h=6)
netVisual_circle(cellchat@net$count, weight.scale = T, label.edge= F, edge.label.cex = 1, edge.curved = 0.2,arrow.size = 0.6,
                 edge.width.max = 12, top=0.7,
                 vertex.label.cex = 1.2,color.use = ScelltypeCols)
dev.off()

mat <- cellchat@net$count
shellPath=paste0(sampleoutpath,"shellPath/")
if(!dir.exists(shellPath)){dir.create(shellPath)}
for (j in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[j, ] <- mat[j, ]
  CairoPDF(paste0(shellPath,"netVisual_circle plot_", rownames(mat)[j],"_strength.pdf"),w=5,h=6)#
  netVisual_circle(mat2,  weight.scale = T, label.edge= F, edge.label.cex = 1, edge.curved = 0.2,arrow.size = 0.6,
                   edge.width.max = 8, top=1,
                   vertex.label.cex = 1.2,color.use = ScelltypeCols,
                    title.name = rownames(mat)[j])
  dev.off()
}



# Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
gg <- netAnalysis_contribution(cellchat, signaling = "CCL")
ggsave(filename=paste0(pathway_path,pathways.show, "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
par(mfrow=c(1,2))
pdf(file = paste0("allsigPathways","_netAnalysis_signalingRole_heatmap.pdf"),w=10,h=6)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",color.use = ScelltypeCols,signaling = "CCL")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",color.use = ScelltypeCols)
# print(list(ht1 ,ht2))
print(ht1)
print(ht2)
dev.off()

netVisual_aggregate(cellchat, signaling = "CCL",top=0.5, color.use = ScelltypeCols,layout = "circle")
netVisual_aggregate(cellchat, signaling = "CXCL",top=0.5, color.use = ScelltypeCols,layout = "circle")

plotGeneExpression(cellchat, signaling = "CCL", enriched.only = FALSE)
ggsave(paste0(sampleoutpath,"CCLgeneExp_all.pdf"),w=6,h=8)
plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = FALSE)
ggsave(paste0(sampleoutpath,"CXCLgeneExp_all.pdf"),w=6,h=8)

plotGeneExpression(cellchat, features = c("Ccl3","Ccl4","Ccl5","Ccl8","Ccl12","Ccr1","Ccr2","Ccr5"))
ggsave(paste0(sampleoutpath,"CCLgeneExp.pdf"),w=6,h=5)


# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
pathways.show="CCL"
pdf(paste0(sampleoutpath,"CCL_signalingRole.pdf"),w=4.5,h=3)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10, color.use = ScelltypeCols)
dev.off()
pathways.show="CXCL"
pdf(paste0(sampleoutpath,"CXCL_signalingRole.pdf"),w=4.5,h=3)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10, color.use = ScelltypeCols)
dev.off()

source("/home/wangxf/data/projects/CRCprojects/scripts/mynetAnalysis_signalingRole_scatter.R")
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat, color.use = ScelltypeCols)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- mynetAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"),label.size = 4,font.size = 12,color.use = ScelltypeCols)
gg2
ggsave(paste0(sampleoutpath,"CCL-CXCL_signalingRole_scatter.pdf"),w=4.5,h=3.8)

gg2 <- mynetAnalysis_signalingRole_scatter(cellchat, 
                                           # signaling = c("MHC-II",'CD86'),
                                           # signaling = c("MHC-I","MHC-II",'CD86','CD80'),
                            label.size = 4,font.size = 12,color.use = ScelltypeCols)
gg2
ggsave(paste0(sampleoutpath,"allCell-cellContact_signalingRole_scatter.pdf"),w=4.5,h=3.8)

cellchat@netP$pathways #查看所有的信号通路


