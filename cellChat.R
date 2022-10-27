require(devtools)
devtools::install_github("sqjin/CellChat")
install.packages("Seurat", repos = c("https://seurat.nygenome.org/", "https://cloud.r-project.org"))
install.packages("dplyr")

library(CellChat)
library(Seurat)
library(patchwork)
library(NMF)
library(ggalluvial)
library(ComplexHeatmap)
options(stringsAsFactors = FALSE)


# expression data is downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110499
project_list=c('T22so1',"T22so2" ,"T22so3","T22so4")


outputdir="pool_reconcat_DEG/"


counts=read.table(paste0(outputdir,"matrix_everything_s.csv"),header=T,sep=",",as.is=T,row.names=1)
meta_everything_anno <- read.table(paste0(outputdir,"meta_everything_anno_clusters.csv"), 
	                         header=T, sep=",", as.is=T, row.names=1)
meta_everything <- unclass(meta_everything_anno)
meta_everything_factor <- as.factor(meta_everything$anno_clusters)

data <- CreateSeuratObject(counts = counts, project = "all")
data<-AddMetaData(data, meta_everything_factor, col.name="anno_clusters")

meta_everything_treatment <- read.table(paste0(outputdir,"meta_everything_treatment.csv"), 
	                         header=T, sep=",", as.is=T, row.names=1)
data <- AddMetaData(data, meta_everything_treatment, col.name="treatment")
data <- NormalizeData(data , scale.factor = 1e6, normalization.method = "LogNormalize")
data@meta.data$anno_clusters = factor(data@meta.data$anno_clusters, levels =levels(data@meta.data$anno_clusters))

for (i in project_list) {
project=paste0(i,"so")


data_sub=subset(x = data, treatment==i)

data_sub <- FindVariableFeatures(data_sub, mean.function = ExpMean, dispersion.function = LogVMR,
	                          selection.method = "vst", nfeatures = 20000)
data_sub <- ScaleData(data_sub)
cellchat <- createCellChat(object = data_sub,  group.by = "anno_clusters")
#cellchat@meta$anno_clusters = factor(cellchat@meta$anno_clusters, levels =levels(cellchat@idents))
cellchat <- setIdent(cellchat, ident.use = "anno_clusters") 
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.human
#showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)

#future::plan("multiprocess", workers = 1)

cellchat <- identifyOverExpressedGenes(cellchat, thresh.pc = 0.2, thresh.fc = 0.1, thresh.p = 1)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat@idents = droplevels(cellchat@idents, exclude = setdiff(levels(cellchat@idents),unique(cellchat@idents)))
cellchat <- computeCommunProb(cellchat, raw.use=T, type="truncatedMean",
	                           trim=0.2, do.fast =  T, nboot = 24, population.size = T )
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,1), xpd=TRUE)

pdf(paste0(outputdir,project, "_","Visual_circle", ".pdf"),height =8, width = 12)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()


saveRDS (cellchat, file=paste0(outputdir,project,"_cellchat.rds"))
dev.off()
}

######################################################################################################



######################################################################################################

for ( i in c ( 'T22so1',"T22so2" ,"T22so3","T22so4")){
project_v=paste0(i, "co")
cellchat.v <- readRDS(file=paste0(outputdir,project_v,"_cellchat.rds"), refhook = NULL)
	for (k in c ('coo1',"coo2" ,"coo3","coo4")){
project_o=paste0(k,"co")
cellchat.o <- readRDS(file=paste0(outputdir,project_o,"_cellchat.rds"), refhook = NULL)


group.new=union(levels(cellchat.v@idents), levels (cellchat.o@idents))
#cellchat.v <- liftCellChat(cellchat.v, group.new=group.new)
#cellchat.o <- liftCellChat(cellchat.o, group.new=group.new)

cellchat.v <- netAnalysis_computeCentrality(cellchat.v)
cellchat.o <- netAnalysis_computeCentrality(cellchat.o)


object.list <- list(single = cellchat.v, co = cellchat.o)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = T)
#cellchat <- liftCellChat(cellchat)

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg <- gg1 + gg2
ggsave(filename=paste0(outputdir,i, "_", k,"_bar", "_comparison.pdf"), 
  	plot=gg, width = 10, height = 8, units = 'in', dpi = 300)

source_group="Malign"
list_source1 <- c()
for (m in group.new) {
   list_ <- grepl(source_group,m, fixed=T)
	 list_source1 <- append(list_source1, list_)
}
list_source2 <- c()
for (m in group.new) {
   list_ <- grepl("Macro-M2",m, fixed=T)
	 list_source2 <- append(list_source2, list_)
}
list_source3 <- c()
for (m in group.new) {
   list_ <- grepl("NK",m, fixed=T)
	 list_source3 <- append(list_source3, list_)
}


target_group1="Macro-M2"
list_target1 <- c()
for (m in group.new) {
   list_ <- grepl(target_group1,m, fixed=T)
	 list_target1 <- append(list_target1, list_)
}
target_group2="Macro-M2"
list_target2 <- c()
for (m in group.new) {
   list_ <- grepl(target_group2,m, fixed=T)
	 list_target2 <- append(list_target2, list_)
}


sources_list = group.new[list_source1]#,group.new[list_source2]),group.new[list_source3])
targets_list = group.new[list_target1]#,group.new[list_target2])



gg1<-ggplot()
cellchat <- identifyOverExpressedGenes(cellchat,  group.dataset = "datasets", 
	 pos.dataset = "co", only.pos = FALSE, thresh.pc = 0.01, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = "features")
net <- subset(net, prob > 0.0001)
net<- subset (net, pathway_name != "TGFb")
net<- subset (net, pathway_name != "COLLAGEN")
net<- subset (net, pathway_name != "SPP1")
net<- subset (net, pathway_name != "LAMININ")
net<- subset (net, pathway_name != "ICAM")
net<- subset (net, pathway_name != "JAM")
net.up <- subsetCommunication(cellchat, net = net, datasets = "co",  sources.use = sources_list, targets.use = targets_list,)
tryCatch({
gg1 <- netVisual_bubble(cellchat, pairLR.use = net.up[, "interaction_name", drop = F],
	sources.use = sources_list, targets.use = targets_list,  comparison = c(1,2),
	max.dataset = 2, title.name = "Increased signaling in co", angle.x = 45, remove.isolate = T, thresh=1)
#> Comparing communications on a merged object
 }, error=function(e){})


gg2<-ggplot()
cellchat <- identifyOverExpressedGenes(cellchat,  group.dataset = "datasets", 
	 pos.dataset = "single", only.pos = FALSE, thresh.pc = 0.01, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = "features")
net <- subset(net, prob > 0.0001)
net<- subset (net, pathway_name != "TGFb")
net<- subset (net, pathway_name != "COLLAGEN")
net<- subset (net, pathway_name != "SPP1")
net<- subset (net, pathway_name != "LAMININ")
net<- subset (net, pathway_name != "ICAM")
net<- subset (net, pathway_name != "JAM")
net.up <- subsetCommunication(cellchat, net = net, datasets = "single",  sources.use = sources_list, targets.use = targets_list,)
tryCatch({
gg2 <- netVisual_bubble(cellchat, pairLR.use = net.up[, "interaction_name", drop = F],
	sources.use = sources_list, targets.use = targets_list,  comparison = c(1,2),
	max.dataset = 1, title.name = "Decreased signaling in single", angle.x = 45, remove.isolate = T, thresh=1)
#> Comparing communications on a merged object
}, error=function(e){})
gg <- gg1 + gg2
ggsave(filename=paste0(outputdir,i, "_",k,"_Visual_bubble", "_comparison_",source_group,"_",target_group1,".pdf"), 
  	plot=gg, width = 5, height = 6, units = 'in', dpi = 300)

}
}

######################################################################################################
######################################################################################################

for ( i in c ( 'T127_V_1', 'T127_V_2')){
project_v=paste0(i, "short")
cellchat.v <- readRDS(file=paste0(outputdir,project_v,"_cellchat.rds"), refhook = NULL)
	for (k in c ('T127_O8_1',"T127_O8_2" ,"T127_O8_3")){
project_o=paste0(k,"short")
cellchat.o <- readRDS(file=paste0(outputdir,project_o,"_cellchat.rds"), refhook = NULL)


group.new=union(levels(cellchat.v@idents), levels (cellchat.o@idents))
#cellchat.v <- liftCellChat(cellchat.v, group.new=group.new)
#cellchat.o <- liftCellChat(cellchat.o, group.new=group.new)

cellchat.v <- netAnalysis_computeCentrality(cellchat.v)
cellchat.o <- netAnalysis_computeCentrality(cellchat.o)


object.list <- list(veh = cellchat.v, ola = cellchat.o)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = T)
#cellchat <- liftCellChat(cellchat)

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg <- gg1 + gg2
ggsave(filename=paste0(outputdir,i, "_", k,"_bar", "_comparison.pdf"), 
  	plot=gg, width = 10, height = 8, units = 'in', dpi = 300)

source_group="Malign"
list_source1 <- c()
for (m in group.new) {
   list_ <- grepl(source_group,m, fixed=T)
	 list_source1 <- append(list_source1, list_)
}
list_source2 <- c()
for (m in group.new) {
   list_ <- grepl("Macro-M2",m, fixed=T)
	 list_source2 <- append(list_source2, list_)
}
list_source3 <- c()
for (m in group.new) {
   list_ <- grepl("NK",m, fixed=T)
	 list_source3 <- append(list_source3, list_)
}


target_group1="Macro-M1"
list_target1 <- c()
for (m in group.new) {
   list_ <- grepl(target_group1,m, fixed=T)
	 list_target1 <- append(list_target1, list_)
}
target_group2="Macro-M2"
list_target2 <- c()
for (m in group.new) {
   list_ <- grepl(target_group2,m, fixed=T)
	 list_target2 <- append(list_target2, list_)
}


sources_list = group.new[list_source1]#,group.new[list_source2]#)#,group.new[list_source3])
targets_list = append(group.new[list_target1],group.new[list_target2])



gg1<-ggplot()
cellchat <- identifyOverExpressedGenes(cellchat,  group.dataset = "datasets", 
	 pos.dataset = "ola", only.pos = FALSE, thresh.pc = 0.01, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = "features")
net <- subset(net, prob > 0.0001)
net<- subset (net, pathway_name != "TGFb")
net<- subset (net, pathway_name != "COLLAGEN")
net<- subset (net, pathway_name != "SPP1")
net<- subset (net, pathway_name != "LAMININ")
net<- subset (net, pathway_name != "ICAM")
net<- subset (net, pathway_name != "JAM")
net.up <- subsetCommunication(cellchat, net = net, datasets = "ola",  sources.use = sources_list, targets.use = targets_list,)
tryCatch({
gg1 <- netVisual_bubble(cellchat, pairLR.use = net.up[, "interaction_name", drop = F],
	sources.use = sources_list, targets.use = targets_list,  comparison = c(1,2),
	max.dataset = 2, title.name = "Increased signaling in ola", angle.x = 45, remove.isolate = T, thresh=1)
#> Comparing communications on a merged object
 }, error=function(e){})


gg2<-ggplot()
cellchat <- identifyOverExpressedGenes(cellchat,  group.dataset = "datasets", 
	 pos.dataset = "veh", only.pos = FALSE, thresh.pc = 0.01, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = "features")
net <- subset(net, prob > 0.0001)
net<- subset (net, pathway_name != "TGFb")
net<- subset (net, pathway_name != "COLLAGEN")
net<- subset (net, pathway_name != "SPP1")
net<- subset (net, pathway_name != "LAMININ")
net<- subset (net, pathway_name != "ICAM")
net<- subset (net, pathway_name != "JAM")
net.up <- subsetCommunication(cellchat, net = net, datasets = "veh",  sources.use = sources_list, targets.use = targets_list,)
tryCatch({
gg2 <- netVisual_bubble(cellchat, pairLR.use = net.up[, "interaction_name", drop = F],
	sources.use = sources_list, targets.use = targets_list,  comparison = c(1,2),
	max.dataset = 1, title.name = "Decreased signaling in ola", angle.x = 45, remove.isolate = T, thresh=1)
#> Comparing communications on a merged object
}, error=function(e){})
gg <- gg1 + gg2
ggsave(filename=paste0(outputdir,i, "_",k,"_Visual_bubble", "_comparison_",source_group,"_",target_group1,".pdf"), 
  	plot=gg, width = 7, height = 5, units = 'in', dpi = 300)

}
}
######################################################################################################

