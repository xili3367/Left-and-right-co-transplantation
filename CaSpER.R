require(devtools)
install_github("akdess/CaSpER")
devtools::install("BiomaRt")
BiocManager::install("maftools")
BiocManager::install("org.Mm.eg.db")

library(sigminer)
library(biomaRt)
library(Seurat)
library(CaSpER)
library(sigminer)
library(org.Mm.eg.db)
library(clusterProfiler)
library(pheatmap)

getwd()
setwd("/mnt/c9b6130c-37e5-4f62-becc-dd4240b42021/T127_T22_ola_2")

# expression data is downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110499
outputdir="pool_reconcat_DEG/"
GEX_ID="coGEX6"
counts=read.table(paste0(outputdir,GEX_ID,"_matrix","_raw.csv"),header=T,sep=",",as.is=T,row.names=1)
counts$X=row.names(counts)
anno<-bitr(row.names(counts), fromType = "SYMBOL", toType = "ENSEMBL",OrgDb = "org.Mm.eg.db")
X<-merge(counts,anno,by.x="X" , by.y = "SYMBOL")
X<- subset(X, select = -c(X))
row.names(X)=make.names(X$ENSEMBL, unique=T)
X<- subset(X, select = -c(ENSEMBL))
meta <- read.table(file=paste0(outputdir,GEX_ID,"_meta","_anno_clusters.csv"),sep=",", as.is=TRUE,
                   row.names=1,header=TRUE)
CellID <- row.names(meta)
meta <- subset.data.frame(meta, select = c(anno_clusters))
counts <- counts[,CellID]
  
data <- CreateSeuratObject(counts = X, project = GEX_ID, min.cells = 0, min.features = 0)
#data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
#data <- subset(data, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 40)
data <- NormalizeData(data , scale.factor = 1e6, normalization.method = "RC")
data <- FindVariableFeatures(data, do.plot = T, nfeatures = 20000)
data <- ScaleData(data)

data <- RunPCA(data, features = VariableFeatures(object = data),npcs = 100)
data <- RunTSNE(data, dims.use = 1:10)
DimPlot(data, reduction = "tsne")

data <- FindNeighbors(data, dims = 1:10)
data <- FindClusters(data, resolution = 0.5)
DimPlot(data, reduction = "tsne", label=T)
Idents(data)<- meta$anno_clusters

log.ge <- as.matrix(data@assays$RNA@data)
control <- names(Idents(data) )[Idents(data) %in% c('NK.cells', 'Tgd', 'Tprf.MKI67','B.cells', 'CD14.IL1B', 'CD14.PF4', 'CD16.CDKN1C', 'CD4_Tcm.TCF7',
       'CD4_Treg.FOXP3', 'CD8_Tn.CCR7')]

genes <- rownames(log.ge)
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=genes)
log.ge <- log.ge[match( annotation$Gene,rownames(log.ge)) , ]
rownames(log.ge) <- annotation$Gene
log.ge <- log2(log.ge +1)


loh <-readBAFExtractOutput ( path=paste0("./baf/",GEX_ID), sequencing.type="bulk", suffix=".snp")
names(loh) <- gsub(".snp", "", names(loh))
loh.name.mapping <- data.frame (loh.name= "baf" , sample.name=colnames(log.ge))


object <- CreateCasperObject(raw.data=log.ge,loh.name.mapping=loh.name.mapping, sequencing.type="single-cell", 
                             cnv.scale=3, loh.scale=3, genomeVersion = "mm10",
                             expr.cutoff=0.1, filter="median", matrix.type="normalized",
                             annotation=annotation, method="iterative", loh=loh, 
                             control.sample.ids=control, cytoband=cytoband)


## runCaSpER

object@annotation$cytoband <- object@annotation$Chr
object@annotation.filt$cytoband <- object@annotation.filt$Chr
final.objects <- runCaSpER(object)


## summarize large scale events 
 
final.objects <- readRDS("coGEX6_CNV.rds")
obj <- final.objects[[9]]

annotation_row=meta
ann_colors = list(
  anno_clusters= c(B.cells = "#8b008b",CD14.IL1B = "#00ff7f" ,CD14.PF4 = "#20b2aa" ,CD16.CDKN1C = "#6b8e23" ,
CD4_Tcm.TCF7 = "#da70d6" ,CD4_Treg.FOXP3 ="#f08080" ,CD8_Tn.CCR7 = "#b22222" ,Endo.cells = "#d2b48c" ,
Macro.M1 = "#00ffff" ,Macro.M2 = "#5f9ea0",Macro.MGP = "#4169e1" ,Macro.MKI67 = "#191970" ,
Malign_0 = "#b0c4de" ,Malign_1 = "#696969" ,Malign_2 = "#708090" ,Malign_3 = "#000000" ,
NK.cells = "#cd853f" ,Tgd = "#dc143c" ,Tprf.MKI67 = "#ff7f50" ,i_CAF = "#bc8f8f" ,my_CAF = "#daa520" ,
pDC = "#00bfff" )
)

obj_order=obj@control.normalized.noiseRemoved[[3]]
meta=cbind(meta,meta)
order<-meta[order(meta$anno_clusters),]
obj_order=obj_order[,row.names(order)]
obj@control.normalized.noiseRemoved[[3]]=obj_order
write.csv(obj_order, "coGEX6_cnv2.csv")
#saveRDS(final.objects, file = "coGEX6_CNV.rds")
p<-plotHeatmap10x(object=obj,
               annotation_row=annotation_row,
               annotation_colors = ann_colors,
               fileName="CaSpER_heatmap_coGEX6_3.png",cnv.scale= 3,
               cluster_cols = F, cluster_rows = F, show_rownames =F, only_soi = F)

#### VISUALIZATION 
finalChrMat <- extractLargeScaleEventsMouse (final.objects, thr=0.75)
chrMat <- finalChrMat
chrMat <- chrMat[row.names(order),]
plot.data <- melt(chrMat)
plot.data$value2 <- "neutral"
plot.data$value2[plot.data$value > 0] <- "amplification"
plot.data$value2[plot.data$value < 0] <- "deletion"
plot.data$value2 <- factor(plot.data$value2, levels = c("amplification", 
                                                        "deletion", "neutral"))
plot.data$X2 <- factor(plot.data$X2, levels = colnames(chrMat))

write.csv(chrMat, "coGEX6_del_amp2.csv")
tiff(filename="CaSpER_coGEX6_del_amp2.tiff", width=1800, height=2000, res=300)
ggplot(aes(x = X2, y = X1,fill = value2), data = plot.data) + 
  geom_tile(colour = "white", size = 0.01) + 
  labs(x = "", 
       y = "") + scale_fill_manual(values = c(amplification = muted("red"), 
                                              deletion = muted("green"), neutral = "white")) + theme_grey(base_size = 6) + 
  theme(legend.position = "right", legend.direction = "vertical", 
        legend.title = element_blank(), strip.text.x = element_blank(), 
        legend.text = element_text(colour = "black", size = 7, 
                                   face = "bold"), legend.key.height = grid::unit(0.8, 
                                                                                  "cm"), legend.key.width = grid::unit(0.5, "cm"), 
        axis.text.x = element_text(size = 5, colour = "black", 
                                  angle = -45, hjust = 0), axis.text.y = element_blank(), 
        plot.title = element_text(colour = "black", hjust = 0, 
                                  size = 6, face = "bold"))

dev.off()

