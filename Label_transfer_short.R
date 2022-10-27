require(devtools)
devtools::install_github("immunogenomics/harmony", force=TRUE)
devtools::install_github("powellgenomicslab/scPred")
library(scPred)
library(magrittr)
library(Seurat)
library(dplyr)
library(tibble)

# expression data is downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110499

counts=read.table("pool_macro_DEG/matrix_s_raw_ortholog.csv",header=T,sep=",",as.is=T,row.names=1)
data <- CreateSeuratObject(counts = counts, project = "ref_T", min.cells = 5, min.features = 5)
data <- FindVariableFeatures(data, do.plot = T, nfeatures =20000,selection_method= "disp" )
data <- ScaleData(data)

lineage="T"
ref.small <- readRDS(file=paste0("/mnt/533ee9c3-18c0-4c72-a09e-d9ce5a10ef9e/T_anno_ref_total/BrCa_",lineage,"_scPred_training.rds"),
	                             refhook = NULL)
get_scpred(ref.small)
plot_probabilities(ref.small)


data <- scPredict(data, ref.small,   threshold = 0.4)
data <- RunPCA(data)
data <- RunUMAP(data, dims = 1:5)
DimPlot(data, group.by="scpred_prediction", reduction = "umap")
scpred_predicted_id <- subset.data.frame(data.frame( data@meta.data), select = c(scpred_prediction))
write.csv(scpred_predicted_id,paste0("pool_",lineage,"_scpred_subtype.csv"),row.names=T)

ref.small <- subset(ref.small, downsample = 1000)

write.table(as.matrix(GetAssayData(object = ref.small, slot = "counts")), 
            paste0("/mnt/533ee9c3-18c0-4c72-a09e-d9ce5a10ef9e/T_anno_ref_total/BrCa_",lineage,"_scPred_training_mtx.csv"), 
            sep = ',', row.names = T, col.names = T, quote = F)
