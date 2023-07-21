# Kisme

library(dplyr)

library(Seurat)

library(patchwork)

library(SoupX)

#remove the ambient RNA form the cells using SoupX for each sample

sc= load10X('D:/xxxx_gh38/..../file')

sc = autoEstCont(sc)

out = adjustCounts(sc)

pbmc <- CreateSeuratObject(counts = out, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")

plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 20)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 3000)

top10 <- head(VariableFeatures(pbmc), 20)

plot1 <- VariableFeaturePlot(pbmc)

plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

plot1 + plot2

all.genes <- rownames(pbmc)

pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:20)

pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc <- RunUMAP(pbmc, dims = 1:20)

DimPlot(pbmc, reduction = "umap")

## generate or use the nf2 file for each sample

library(DropletQC)

nf2 <- nuclear_fraction_annotation(
  annotation_path = "/data/..../humangenome_genes.gtf",
  bam ="/data/..../possorted_genome_bam.bam",
  barcodes = "/data/..../filtered_feature_bc_matrix/barcodes.tsv.gz",
  tiles = 1, cores = 1, verbose = FALSE)

saveRDS(nf2, file = "/data/..../FBxxxx_nf2.rds")


nf2=readRDS("/data/..../FBxxxx_nf2.rds")

pbmc<-AddMetaData(pbmc, nf2)


input=pbmc@meta.data[,c("nuclear_fraction","nCount_RNA")]

output=identify_empty_drops(input)

pbmc<-AddMetaData(pbmc, output)

DimPlot(pbmc, group.by = "cell_status")

saveRDS(pbmc, file = "E:/cellranger3.1_gh38/FB20195_SoupX_nf2_200gene.rds")


#Merge these samples together

library(dplyr)

library(Seurat)

library(patchwork)

library(cowplot)

PCW8 = readRDS("/data/..../FB15_SoupX_nf2.rds")

PCW9 = readRDS("/data/..../FB13_SoupX_nf2.rds")

PCW12 = readRDS("/data/..../FB20195_SoupX_nf2.rds")

PCW13 = readRDS("/data/..../FB20197_SoupX_nf2.rds")

PCW14 = readRDS("/data/..../FB20192_SoupX_nf2.rds")

PCW15 = readRDS("/data/..../FB20191_SoupX_nf2.rds")

PCW16 = readRDS("/data/..../FB20193_SoupX_nf2.rds")

PCW17 = readRDS("/data/..../FB20194_SoupX_nf2.rds")

pbmc = merge(x = PCW8, y = c(PCW9, PCW12, PCW13, PCW14, PCW15, PCW16, PCW17))

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 3000)

all.genes <- rownames(pbmc)

pbmc <- ScaleData(pbmc, features = all.genes)


pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


ElbowPlot(pbmc)


BATCH=c(rep('PCW8',ncol(PCW8)),
        rep('PCW9',ncol(PCW9)),
        rep('PCW12',ncol(PCW12)),
        rep('PCW13',ncol(PCW13)),
        rep('PCW14',ncol(PCW14)),
        rep('PCW15',ncol(PCW15)),
        rep('PCW16',ncol(PCW16)), 
        rep('PCW17',ncol(PCW17)))

pbmc$batch=BATCH


pbmc <- FindNeighbors(pbmc, dims = 1:23)

pbmc <- FindClusters(pbmc, resolution = 0.8)

pbmc <- RunUMAP(pbmc, dims = 1:23)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters

DimPlot(pbmc, reduction = "umap")

DimPlot(pbmc, group.by = "cell_status")

saveRDS(pbmc, file = "/data/xxxxx/pbmc_SoupX_nf2.rds")
