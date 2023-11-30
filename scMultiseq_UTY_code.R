# Analysis of UTY cardiac leukocyte scMultiomics dataset
# Nick Chavkin
# Walsh lab

# Load Modules
library(tidyverse)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)

# load the RNA and ATAC data
counts <- Read10X_h5("tac_aggr/outs/filtered_feature_bc_matrix.h5")
fragpath <- "tac_aggr/outs/atac_fragments.tsv.gz"

# get gene annotations for mm10
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "GRCm38"


# create a Seurat object containing the RNA data
tac_so <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
tac_so[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
tac_so # 17325 cells with 176007 features in 2 assays

DefaultAssay(tac_so) <- "ATAC"

tac_so <- NucleosomeSignal(tac_so)
tac_so <- TSSEnrichment(tac_so)

VlnPlot(
  object = tac_so,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 1
)

# filter out low quality cells
tac_so <- subset(
  x = tac_so,
  subset = nCount_ATAC < 60000 &
    nCount_RNA < 1500 &
    nCount_ATAC > 500 &
    nCount_RNA > 50 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
tac_so # 16341 cells

saveRDS(tac_so, "tac_so.rds")
tac_so <- readRDS("tac_so.rds")

# call peaks using MACS2
peaks <- CallPeaks(tac_so, macs2.path = "~/opt/miniconda3/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(tac_so),
  features = peaks,
  cells = colnames(tac_so)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
tac_so[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

# Gene Expression Processing
DefaultAssay(tac_so) <- "RNA"
tac_so <- SCTransform(tac_so)
tac_so <- RunPCA(tac_so)


# DNA Accessibility Processing
DefaultAssay(tac_so) <- "peaks"
tac_so <- FindTopFeatures(tac_so, min.cutoff = 5)
tac_so <- RunTFIDF(tac_so)
tac_so <- RunSVD(tac_so)

# UMAP on ATAC data only
DefaultAssay(tac_so) <- "peaks"
tac_so <- RunUMAP(tac_so, reduction = 'lsi', dims = 2:30)
tac_so <- FindNeighbors(tac_so, reduction = 'lsi', dims = 2:30)
tac_so <- FindClusters(tac_so, verbose = FALSE, algorithm = 3)
DimPlot(tac_so, label = TRUE) + NoLegend()
DimPlot(tac_so, split.by = "orig.ident") + NoLegend()
tac_so$orig.ident <- substr(Cells(tac_so), 18, 19)
tac_so$celltype.condition <- paste(tac_so$seurat_clusters, tac_so$orig.ident, sep = "_")
tac_so$phenotype.condition <- paste(tac_so$celltype, tac_so$orig.ident, sep = "_")

# Gene Activity Matrix
gene.activities <- GeneActivity(tac_so)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
tac_so[['GAM']] <- CreateAssayObject(counts = gene.activities)
tac_so <- NormalizeData(
  object = tac_so,
  assay = 'GAM',
  normalization.method = 'LogNormalize',
  scale.factor = median(tac_so$nCount_RNA)
)

saveRDS(tac_so, "tac_so.rds")
tac_so <- readRDS("tac_so.rds")

FeaturePlot(
  object = tac_so,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

# MAGIC imputation
library(Rmagic)
tac_so.sct <- as.data.frame(tac_so$SCT@scale.data)
write.csv(tac_so.sct, "tac_so_sct.csv")
tac_so.magicdata <- as.data.frame(t(tac_so.sct))
write.csv(tac_so.magicdata, "tac_so_magicdata.csv")
magic.genes <- row.names(tac_so.sct)
tac_so.magic <- magic(tac_so.magicdata, t = 10, genes = magic.genes)
tac_so.magicresults <- t(as.data.frame(tac_so.magic$result))
tac_so[['MAGIC']] <- CreateAssayObject(counts = tac_so.magicresults)

DefaultAssay(tac_so) <- 'RNA'
tac_so <- ScaleData(tac_so)
tac_so.rna <- as.data.frame(tac_so$RNA@scale.data)
write.csv(tac_so.rna, "tac_so_rna.csv")
tac_so.magicdatarna <- as.data.frame(t(tac_so.rna))
magic.genesrna <- row.names(tac_so.rna)
tac_so.magicrna <- magic(tac_so.magicdatarna, t = 10, genes = magic.genesrna)
tac_so.magicresultsrna <- t(as.data.frame(tac_so.magicrna$result))
tac_so[['MAGICRNA']] <- CreateAssayObject(counts = tac_so.magicresultsrna)
rm(tac_so.magic, tac_so.magicresults, tac_so.magicresultsrna, tac_so.magicrna)
saveRDS(tac_so, "tac_so.rds")
saveRDS(tac_so_motifs, "tac_so_motifs.rds")

# Gene expression by MAGIC
DefaultAssay(tac_so) <- 'MAGIC'
FeaturePlot(tac_so, features = c("Il1b", "Mrc1"))
VlnPlot(tac_so, features = c("Il1b", "Mrc1"), assay = 'MAGIC', pt.size = 0)
DimPlot(tac_so, label = TRUE)

DefaultAssay(tac_so) <- 'MAGICRNA'
tac_so <- SetIdent(tac_so, value = "seurat_clusters")
VlnPlot(tac_so, features = c("Lgals3"), split.by = "orig.ident", group.by = "celltype", assay = 'MAGICRNA', pt.size = 0)
VlnPlot(tac_so, features = c("Lgals3"), split.by = "orig.ident", group.by = "celltype", assay = 'GAM', pt.size = 0)

DefaultAssay(tac_so) <- 'MAGICRNA'
DimPlot(tac_so, label = TRUE) + NoLegend()
VlnPlot(tac_so, features = "Cd68", pt.size = 0) + NoLegend() 
# Mono/Mac/Neut: 1, 4, 5, 7, 8, 9, 10, 14, 16, 18, 19, 20, 21
VlnPlot(tac_so, features = "Cd79a", pt.size = 0) + NoLegend() 
# B cells: 0, 2, 11, 18, 21
DefaultAssay(tac_so) <- 'RNA'
FeaturePlot(tac_so, features = "Csf3r")
VlnPlot(tac_so, features = "Csf3r", pt.size = 0) + NoLegend() 
# Neut: 
VlnPlot(tac_so, features = "Cd3d", pt.size = 0) + NoLegend() 
# T cells: 12, 15
VlnPlot(tac_so, features = "Klrk1", pt.size = 0) + NoLegend() 
# NK cells: 13

DefaultAssay(tac_so) <- 'MAGICRNA'
DotPlot(tac_so, features = c("Cd68", "Ly6c2", "Il1b", "Ccr2", "Mrc1", "Lyve1", "Cd79a", "Cd3d", "Klrk1", "Csf3r", "Cd209a", "Mki67"), group.by = "celltype")
FeaturePlot(tac_so, features = c("Cd68", "Ly6c2", "Ccr2", "Cd79a", "Cd3d", "Klrk1", "Csf3r", "Cd209a"), ncol = 4, min.cutoff = 'q05', max.cutoff = 'q95')
FeaturePlot(tac_so, features = "Mki67", min.cutoff = 'q05', max.cutoff = 'q95')

DimPlot(tac_so)

# Annotation of UMAP clusters
tac_so <- SetIdent(tac_so, value = "seurat_clusters")
DimPlot(tac_so, label = TRUE)
tac_so <-RenameIdents(tac_so,
                      "0" = "B Cells",
                      "1" = "Macrophages",
                      "2" = "B Cells",
                      "3" = "B Cells",
                      "4" = "Macrophages",
                      "5" = "Neutrophils",
                      "6" = "Macrophages",
                      "7" = "Macrophages",
                      "8" = "Monocytes",
                      "9" = "Macrophages",
                      "10" = "Monocytes",
                      "11" = "B Cells",
                      "12" = "T Cells",
                      "13" = "NK Cells",
                      "14" = "Proliferative Cells",
                      "15" = "T Cells",
                      "16" = "Macrophages",
                      "17" = "Unidentified",
                      "18" = "Unidentified",
                      "19" = "Unidentified",
                      "20" = "Basophils",
                      "21" = "Unidentified")
tac_so[["leukocyte"]] <- Idents(tac_so)

# Load the pre-processed scRNA-seq data for PBMCs
mLOY_rna <- readRDS("mLOY2_so.rds")

DefaultAssay(tac_so) <- 'GAM'
DefaultAssay(mLOY_rna) <- 'integrated'

transfer.anchors <- FindTransferAnchors(
  reference = mLOY_rna,
  query = tac_so,
  reduction = 'cca',
  normalization.method = 'SCT'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = mLOY_rna$seurat_clusters,
  weight.reduction = tac_so[['lsi']],
  dims = 2:30
)

tac_so <- SetIdent(tac_so, value = "celltype.condition")

FeaturePlot(tac_so, features = "prediction.score.2", split.by = "orig.ident")
FeaturePlot(tac_so, features = "prediction.score.3", split.by = "orig.ident")

tac_so <- AddMetaData(object = tac_so, metadata = predicted.labels)

tac_so <- RenameIdents(tac_so,
                       "Monocytes" = "Monocytes",
                       "Macrophages" = "Macrophages",
                       "Neutrophils" = "Neutrophils",
                       "Proliferative Cells" = "Proliferative Cells",
                       "B Cells" = "B Cells",
                       "T Cells" = "T Cells",
                       "NK Cells" = "NK Cells",
                       "Basophils" = "Basophils",
                       "Unidentified" = "Unidentified")

tac_so$celltype <- tac_so@active.ident
tac_so$celltype.condition <- paste(tac_so$seurat_clusters, tac_so$orig.ident, sep = "_")

DimPlot(tac_so, group.by = "seurat_clusters", label = TRUE)
DimPlot(tac_so, label = TRUE)
DimPlot(tac_so, label = FALSE, split.by = "orig.ident")

# significantly enriched peaks
DefaultAssay(tac_so) <- 'peaks'
tac_so <- SetIdent(tac_so, value = "orig.ident")

da_peaks <- FindMarkers(
  object = tac_so,
  ident.1 = "1",
  ident.2 = "2",
  min.pct = 0.05,
  test.use = 'LR'
)

# Motif enrichment
# Get a list of motif position frequency matrices from the JASPAR database
library(JASPAR2020)
library(TFBSTools)
library(patchwork)

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
DefaultAssay(tac_so) <- 'peaks'
# issue solved: https://github.com/timoast/signac/issues/486
main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
keep.peaks <- as.logical(seqnames(granges(tac_so)) %in% main.chroms)
tac_so_motifs <- tac_so[keep.peaks, ]
tac_so_motifs <- AddMotifs(
  object = tac_so_motifs,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)



top_uty_peaks_fibrotic <- rownames(markers_fibrotic_uty_peaks[markers_fibrotic_uty_peaks$p_val < 0.005 & markers_fibrotic_uty_peaks$avg_log2FC > 0, ])

# test enrichment
tac_so_motifs <- SetIdent(tac_so_motifs, value = "seurat_clusters")
open.peaks <- AccessiblePeaks(tac_so_motifs, idents = c("6", "7", "9", "16"))
meta.feature <- GetAssayData(tac_so_motifs, assay = 'peaks', slot = 'meta.features')
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top_uty_peaks_fibrotic, ],
  n = 50000
)
enriched.motifs <- FindMotifs(
  object = tac_so_motifs,
  features = top_uty_peaks_fibrotic
)
write.csv(enriched.motifs, "motifs_uty_fibrotic_up.csv")

tac_so <- SetIdent(tac_so, value = "celltype.condition")
markers_monomac_uty_GAM_all <- FindMarkers(
  tac_so,
  assay = "GAM",
  ident.1 = c("1_1", "4_1","6_1", "7_1", "8_1", "9_1", "10_1", "16_1"),
  ident.2 = c("1_2", "4_2","6_2", "7_2", "8_2", "9_2", "10_2", "16_2"),
  test.use = 'LR',
  min.pct = '0.05',
  logfc.threshold = '0.00'
)
write.csv(markers_monomac_uty_GAM_all, "markers_monomac_uty_GAM_all.csv")

top_uty_peaks_monomac <- rownames(markers_monomac_uty_peaks[markers_monomac_uty_peaks$p_val < 0.005 & markers_monomac_uty_peaks$avg_log2FC > 0, ])

bottom_uty_peaks_monomac <- rownames(markers_monomac_uty_peaks[markers_monomac_uty_peaks$p_val < 0.005 & markers_monomac_uty_peaks$avg_log2FC < 0, ])

# test enrichment
tac_so_motifs <- SetIdent(tac_so_motifs, value = "seurat_clusters")
open.peaks <- AccessiblePeaks(tac_so_motifs, idents = c("1", "4", "6", "7", "8", "9", "10", "16"))
meta.feature <- GetAssayData(tac_so_motifs, assay = 'peaks', slot = 'meta.features')
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top_uty_peaks_monomac, ],
  n = 50000 
)
enriched.motifs_monomac <- FindMotifs(
  object = tac_so_motifs,
  features = top_uty_peaks_monomac
)
write.csv(enriched.motifs_monomac, "motifs_uty_monomac_up.csv")

enriched.motifs_monomac_down <- FindMotifs(
  object = tac_so_motifs,
  features = bottom_uty_peaks_monomac
)
write.csv(enriched.motifs_monomac_down, "motifs_uty_monomac_down.csv")

MotifPlot(
  object = tac_so_motifs,
  motifs = head(rownames(enriched.motifs_monomac))
)

tac_so_motifs <- RunChromVAR(
  object = tac_so_motifs,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

saveRDS(tac_so_motifs, "tac_so_motifs.rds")

markers_chromvar_uty <- FindMarkers(
  tac_so_motifs_monomac,
  assay = 'chromvar',
  ident.1 = "1",
  ident.2 = "2"
)

markers_monomac_uty_peaks <- FindMarkers(
  tac_so,
  assay = "peaks",
  ident.1 = c("1_1", "4_1", "6_1", "7_1", "8_1", "9_1", "10_1", "16_1"),
  ident.2 = c("1_2", "4_2", "6_2", "7_2", "8_2", "9_2", "10_2", "16_2"),
  min.pct = '0.05',
  test.use = 'LR',
  logfc.threshold = '0.05'
) # macrophage clusters are 1, 4, 6, 7, 9, 16, monocytes are 8, 10


