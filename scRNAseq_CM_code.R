# Analysis of publicly availabile Cardiomyopathy scRNAseq datasets
# Nick Chavkin
# Walsh lab

### Seurat Object Generation ###
Hil_meta <- read.csv("/project/walsh_rivanna_paid/cm/Hill/GSE203274/suppl/GSE203274_AllNuclei_snRNA_metadata.csv", header = TRUE)
Hil_counts <- Read10X(data.dir = "/project/walsh_rivanna_paid/cm/Hill/GSE203274/suppl/AllNuclei_snRNA_rawCount/", gene.column = 1)
Hil_so <- CreateSeuratObject(counts = Hil_counts, meta.data = Hil_meta, project = "Hil")
saveRDS(Hil_so, "/project/walsh_rivanna_paid/cm/Analysis/Hil_so.rds")

DCM1 <- Read10X("/project/walsh_rivanna_paid/cm/Mehdiabadi/GSE185100/suppl/DCM1")
DCM2 <- Read10X(data.dir = "/project/walsh_rivanna_paid/cm/Mehdiabadi/GSE185100/suppl/DCM2")
DCM3 <- Read10X(data.dir = "/project/walsh_rivanna_paid/cm/Mehdiabadi/GSE185100/suppl/DCM3")
DCM4 <- Read10X(data.dir = "/project/walsh_rivanna_paid/cm/Mehdiabadi/GSE185100/suppl/DCM4")
Fetal1 <- Read10X(data.dir = "/project/walsh_rivanna_paid/cm/Mehdiabadi/GSE185100/suppl/Fetal1")
Fetal2 <- Read10X(data.dir = "/project/walsh_rivanna_paid/cm/Mehdiabadi/GSE185100/suppl/Fetal2")
Fetal3 <- Read10X(data.dir = "/project/walsh_rivanna_paid/cm/Mehdiabadi/GSE185100/suppl/Fetal3")
Young1 <- Read10X(data.dir = "/project/walsh_rivanna_paid/cm/Mehdiabadi/GSE185100/suppl/Young1")
Young2 <- Read10X(data.dir = "/project/walsh_rivanna_paid/cm/Mehdiabadi/GSE185100/suppl/Young2")
Young3 <- Read10X(data.dir = "/project/walsh_rivanna_paid/cm/Mehdiabadi/GSE185100/suppl/Young3")
DCM1 <- CreateSeuratObject(DCM1)
DCM2 <- CreateSeuratObject(DCM2)
DCM3 <- CreateSeuratObject(DCM3)
DCM4 <- CreateSeuratObject(DCM4)
Fetal1 <- CreateSeuratObject(Fetal1)
Fetal2 <- CreateSeuratObject(Fetal2)
Fetal3 <- CreateSeuratObject(Fetal3)
Young1 <- CreateSeuratObject(Young1)
Young2 <- CreateSeuratObject(Young2)
Young3 <- CreateSeuratObject(Young3)
Meh_so <- merge(x = DCM1, y = list(DCM2,DCM3,DCM4,Fetal1,Fetal2,Fetal3,Young1,Young2,Young3))
saveRDS(Meh_so, "/project/walsh_rivanna_paid/cm/Analysis/Meh_so.rds")
rm(DCM1,DCM2,DCM3,DCM4,Fetal1,Fetal2,Fetal3,Young1,Young2,Young3)

Koe_so <- load("/project/walsh_rivanna_paid/cm/Koenig/GSE183852/suppl/GSE183852_DCM_Integrated.Robj")
saveRDS(Koe_so, "/project/walsh_rivanna_paid/cm/Analysis/Koe_so.rds")

### LOY determination ###
Koe_so <- readRDS("/project/walsh_rivanna_paid/cm/Analysis/Koe_so.rds")
Koe_so # 269,794 cells
Koe_m_so <- subset(Koe_so, subset = Sex == "Male")
Koe_m_so <- subset(Koe_m_so, subset = nCount_RNA > 2000)
Koe_m_so # 94,328
age <- data.frame(orig.ident = Koe_m_so$orig.ident)
Koe_age <- read.csv("Koe_age.csv")
age <- left_join(age, Koe_age)
Ygenes <- read.csv("Ygenes.csv", header = TRUE)
Koe_m_so[["percent.y"]] <- PercentageFeatureSet(Koe_m_so, features = intersect(Ygenes$external_gene_name, rownames(Koe_m_so)))
Koe_m_so[["LOY"]] <- ifelse(Koe_m_so$percent.y == 0, "LOY", "WT")
Koe_meta <- data.frame(orig.ident = Koe_m_so$orig.ident,
                       nCount_RNA = Koe_m_so$nCount_RNA,
                       nFeature_RNA = Koe_m_so$nFeature_RNA,
                       age = age$age,
                       sex = Koe_m_so$Sex,
                       condition = Koe_m_so$condition,
                       percent.y = Koe_m_so$percent.y,
                       LOY = Koe_m_so$LOY)
Koe_meta <- mutate(Koe_meta, condition = case_when(condition == "Donor" ~ "Healthy",
                                                   condition == "DCM" ~ "DCM"))
Koe_so.counts <- GetAssayData(Koe_m_so, slot = "counts")
Koe_c_so <- CreateSeuratObject(counts = Koe_so.counts)
Koe_c_so # 94,328
Koe_c_so <- AddMetaData(Koe_c_so, Koe_meta)

Ygenes <- read.csv("Ygenes.csv", header = TRUE)
Koe_m_so[["percent.y"]] <- PercentageFeatureSet(Koe_m_so, features = intersect(Ygenes$external_gene_name, rownames(Koe_m_so)))
Koe_m_so[["LOY"]] <- ifelse(Koe_m_so$percent.y == 0, "LOY", "WT")
Koe_m_so.meta <- Koe_m_so[[]]
write.csv(Koe_m_so[[]], "Koe_m_so_meta.csv")

unique(Koe_m_so$Names)
Koe_m_so[["Leukocyte"]] <- Koe_m_so$Names
Koe_m_so[["Leukocyte"]] <- ifelse(Koe_m_so$Names == "Myeloid", "Leukocyte", Koe_m_so$Leukocyte)
Koe_m_so[["Leukocyte"]] <- ifelse(Koe_m_so$Names == "NK/T-Cells", "Leukocyte", Koe_m_so$Leukocyte)
Koe_m_so[["Leukocyte"]] <- ifelse(Koe_m_so$Names == "Mast", "Leukocyte", Koe_m_so$Leukocyte)
Koe_m_so[["Leukocyte"]] <- ifelse(Koe_m_so$Names == "B_Cells", "Leukocyte", Koe_m_so$Leukocyte)
unique(Koe_m_so$Leukocyte)

Koe_LOY_counts <- data.frame(matrix(ncol = 3, nrow = 0))
for (x in unique(Koe_m_so$orig.ident)) {
  Koe_LOY_counts <- rbind(Koe_LOY_counts, c(
    sum(Koe_m_so$LOY == "LOY" & Koe_m_so$orig.ident == x & Koe_m_so$Leukocyte == "Leukocyte"),
    sum(Koe_m_so$orig.ident == x & Koe_m_so$Leukocyte == "Leukocyte"),
    sum(Koe_m_so$LOY == "LOY" & Koe_m_so$orig.ident == x & Koe_m_so$Leukocyte == "Leukocyte")/sum(Koe_m_so$orig.ident == x & Koe_m_so$Leukocyte == "Leukocyte")
  ))
}
colnames(Koe_LOY_counts) <- c("LOY", "Total", "Percent_LOY")
rownames(Koe_LOY_counts) <- unique(Koe_m_so$orig.ident)
write.csv(Koe_LOY_counts, "Koe_LOY_counts.csv")
rm(Koe_so)

Cha_so <- readRDS("/project/walsh_rivanna_paid/cm/Analysis/Cha_so.rds")
Cha_so <- subset(Cha_so, subset = sex == "male")
Cha_so <- subset(Cha_so, subset = nCount_cellranger_raw > 2000)
Cha_so[["percent.y"]] <- PercentageFeatureSet(Cha_so, features = intersect(Ygenes$external_gene_name, rownames(Cha_so)))
Cha_so[["LOY"]] <- ifelse(Cha_so$percent.y == 0, "LOY", "WT")
Cha_meta <- data.frame(orig.ident = Cha_so$donor_id,
                       nCount_RNA = Cha_so$nCount_cellranger_raw,
                       nFeature_RNA = Cha_so$nFeature_cellranger_raw,
                       age = Cha_so$age,
                       sex = Cha_so$sex,
                       condition = Cha_so$disease,
                       percent.y = Cha_so$percent.y,
                       LOY = Cha_so$LOY)
Cha_meta <- mutate(Cha_meta, condition = case_when(condition == "NF" ~ "Healthy",
                                                   condition == "DCM" ~ "DCM",
                                                   condition == "HCM" ~ "HCM"))
Cha_so.counts <- GetAssayData(Cha_so, slot = "counts")
Cha_c_so <- CreateSeuratObject(counts = Cha_so.counts)
Cha_c_so <- AddMetaData(Cha_c_so, Cha_meta)





Cha_so <- subset(Cha_so, subset = nCount_cellranger_raw > 2500)
Cha_so # 114,149 cells
Cha_so.meta <- Cha_so[[]]
Cha_m_so <- subset(Cha_so, subset = sex == "male")
Cha_m_so[["percent.y"]] <- PercentageFeatureSet(Cha_m_so, features = intersect(Ygenes$external_gene_name, rownames(Cha_m_so)))
Cha_m_so[["LOY"]] <- ifelse(Cha_m_so$percent.y == 0, "LOY", "WT")

Cha_so <- subset(Cha_so, subset = nCount_cellranger_raw > 2500)
Cha_so # 114,149 cells
Cha_so.meta <- Cha_so[[]]
Cha_m_so <- subset(Cha_so, subset = sex == "male")
Cha_m_so[["percent.y"]] <- PercentageFeatureSet(Cha_m_so, features = intersect(Ygenes$external_gene_name, rownames(Cha_m_so)))
Cha_m_so[["LOY"]] <- ifelse(Cha_m_so$percent.y == 0, "LOY", "WT")



unique(Cha_m_so$cell_type_leiden0.6)
Cha_m_so[["Leukocyte"]] <- ifelse(Cha_m_so$cell_type_leiden0.6 == "Macrophage", "Leukocyte", Cha_m_so$cell_type_leiden0.6)
Cha_m_so[["Leukocyte"]] <- ifelse(Cha_m_so$cell_type_leiden0.6 == "Mast_cell", "Leukocyte", Cha_m_so$Leukocyte)
Cha_m_so[["Leukocyte"]] <- ifelse(Cha_m_so$cell_type_leiden0.6 == "Proliferating_macrophage", "Leukocyte", Cha_m_so$Leukocyte)
Cha_m_so[["Leukocyte"]] <- ifelse(Cha_m_so$cell_type_leiden0.6 == "Lymphocyte", "Leukocyte", Cha_m_so$Leukocyte)
unique(Cha_m_so$Leukocyte)

Cha_m_so.meta <- Cha_m_so[[]]
write.csv(Cha_m_so[[]], "Cha_so_meta.csv")

Cha_LOY_counts <- data.frame(matrix(ncol = 3, nrow = 0))
for (x in unique(Cha_m_so$donor_id)) {
  Cha_LOY_counts <- rbind(Cha_LOY_counts, c(
    sum(Cha_m_so$LOY == "LOY" & Cha_m_so$donor_id == x & Cha_m_so$Leukocyte == "Leukocyte"),
    sum(Cha_m_so$donor_id == x & Cha_m_so$Leukocyte == "Leukocyte"),
    sum(Cha_m_so$LOY == "LOY" & Cha_m_so$donor_id == x & Cha_m_so$Leukocyte == "Leukocyte")/sum(Cha_m_so$donor_id == x & Cha_m_so$Leukocyte == "Leukocyte")
  ))
}
colnames(Cha_LOY_counts) <- c("LOY", "Total", "Percent_LOY")
rownames(Cha_LOY_counts) <- unique(Cha_m_so$donor_id)
write.csv(Cha_LOY_counts, "Cha_LOY_counts.csv")

Rao_so <- readRDS("/project/walsh_rivanna_paid/cm/Analysis/Rao_so.rds")
Rao_so <- subset(Rao_so, subset = nCount_RNA > 2000)
Rao_so[["percent.y"]] <- PercentageFeatureSet(Rao_so, features = intersect(Ygenes$external_gene_name, rownames(Rao_so)))
Rao_so[["LOY"]] <- ifelse(Rao_so$percent.y == 0, "LOY", "WT")
Rao_so.counts <- GetAssayData(Rao_so, slot = "counts")
Rao_c_so <- CreateSeuratObject(counts = Rao_so.counts)
Rao_patients <- data.frame(id = as.character(c(1:24)),
                           orig.ident = c("Normal-1", "Normal-1", "Normal-1", "Normal-1",
                                          "DCM-2", "DCM-2", "DCM-2", "DCM-2",
                                          "DCM-3", "DCM-3", "DCM-3", "DCM-3",
                                          "ICM-1", "ICM-1", "ICM-1", "ICM-1",
                                          "ICM-2", "ICM-2", "ICM-2", "ICM-2",
                                          "ICM-3", "ICM-3", "ICM-3", "ICM-3"),
                           age = c(50, 50, 50, 50, 60, 60, 60, 60, 62, 62, 62, 62,
                                   40, 40, 40, 40, 54, 54, 54, 54, 45, 45, 45, 45),
                           sex = c("male", "male", "male", "male", "male", "male",
                                   "male", "male", "male", "male", "male", "male",
                                   "male", "male", "male", "male", "male", "male",
                                   "male", "male", "male", "male", "male", "male"),
                           condition = c("Healthy", "Healthy", "Healthy", "Healthy",
                                         "DCM", "DCM", "DCM", "DCM",
                                         "DCM", "DCM", "DCM", "DCM",
                                         "ICM", "ICM", "ICM", "ICM",
                                         "ICM", "ICM", "ICM", "ICM",
                                         "ICM", "ICM", "ICM", "ICM"))
Rao_meta <- Rao_so[[]]
Rao_meta <- mutate(Rao_meta, id = str_sub(row.names(Rao_meta),-2,-1))
Rao_meta$id <- str_replace(Rao_meta$id,"_","")
Rao_meta <- select(Rao_meta, -c("orig.ident"))
Rao_meta <- left_join(Rao_meta, Rao_patients)
row.names(Rao_meta) <- Cells(Rao_c_so)
Rao_meta <- select(Rao_meta, -c("id"))
Rao_c_so <- AddMetaData(Rao_c_so, Rao_meta)
Rao_c_meta <- Rao_c_so[[]]

Rao_m_so <- readRDS("/project/walsh_rivanna_paid/cm/Analysis/Rao_so.rds")
VlnPlot(Rao_m_so, features = "nCount_RNA")
Rao_m_so <- subset(Rao_m_so, subset = nCount_RNA > 2500)
Rao_m_so # 45,656 cells
Rao_m_so.meta <- Rao_m_so[[]]
Rao_m_so <- NormalizeData(Rao_m_so)
Rao_m_so <- FindVariableFeatures(Rao_m_so)
Rao_m_so <- ScaleData(Rao_m_so)
Rao_m_so <- RunPCA(Rao_m_so)
Rao_m_so <- FindNeighbors(Rao_m_so)
Rao_m_so <- FindClusters(Rao_m_so)
Rao_m_so <- RunUMAP(Rao_m_so, dims = 1:30)
DimPlot(Rao_m_so, label = TRUE)
FeaturePlot(Rao_m_so, features = "PTPRC")

Rao_m_so[["percent.y"]] <- PercentageFeatureSet(Rao_m_so, features = intersect(Ygenes$external_gene_name, rownames(Rao_m_so)))
Rao_m_so[["LOY"]] <- ifelse(Rao_m_so$percent.y == 0, "LOY", "WT")
write.csv(Rao_m_so[[]], "Rao_m_so_meta.csv")

### Dataset integration ###
library(sctransform)
cm.list <- list(Koe_m_so, Cha_m_so, Rao_m_so)
cm.list <- lapply(X = cm.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = cm.list, nfeatures = 3000)
cm.list <- PrepSCTIntegration(object.list = cm.list, anchor.features = features)
cm.anchors <- FindIntegrationAnchors(object.list = cm.list, normalization.method = "SCT",
                                     anchor.features = features)
cm_im_so <- IntegrateData(anchorset = cm.anchors, normalization.method = "SCT")
rm(cm.list)

cm_im_so <- readRDS("cm_im_so.rds")
cm_im_so <- FindVariableFeatures(cm_im_so)
cm_cm_im_soso <- ScaleData(cm_im_so)
cm_im_so <- RunPCA(cm_so)
cm_im_so <- FindNeighbors(cm_im_so)
cm_im_so <- FindClusters(cm_im_so)

saveRDS(cm_im_so, "cm_im_so.rds")

### Dataset analysis ###
library(Seurat)
library(tidyverse)
cm_im_so <- readRDS("cm_im_so.rds")
cm_im_so # 281,980 cells
cm_im_so$condition <- factor(cm_im_so$condition, levels = c("Healthy", "DCM", "HCM", "ICM"))
DefaultAssay(cm_im_so) <- "RNA"
FeaturePlot(cm_im_so, features = c("PECAM1", "PTPRC", "PDGFRA", "NKX2-5", "ACTA2"), ncol = 5)
cm_im_so <- SetIdent(cm_im_so, value = "seurat_clusters")
cm_im_so <- RenameIdents(cm_im_so, 
                         '0' = "Cardiomyocyte",
                         '1' = "Cardiomyocyte",
                         '2' = "Fibroblast",
                         '3' = "Fibroblast",
                         '4' = "Leukocyte",
                         '5' = "Leukocyte",
                         '6' = "Endothelial",
                         '7' = "Leukocyte",
                         '8' = "VSMC",
                         '9' = "Fibroblast",
                         '10' = "Leukocyte",
                         '11' = "Endothelial",
                         '12' = "Fibroblast",
                         '13' = "Endothelial",
                         '14' = "Cardiomyocyte",
                         '15' = "Endothelial",
                         '16' = "VSMC",
                         '17' = "Leukocyte",
                         '18' = "VSMC",
                         '19' = "Endothelial",
                         '20' = "Endothelial",
                         '21' = "Leukocyte",
                         '22' = "Leukocyte",
                         '23' = "Leukocyte",
                         '24' = "Fibroblast",
                         '25' = "Endothelial",
                         '26' = "Leukocyte",
                         '27' = "Unclassified",
                         '28' = "Unclassified",
                         '29' = "Unclassified",
                         '30' = "Cardiomyocyte",
                         '31' = "Endothelial",
                         '32' = "Leukocyte")
cm_im_so[["celltype"]] <- Idents(cm_im_so)
cm_im_so$celltype <- factor(cm_im_so$celltype, levels = c("Leukocyte", "Fibroblast", "Endothelial", "Cardiomyocyte", "VSMC", "Unclassified"))
saveRDS(cm_im_so, "cm_im_so.rds")

### Leukocyte analysis ###
library(Seurat)
library(tidyverse)
cm_im_so <- readRDS("cm_im_so.rds")
cm_im_so <- SetIdent(cm_im_so, value = "celltype")
cm_im_so_l <- subset(cm_im_so, idents = "Leukocyte")
DefaultAssay(cm_im_so_l) <- 'integrated'
cm_im_so_l <- FindVariableFeatures(cm_im_so_l)
cm_im_so_l <- RunPCA(cm_im_so_l)
cm_im_so_l <- FindNeighbors(cm_im_so_l)
cm_im_so_l <- FindClusters(cm_im_so_l)
cm_im_so_l <- RunUMAP(cm_im_so_l, dims = 1:30)
DotPlot(cm_im_so_l, features = c("CD68", "CSF1R", "LYVE1", "IL1B", "CCR2", "MKI67", "CD79A", "CD3G", "NCAM1", "CD80"), group.by = "leukocyte", assay = "RNA")
cm_im_so_l <- SetIdent(cm_im_so_l, value = "seurat_clusters")
cm_im_so_l <- RenameIdents(cm_im_so_l, 
                           '0' = "Tissue-Resident Macrophage",
                           '1' = "T Cell",
                           '2' = "Tissue-Resident Macrophage",
                           '3' = "Macrophage",
                           '4' = "T Cell",
                           '5' = "NK Cell",
                           '6' = "Macrophage",
                           '7' = "Macrophage",
                           '8' = "Monocyte",
                           '9' = "T Cell",
                           '10' = "Neutrophil",
                           '11' = "Macrophage",
                           '12' = "Macrophage",
                           '13' = "T Cell",
                           '14' = "Macrophage",
                           '15' = "T Cell",
                           '16' = "B Cell",
                           '17' = "Macrophage",
                           '18' = "Monocyte",
                           '19' = "T Cell",
                           '20' = "B Cell",
                           '21' = "Basophil",
                           '22' = "T Cell",
                           '23' = "B Cell")
cm_im_so_l[["leukocyte"]] <- Idents(cm_im_so_l)
cm_im_so_l[["leukocyte_LOY"]] <- paste(cm_im_so_l$leukocyte, cm_im_so_l$LOY, sep = "_")
cm_im_so_l <- SetIdent(cm_im_so_l, value = "condition")
cm_im_so_l_dcm <- subset(cm_im_so_l, idents = "DCM")
cm_im_so_l_dcm # 23,734 cells
cm_im_so_l_dcm <- SetIdent(cm_im_so_l_dcm, value = "leukocyte")
DimPlot(cm_im_so_l_dcm)
DimPlot(cm_im_so_l_dcm, group.by = "LOY", cols = c("WT" = "Grey", "LOY" = "Red"))
cm_im_so_l_dcm <- SetIdent(cm_im_so_l_dcm, value = "leukocyte_LOY")
DCM_macLOY_markers <- FindMarkers(cm_im_so_l_dcm, ident.1 = "Macrophage_LOY", ident.2 = "Macrophage_WT", logfc.threshold = 0, min.pct = 0.001)

### Fibroblast analysis ###
cm_im_so <- SetIdent(cm_im_so, value = "celltype")
cm_im_so_f <- subset(cm_im_so, idents = "Fibroblast")
DefaultAssay(cm_im_so_f) <- 'integrated'
cm_im_so_f <- FindVariableFeatures(cm_im_so_f)
cm_im_so_f <- RunPCA(cm_im_so_f)
cm_im_so_f <- FindNeighbors(cm_im_so_f)
cm_im_so_f <- FindClusters(cm_im_so_f)
cm_im_so_f <- RunUMAP(cm_im_so_f, dims = 1:30)
cm_im_so_f <- SetIdent(cm_im_so_f, value = "condition")
DotPlot(cm_im_so_f, features = c("IL6", "COL1A1", "POSTN", "VEGFA"), group.by = "seurat_clusters")
cm_im_so_f <- SetIdent(cm_im_so_f, value = "seurat_clusters")
cm_im_so_f <- RenameIdents(cm_im_so_f, 
                           '0' = "Fibroblast",
                           '1' = "Fibroblast",
                           '2' = "Myofibroblast",
                           '3' = "Fibroblast",
                           '4' = "Fibroblast",
                           '5' = "Fibroblast",
                           '6' = "Angiogenic",
                           '7' = "Angiogenic",
                           '8' = "Myofibroblast",
                           '9' = "Myofibroblast",
                           '10' = "Fibroblast",
                           '11' = "Fibroblast",
                           '12' = "Angiogenic",
                           '13' = "Myofibroblast",
                           '14' = "Inflammatory")
cm_im_so_f[["celltype_fibroblast"]] <- Idents(cm_im_so_f)
DimPlot(cm_im_so_f, group.by = "celltype_fibroblast")
DotPlot(cm_im_so_f, features = c("POSTN", "ACTA2", "COL1A1", "COL3A1"), group.by = "celltype_fibroblast")
DefaultAssay(cm_im_so_f) <- "RNA"
cm_im_so_f_genes <- data.frame(POSTN = FetchData(cm_im_so_f, vars = "POSTN"), 
                               ACTA2 = FetchData(cm_im_so_f, vars = "ACTA2"), 
                               COL1A1 = FetchData(cm_im_so_f, vars = "COL1A1"),
                               COL3A1 = FetchData(cm_im_so_f, vars = "COL3A1"),
                               orig.ident = FetchData(cm_im_so_f, vars = "orig.ident"),
                               age = FetchData(cm_im_so_f, vars = "age"),
                               condition = FetchData(cm_im_so_f, vars = "condition"))