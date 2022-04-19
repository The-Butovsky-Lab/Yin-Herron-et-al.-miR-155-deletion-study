if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")
list.of.packages <- c('tidyverse', 'Seurat', 'Signac', 'enrichR', 'openxlsx', 'patchwork', 'data.table', 'dittoSeq', 'ggplot2','EnhancedVolcano')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

projectname = ''

options(future.globals.maxSize = 4000 * 1024^5)

# Read in the data
data1 <- Read10X(data.dir = paste0('data_files/2905_control'))
data2 <- Read10X(data.dir = paste0('data_files/2908_control'))
data3 <- Read10X(data.dir = paste0('data_files/2904_KO'))
data4 <- Read10X(data.dir = paste0('data_files/3125_KO'))


sobj1 <- CreateSeuratObject(counts = data1, project = '2905_control')
sobj2 <- CreateSeuratObject(counts = data2, project = '2908_control')
sobj3 <- CreateSeuratObject(counts = data3, project = '2904_KO')
sobj4 <- CreateSeuratObject(counts = data4, project = '3125_KO')

############################  MERGING ###########################################
scrna_merged <- merge(sobj1, y = c(sobj2, sobj3,sobj4), add.cell.ids = c('2905_control','2908_control','2904_KO','3125_KO'), project = projectname)
setable(scrna_merged$orig.ident)
table(Idents(scrna_merged))

######################### MITOCHONDRIAL REGRESSION ############################################

# Store mitochondrial gene statistics in your Seurat object
scrna_merged[['percent_mt']] <- PercentageFeatureSet(scrna_merged, pattern = '^mt-')

# Basic QC plot to set cutoffs
VlnPlot(scrna_merged, features = c('nFeature_RNA', 'nCount_RNA', 'percent_mt'), ncol = 3)
FeatureScatter(scrna_merged, feature1 = "nCount_RNA", feature2 = "percent_mt")
FeatureScatter(scrna_merged, feature1 = "nFeature_RNA", feature2 = "percent_mt")
FeatureScatter(scrna_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Filter data to remove unwanted cells
scrna_final <- subset(scrna_merged, nFeature_RNA > 1000 & nCount_RNA < 25000 & percent_mt < 20)

Idents(scrna_final) <- 'orig.ident'
scrna_final$orig.ident<- factor(x = scrna_final$orig.ident, levels = c('2905_control', '2908_control', '2904_KO','3125_KO'))
sample_counts <- table(scrna_final$orig.ident)
write.xlsx(sample_counts, file = paste0(''), overwrite = T)

dev.off()
pdf(file = paste0(''), pointsize = 10)
VlnPlot(scrna_final, features = c('nFeature_RNA', 'nCount_RNA', 'percent_mt'), ncol = 3)
dev.off()

pdf(file = paste0(''), pointsize = 10)
scat1 <- FeatureScatter(scrna_final, feature1 = 'nCount_RNA', feature2 = 'percent_mt')
scat2 <- FeatureScatter(scrna_final, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
scat1 + scat2
dev.off()

########################### NORMALIZATION AND SCALING #######################################
# Normalize RNA expression data and scale to variable features
scrna_final <- NormalizeData(scrna_final)
scrna_final <- FindVariableFeatures(scrna_final, selection.method = 'vst')
scrna_final <- ScaleData(scrna_final, features = VariableFeatures(scrna_final), vars.to.regress = 'percent_mt')


#################################### DIMENSIONAL REDUCTIONS##############################
# Run principal component analysis
scrna_final <- RunPCA(scrna_final, features = VariableFeatures(object = scrna_final))
scrna_final <- JackStraw(object = scrna_final, num.replicate = 50, prop.freq=0.025, dims = 50)
scrna_final <- ScoreJackStraw(scrna_final, dims = 1:50)


# Visualize PCA to ensure merged samples are comparable
dev.off()
Idents(scrna_final) <- 'orig.ident'
pdf(file = paste0(''), pointsize = 10)
DimPlot(scrna_final, reduction = 'pca')
dev.off()

# Visualize component strengths to decide how many to use
pdf(file = paste0(''), pointsize = 10)
JackStrawPlot <- JackStrawPlot(object = scrna_final, dims = 1:50, xmax = 0.05) + guides(col = guide_legend(ncol = 1)) + theme(legend.text = element_text(size = 6), legend.key.size = unit(0.02, "cm"))
ElbowPlot <- ElbowPlot(object = scrna_final, ndims = 50)
JackStrawPlot + ElbowPlot
dev.off()


# Cluster cells according to elbow plot dimension choice
scrna_final <- FindNeighbors(scrna_final, dims = 1:26)
scrna_final <- FindClusters(scrna_final, resolution = res)

# Run UMAP reduction to visualize clusters
scrna_final <- RunUMAP(scrna_final, dims = 1:26)

# Plot UMAP
pdf(file = paste0(''), width = 12, height = 9)
dittoDimPlot(scrna_final, var = 'seurat_clusters', 
             split.ncol = 3, size = 0.8, opacity = 0.9,
             do.label =T, labels.size = 2.5, labels.repel = F) 
dev.off()

# Save final object
saveRDS(scrna_final, file = paste0(''))

# Load Seurat data set
dataset <- readRDS(paste0('seurat_objects/SeuratObject_4samples_', dims2, res2, projectname, '.rds'))


Idents(dataset) <- 'orig.ident'
new_mappings <- c('KI','KI','KO','KO')
names(new_mappings) <- levels(dataset)
dataset <- RenameIdents(dataset, new_mappings)
dataset$Condition <- Idents(dataset)

Idents(dataset) <- 'orig.ident'
new_mappings <- c('KI1','KI2','KO1','KO2')
names(new_mappings) <- levels(dataset)
dataset <- RenameIdents(dataset, new_mappings)
dataset$Sample <- Idents(dataset)


######################## METADATA MANIPULATION AND CELLTYPE ASSIGNMENT ##########################
# Cluster mapping
Idents(dataset) <- 'seurat_clusters'
dataset_subsample <- subset(dataset, idents = c('0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15')) 
dataset <- dataset_subsample

AverageExpression<-(AverageExpression(object = dataset))
write.xlsx(AverageExpression,row.names= T,overwrite =T, '')

Idents(dataset) <- 'seurat_clusters'
new_mappings <- c('M0','intermediate','M0','MGnD','intermediate',
                  'intermediate','MGnD','M0','Proliferating','Psuedo-pro',
                  'M0','MGnD','12','MGnD','MGnD',
                  '15')
names(new_mappings) <- levels(dataset)
dataset <- RenameIdents(dataset, new_mappings)
dataset$CellType <- Idents(dataset)


Idents(dataset) <- 'CellType'
dataset_microglialcloud_subsample <- subset(dataset, idents = c('M0', 'intermediate','MGnD'))


################################## DIFFERENTIAL EXPRESSION #######################################
Idents(dataset) = 'seurat_clusters'
Idents(dataset) = 'orig.ident'
Idents(dataset) = 'Condition'
Idents(dataset) = 'CellType'
Idents(dataset) = 'Sample'

unwanted_genes <- paste(c('^mt-', '^Rp', '^Gm'), collapse = '|')

# Overall markers
sample_markers <- FindAllMarkers(dataset, 
                                 only.pos = T, 
                                 min.pct = 0.15, 
                                 logfc.threshold = 0.1) %>%
                                 filter(!str_detect(gene, unwanted_genes))
                                 top_markers <- sample_markers %>% group_by(cluster) %>% top_n(n=20, wt = avg_log2FC)

# For heatmap
heatmap_markers <- sample_markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC) %>%
  ungroup() %>% distinct(., gene, .keep_all=T) 

# For volcano plot
cluster_volcMarkers <- FindMarkers(dataset, ident.1 = '', ident.2 = '', 
                                   min.pct = 0.15, logfc.threshold = 0) %>%
                                   rownames_to_column('gene') %>% 
                                   filter(!str_detect(gene, unwanted_genes))

write.xlsx(sample_markers, '')
write.xlsx(top_markers, '')
write.xlsx(cluster_volcMarkers,'')

################################### VISUALIZATION ##############################################

#UMAP

pdf(file = '', pointsize = 10, width = 7, height = 6)
  dittoDimPlot(dataset, var = 'CellType', #split.by = 'Condition', 
             split.ncol = 2, size = 1.5, opacity = 0.9,
             do.label = F, labels.size = 2.5, labels.repel = T,color.panel=color_pallette3, show.others = FALSE,
             boxplot.show.outliers= TRUE)
dev.off()

# Violin plots
pdf(file = '', pointsize = 10, width = 5.5, height = 5)
  vln<- dittoPlot(dataset_sub, var = 'Cebpb', group.by = "Condition", 
              plots = c('boxplot'), 
              jitter.size = 0.5, split.ncol = 1)

              vln+stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", color="red", width=0.2) +
              stat_summary(fun=mean, geom="point", color="red")
dev.off()

# Feature plots
pdf(file = ''), pointsize = 10, width = 11, height =10)
  FeaturePlot(dataset, features = c('Ifitm3'), order = TRUE, pt.size=3,  min.cutoff = 0.25)
dev.off()


# Heatmap
pdf(file = ''), pointsize = 10, width = 9, height = 7)
  dittoHeatmap(dataset, sample(colnames(dataset), size = ncol(dataset)/3, replace = F)], 
             annot.by = c('Condition'), 
             order.by = 'Condition', 
             scale = 'row', 
             breaks = seq(-2, 2, by = 0.1),
             show_rownames = T,
             fontsize_row = 8,
             use_raster = T)
dev.off()

# Population Percentages Bar Plot
pdf(file = ''), pointsize = 10, width =3, height = 6)
  popProps <- dittoBarPlot(dataset, 'CellType', 
                         group.by = 'Condition', 
                         retain.factor.levels = T)
  popProps
dev.off()


# Volcano Plot
cluster_volcMarkers_file <- read.xlsx(''))


outlier_max_pval = 1.0e-150
outlier_max_logFC = 1
outlier_min_logFC = -1

cluster_volcMarkers_file <- cluster_volcMarkers_file %>%
  mutate(p_val= ifelse(p_val < outlier_max_pval, outlier_max_pval, p_val)) %>%
  mutate(avg_log2FC = ifelse(avg_log2FC < outlier_min_logFC, outlier_min_logFC, avg_log2FC)) %>%
  mutate(avg_log2FC = ifelse(avg_log2FC > outlier_max_logFC, outlier_max_logFC, avg_log2FC)) #%>%

vplot <- EnhancedVolcano(cluster_volcMarkers_file,
                         lab = cluster_volcMarkers_file$gene,
                         x = 'avg_log2FC',
                         y = 'p_val',
                         xlim = c(-1,1),
                         ylim = c(0,150),
                         title = '',
                         pCutoff = 0.05,
                         FCcutoff = 0.15,
                         pointSize = 2,
                         labSize = 0,
                         legendPosition = 'top',
                         legendLabSize = 18,
                         drawConnectors = F,
                         gridlines.major = FALSE,
                         gridlines.minor = FALSE,
                         #, selectLab = c('')
                         )

pdf(file = '', pointsize = 10, width = 8, height = 9)
  vplot
dev.off()


