'''
BINF6110 Assignment2 
By: Jacob Hambly

'''

##Loading in the necessary libraries 
library(Rsubread)
library(edgeR)
library(ggplot2)
library(ComplexHeatmap)
library(ggpubr)

##Assigning the path to the yeast annotated genome file to the variable 'gtf.file'
gtf.file <- file.path("genomic.gtf")

##Assigning the names of the nine .bam files in the current directory to the variable name 'bam.files'
bam.files <- list.files(path = ".", pattern = "*.bam", full.names = T)

##Assigning the aligned .bam files to a specific genomic features using the yeast reference genome
feature_counts <- featureCounts(files = bam.files, annot.ext = gtf.file, isGTFAnnotationFile = T, isPairedEnd = F)

##Looking at feature_counts

#getting the summary of the feature counts
summary(feature_counts)

#looking at the first few rows of feature_counts$annotation
head(feature_counts$annotation)

#number of matches for each feature in the yeast genome
rawdata <- feature_counts$counts
rawdata

##Using feature_counts to create a DGEList object.
rnaseq_reads <- DGEList(counts = feature_counts$counts, genes = feature_counts$annotation)
rnaseq_reads 

##Assigning group labels for each of the samples based stage for the sample. E: Early biofilm, T: Thin biofilm, M: Mature biofilm
rnaseq_reads$samples$group <- c('M', 'M', 'M', 'T', 'T', 'T', 'E', 'E', 'E')
rnaseq_reads$samples
rnaseq_reads$counts

##Filtering out genes that have low expression across each sample. 
keep <- filterByExpr(rnaseq_reads)

#looking at a table for the number of genes kept or dropped from rna_reads
table(keep)

#keeping genes that are TRUE (there is enough gene expression across the samples).
rnaseq_reads <- rnaseq_reads[keep,]


## Adjusting the read counts in each sample using the normLibSizes function to account for differences in library size across samples.
rnaseq_reads <- normLibSizes(rnaseq_reads)

#Looking at the norm.factor column
rnaseq_reads$samples

##Plotting the samples from rnaseq_reads on a two-dimensional scatter plot so that distances on the plot approximate the typical log2 fold changes between the samples.
plotMDS(rnaseq_reads, labels =rnaseq_reads$samples$group)

## Creating a design matrix

#creating factors that represent the group for each sample
group <- factor(c('Mature_BioFilm','Mature_BioFilm', 'Mature_BioFilm', 'Thin_BioFilm', 'Thin_BioFilm', 'Thin_BioFilm','Early_BioFilm', 'Early_BioFilm', 'Early_BioFilm'))


#creating a model (design) matrix called design where where each column represents a specific group, and each row represents a sample. Intercept is left out so each group is a feature.
design <- model.matrix(~0+group)
rownames(design) <- colnames(rnaseq_reads)
design

##Estimating the dispersion parameter for each gene across the samples.
rnaseq_reads <- estimateDisp(rnaseq_reads, design)

##Fitting a negative binomial generalized log-linear model for each gene to estimate total abundance using the features from the design matrix as the features.
fitQ <- glmQLFit(rnaseq_reads, design)

#looking at the linear models created for each gene
coef(fitQ)

## Conducting a quasi-likelihood f test to look for genes that are differentially expressed between two of the groups. Here we are comparing early to mature, 1 early 0 thin -1 mature
Qlrt_early_mature <- glmQLFTest(fitQ, contrast = c(1,-1,0)) 


## Conducting a quasi-likelihood f test to look for genes that are differentially expressed between two of the groups. Here we are comparing early to thin, 1 early -1 thin 0 mature

Qlrt_early_thin <- glmQLFTest(fitQ, contrast = c(1,0,-1))

## Conducting a lquasi-likelihood f test to look for genes that are differentially expressed between two of the groups. Here we are comparing thin to mature,early 1 thin -1 mature
Qlrt_thin_mature <- glmQLFTest(fitQ, contrast = c(0,-1,1))


##Identifying the number of genes that are significantly downregulated or upregulated in the early group compared to the mature group as well as the number of genes that don't show significant differences in expression levels between the two groups.
summary(decideTests(Qlrt_early_mature))

##Identifying the number of genes that are significantly downregulated or upregulated in the early group compared to the thin group as well as the number of genes that don't show significant differences in expression levels between the two groups.
summary(decideTests(Qlrt_early_thin))

##Identifying the number of genes that are significantly downregulated or upregulated in the thin group compared to the mature group as well as the number of genes that don't show significant differences in expression levels between the two groups.
summary(decideTests(Qlrt_thin_mature))

#extracting only the significantly differentially expressed genes (adjusted p-value < 0.05) for early vs mature and saving the genes to a new data frame called de_early_mature
de_early_mature <- topTags(object = Qlrt_early_mature, n = Inf, p.value = 0.05, adjust.method = 'BH')[,c(1,7,10)]
de_early_mature <- as.data.frame(de_early_mature)
de_early_mature["Treatments"] <- "Early vs Mature"
dim(de_early_mature)


#extracting only the significant differentially expressed genes (adjusted p-value < 0.05) for early vs thin and saving the genes to a new data frame called de_early_thin
de_early_thin <- topTags(object = Qlrt_early_thin, n = Inf, p.value = 0.05, adjust.method = 'BH')[,c(1,7,10)]
de_early_thin <- as.data.frame(de_early_thin)
de_early_thin["Treatments"] <- "Early vs Thin" 
dim(de_early_thin)

#extracting only the significantly differentially expressed genes (adjusted p-value < 0.05) for thin vs mature and saving the genes to a new data frame called de_thin_mature
de_thin_mature <- topTags(object = Qlrt_thin_mature, n = Inf, p.value = 0.05, adjust.method = 'BH' )[,c(1,7,10)]
de_thin_mature <- as.data.frame(de_thin_mature)
de_thin_mature["Treatments"] <- "Thin vs Mature"
dim(de_thin_mature)

#combining de_early_mature, de_early_thin, de_thin_mature into one data frame called de_any_two_treatments
de_any_two_treatments <- rbind(de_early_mature,de_early_thin, de_thin_mature)

#saving de_any_two_treatments as a csv file called "gene_diffexp_any.csv"
write.csv(de_any_two_treatments, file = "gene_diffexp_any.csv", row.names = FALSE)


##Comparing early with both thin and mature

#Creating a new design matrix where the features are "Early" for the early group and "Not_Early" for the thin and mature group
group2 <- factor(c("Not_Early", "Not_Early", "Not_Early", "Not_Early", "Not_Early", "Not_Early", "Early", "Early", "Early" ))
design2 <- model.matrix(~0+group2)
rownames(design2) <- colnames(rnaseq_reads)
design2

#creating new groups
rnaseq_reads$samples$group <- c("NE", "NE", "NE", "NE", "NE", "NE", "E", "E", "E")

##Estimating the dispersion parameter for each gene across the samples.
rnaseq_reads2 <- estimateDisp(rnaseq_reads, design2)

##Fitting a negative binomial generalized log-linear model for each gene to estimate total abundance using the features from the design matrix as the features. 
fitQ2 <- glmQLFit(rnaseq_reads2, design2)

##Looking at the linear models created for each gene
coef(fitQ2)

## Conducting a quasi-likelihood f-test to look for genes that are differentially expressed between the two groups. Here we are comparing early to not_early, 1 early -1 not_early.
Qlrt_early_other <- glmQLFTest(fitQ2, contrast = c(1,-1)) 

##Identifying the number of genes that are significantly downregulated or upregulated in the early group compared to the not_early group as well as the number of genes that don't show significant differences in expression levels between the two groups.
summary(decideTests(Qlrt_early_other))

#extracting only the significantly differentially expressed genes (adjusted p-value < 0.05) for early vs not_early and saving the genes to a new data frame called de_early_other
de_early_other <- topTags(Qlrt_early_other, n = Inf, p.value = 0.05, adjust.method = 'BH')[,c(1,7,10)]
de_early_other <- as.data.frame(de_early_other)
dim(de_early_other)


#saving de_early_other as a csv file called gene_diffexp_earlyvsother.csv.
write.csv(de_early_other, file = "gene_diffexp_earlyvsother.csv")



###making volcano plot
##Early vs Mature
#creating a data frame with all the gene expression information between early and mature and creating a new feature with the -log10( adjusted p-value)
de_early_mature_all <- topTags(object = Qlrt_early_mature, n = Inf, adjust.method = "BH")[,c(1,7,10)]
de_early_mature_all <- as.data.frame(de_early_mature_all)
de_early_mature_all['-Log10(p-value)'] <- -log10(de_early_mature_all[,3])
de_early_mature_all['differential_expression'] <- NA
de_early_mature_all[de_early_mature_all$logFC >= 1, 5] <- 'UP'
de_early_mature_all[de_early_mature_all$logFC <= -1, 5] <- 'DOWN'
de_early_mature_all[de_early_mature_all$logFC > -1 & de_early_mature_all$logFC < 1 , 5] <- 'NO'

#looking at the number of upregulated and downregulated genes based on fold change
table(de_early_mature_all$differential_expression)

#creating volcano plot
vp_early_mature <- ggplot(data = de_early_mature_all, aes(x = de_early_mature_all$logFC, y =de_early_mature_all$`-Log10(p-value)`, col = de_early_mature_all$differential_expression)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), col = 'grey', linetype = 'dashed') +
  geom_vline(xintercept = c(-1,1), col = 'grey', linetype = 'dashed') +
  scale_color_manual(name = "Gene Expression", values = c("#1B9E77", "#666666", "#D95F02"), labels = c("Downregulated", "No signifance", "Upregulated")) +
  ggtitle("Early vs Mature") +
  labs(x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value"))
##################################################################################
##Early vs Thin
#creating a data frame with all the gene expression information between early and thin and creating a new feature with the -log10(adjusted p-value)
de_early_thin_all <- topTags(object = Qlrt_early_thin, n = Inf, adjust.method = "BH")[,c(1,7,10)]
de_early_thin_all <- as.data.frame(de_early_thin_all)
de_early_thin_all['-Log10(p-value)'] <- -log10(de_early_thin_all[,3])
de_early_thin_all['differential_expression'] <- NA
de_early_thin_all[de_early_thin_all$logFC >= 1, 5] <- 'UP'
de_early_thin_all[de_early_thin_all$logFC <= -1, 5] <- 'DOWN'
de_early_thin_all[de_early_thin_all$logFC > -1 & de_early_thin_all$logFC < 1 , 5] <- 'NO'

#looking at the number of genes that are upregulated or downregulated based of fold change
table(de_early_thin_all$differential_expression)

#creating volcano plot
vp_early_thin <- ggplot(data = de_early_thin_all, aes(x = de_early_thin_all$logFC, y =de_early_thin_all$`-Log10(p-value)`, col = de_early_thin_all$differential_expression)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), col = 'grey', linetype = 'dashed') +
  geom_vline(xintercept = c(-1,1), col = 'grey', linetype = 'dashed') +
  scale_color_manual(name = "Gene Expression", values = c("#1B9E77", "#666666", "#D95F02"), labels = c("Downregulated", "No signifance", "Upregulated")) +
  ggtitle("Early vs Thin") +
  labs(x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value"))
##################################################################################
##Thin vs Mature
#creating a data frame with all the gene expression information between Thin and mature and creating a new feature with the -log10(adjusted p-value)
de_thin_mature_all <- topTags(object = Qlrt_thin_mature, n = Inf, adjust.method = "BH")[,c(1,7,10)]
de_thin_mature_all <- as.data.frame(de_thin_mature_all)
de_thin_mature_all['-Log10(p-value)'] <- -log10(de_thin_mature_all[,3])
de_thin_mature_all['differential_expression'] <- NA
de_thin_mature_all[de_thin_mature_all$logFC >= 1, 5] <- 'UP'
de_thin_mature_all[de_thin_mature_all$logFC <= -1, 5] <- 'DOWN'
de_thin_mature_all[de_thin_mature_all$logFC > -1 & de_thin_mature_all$logFC < 1 , 5] <- 'NO'

table(de_thin_mature_all$differential_expression)

#creating volcano plot
vp_thin_mature <- ggplot(data = de_thin_mature_all, aes(x = de_thin_mature_all$logFC, y =de_thin_mature_all$`-Log10(p-value)`, col = de_thin_mature_all$differential_expression)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), col = 'grey', linetype = 'dashed') +
  geom_vline(xintercept = c(-1,1), col = 'grey', linetype = 'dashed') +
  scale_color_manual(name = "Gene Expression", values = c("#1B9E77", "#666666", "#D95F02"), labels = c("Downregulated", "No signifance", "Upregulated")) +
  ggtitle("Thin vs Mature") +
  labs(x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value"))
##################################################################################
##Early vs others
#creating a data frame with all the gene expression information between early and other and creating a new feature with the -log10(adjusted p-value)
de_early_other_all <- topTags(object = Qlrt_early_other, n = Inf, adjust.method = "BH")[,c(1,7,10)]
de_early_other_all <- as.data.frame(de_early_other_all)
de_early_other_all['-Log10(p-value)'] <- -log10(de_early_other_all[,3])
de_early_other_all['differential_expression'] <- NA
de_early_other_all[de_early_other_all$logFC >= 1, 5] <- 'UP'
de_early_other_all[de_early_other_all$logFC <= -1, 5] <- 'DOWN'
de_early_other_all[de_early_other_all$logFC > -1 & de_early_other_all$logFC < 1 , 5] <- 'NO'

#looking at the number of genes identified as upregulated or downregulated based on fold change
table(de_early_other_all$differential_expression)

#creating volcano plot
vp_early_other <- ggplot(data = de_early_other_all, aes(x = de_early_other_all$logFC, y =de_early_other_all$`-Log10(p-value)`, col = de_early_other_all$differential_expression)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), col = 'grey', linetype = 'dashed') +
  geom_vline(xintercept = c(-1,1), col = 'grey', linetype = 'dashed') +
  scale_color_manual(name = "Gene Expression", values = c("#1B9E77", "#666666", "#D95F02"), labels = c("Downregulated", "No signifance", "Upregulated")) +
  ggtitle("Early vs Later (Thin and Mature)") +
  labs(x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value"))
##################################################################################

##Arranging the 4 volcano plots into 2, 1x2, format
ggarrange(vp_early_mature, vp_early_thin, nrow = 1, ncol = 2, common.legend = TRUE)
ggarrange(vp_thin_mature, vp_early_other, nrow = 1, ncol = 2, legend = 'none')

###PCA 

#getting the normalzied counts data
normalized_counts <- cpm(rnaseq_reads, normalized.lib.sizes = TRUE)
#preparing the counts data frame
normalized_counts_pca <- as.data.frame(t(normalized_counts))
#using prcomp to calculate the principle components
sample_PCA <- prcomp(normalized_counts_pca)
pc_scores <- as.data.frame(sample_PCA$x)
pc_scores["Stage"] <- c("Mature", "Mature", "Mature", "Thin", "Thin", "Thin", "Early", "Early", "Early")

#extracting the proportion of variance for pc1 and pc2
proportion_of_variance_pc1 <- summary(sample_PCA)$importance[2,1] * 100
proportion_of_variance_pc2 <- summary(sample_PCA)$importance[2,2] * 100

ggplot(data = pc_scores, aes(x = PC1, y = PC2 , col = Stage)) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = paste("PC1 (", proportion_of_variance_pc1, "% variance)"), 
       y = paste("PC2 (", proportion_of_variance_pc2, "% variance)"))

###Heatmap

#getting the z-score for the normalized counts
counts.z <- scale(normalized_counts)

#extracting the top 100 differentially expressed genes between early and not early (thin and mature)
de_early_other_100 <- topTags(Qlrt_early_other, n = 100, sort.by = 'PValue', adjust.method = 'BH')

#creating a data frame representing the count data for the top 100 differentially expressed genes between early and not early (thin and mature)
noralized_counts_early_other <- counts.z[row.names(de_early_other_100),]

#plotting the heatmap
Heatmap(noralized_counts_early_other, column_labels = c( "Not Early(Mature)", "Not Early(Mature)", "Not Early(Mature)", "Not Early(Thin)", "Not Early(Thin)", "Not Early(Thin)", "Early", "Early", "Early"), row_labels = rep(" ", 100), heatmap_legend_param = list(title = "z-score of the normalized counts"))
