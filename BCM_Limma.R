setwd("~/Box Sync/Memory_B_Pepper")
#Load packages needed for Limma linear modeling and DE analysis
library(limma)
library(xlsx)
library(RColorBrewer) 
library(gplots)
library(plyr)
library(ggplot2); theme_set(theme_bw(20) + theme(panel.grid.major = element_blank(), 
                                                 panel.grid.minor = element_blank()))
library(edgeR)
library(ggthemes)
library(WGCNA)
library(statmod)
options(stringsAsFactors = FALSE);
#Set seed for reproducability 
set.seed(123456)
#Load data from QC.R
lnames = load(file = "R_Analysis/counts_design_with_RAinfo.RData")

#Make sure the ordering of libs in counts matches that in design
libOrder <- match(design$Lib.ID, colnames(DGEcounts$counts))
DGEcounts <- DGEcounts[,libOrder]

#merge metrics with design
design <- merge(design, Metrics, by.x = "Lib.ID", by.y = "Lib.ID")

#Make sure design_metrics libs are in the same order as in counts
lib_order <- match(colnames(DGEcounts), design$Lib.ID)
design <- design[lib_order,]

#Prepare design matrix columns for analysis
design$Mouse <- as.factor(design$Mouse)
design$Ig <- as.factor(design$Ig)

#Subset the data for the "best of" analysis
#KeepDesign <- c("MS2_D", "MS3_D", "MS4_D", "MS3_M", "MS4_M", "MS5_M", "MS1_G", "MS2_G", "MS5_G")
#wholedesign <- design
#KeepRows <- which(design$Sample.ID %in% KeepDesign)
#design <- design[KeepRows,]
#libOrder <- match(design$Lib.ID, colnames(DGEcounts$counts))
#DGEcounts <- DGEcounts[,libOrder]

#Design & select models
design_mat1 <- model.matrix(~design$Ig)
design_mat2 <- model.matrix(~design$Ig + design$Mouse)
design_mat3 <- model.matrix(~design$Ig + design$PF_INDEL_RATE)
design_mat4 <- model.matrix(~design$Ig + design$MEDIAN_5PRIME_TO_3PRIME_BIAS)
design_mat5 <- model.matrix(~design$Ig + design$MEDIAN_CV_COVERAGE)
design_mat6 <- model.matrix(~design$Ig + design$MEDIAN_CV_COVERAGE + design$PF_INDEL_RATE)
design_mat7 <- model.matrix(~design$Ig + design$MEDIAN_CV_COVERAGE + design$Mouse)
designlist <- list(design_mat1, design_mat2, design_mat3, design_mat4, design_mat5, design_mat6, design_mat7)
vwts <- voomWithQualityWeights(DGEcounts, design=design_mat1, plot=TRUE, span=0.1)
out <- selectModel(vwts, designlist, criterion = "bic")
table(out$pref)
##BIC doesn't show anything beating the basic model.
##As far as quality metrics go, Median CV is the best to include in the model, and using that alone is better than using it with PF_Indel_Rate.
##Adding mouse into the model doesn't strengthen it, as per BIC.  (I might still use it.)

#Show DE genes with FC > 1 and p value <0.01
get_sig_genes <- function(top_genes){
  #Threshold 
  pCut <- 0.01
  foldCut <- 1
  sig_genes = top_genes[top_genes$adj.P.Val<pCut & abs(top_genes$logFC) > foldCut, ]
  return(sig_genes)
} # end get_genes function

#Ease up on significance calling
get_less_sig_genes <- function(top_genes){
  #Threshold 
  pCut <- 0.05
  foldCut <- 0.5849
  sig_genes = top_genes[top_genes$adj.P.Val<pCut & abs(top_genes$logFC) > foldCut, ]
  return(sig_genes)
} # end get_genes function

##Start with basic model
#Get model fits
design_mat1 <- model.matrix(~0+design$Ig)
colnames(design_mat1) <- c("IgD","IgG","IgM")
contrast.matrix <- makeContrasts(IgD-IgM, IgD-IgG, IgM-IgG, levels=design_mat1)
vwts1 <- voomWithQualityWeights(DGEcounts, design=design_mat1, plot=TRUE, span=0.1)

vfit1a <- lmFit(vwts1, design = design_mat1)
vfit1 <- contrasts.fit(vfit1a, contrast.matrix)
vfit_eb1 <- eBayes(vfit1)
design_batches <- model.matrix(~design$Ig)
counts_batch_corr1 <- vwts1$E
top_genes_DvM1 <-topTable (vfit_eb1, coef = 1, number=Inf, sort.by="P")
top_genes_DvG1 <-topTable (vfit_eb1, coef = 2, number=Inf, sort.by="P")
top_genes_MvG1 <-topTable (vfit_eb1, coef = 3, number=Inf, sort.by="P")
sig_genes_DvM1 <- get_sig_genes(top_genes_DvM1) #0 genes
sig_genes_DvG1 <- get_sig_genes(top_genes_DvG1) #0 genes
sig_genes_MvG1 <- get_sig_genes(top_genes_MvG1) #0 genes (3 genes when mouse 1 is excluded - Zdhhc2, Wbscr17, Kctd17)
less_sig_genes_DvM1 <- get_less_sig_genes(top_genes_DvM1) #8 genes (0 when mouse 1 is excluded)
less_sig_genes_DvG1 <- get_less_sig_genes(top_genes_DvG1) #10 genes (12 when mouse 1 is excluded)
less_sig_genes_MvG1 <- get_less_sig_genes(top_genes_MvG1) #3 genes (14 when mouse 1 is excluded)
results <- decideTests(vfit1)
vennDiagram(results)

#Write DE gene lists
write.csv(sig_genes_DvM1,"Limma/Basic/sig_genes_DvM1.csv",
          quote=FALSE, row.names = TRUE)
write.csv(less_sig_genes_DvM1, "Limma/Basic/less_sig_genes_DvM1.csv",
          quote = F, row.names = T)
write.csv(top_genes_DvM1, "Limma/Basic/top_genes_DvM1.csv",
          quote = F, row.names = T)

write.csv(sig_genes_DvG1,"Limma/Basic/sig_genes_DvG1.csv",
          quote=FALSE, row.names = TRUE)
write.csv(less_sig_genes_DvG1, "Limma/Basic/less_sig_genes_DvG1.csv",
          quote = F, row.names = T)
write.csv(top_genes_DvG1, "Limma/Basic/top_genes_DvG1.csv",
          quote = F, row.names = T)

write.csv(sig_genes_MvG1,"Limma/Basic/sig_genes_MvG1.csv",
          quote=FALSE, row.names = TRUE)
write.csv(less_sig_genes_MvG1, "Limma/Basic/less_sig_genes_MvG1.csv",
          quote = F, row.names = T)
write.csv(top_genes_MvG1, "Limma/Basic/top_genes_MvG1.csv",
          quote = F, row.names = T)

###
##Add quality metrics to the model <- This is the model I'm using for Marion and Akshay
design_mat5 <- model.matrix(~0+design$Ig+design$MEDIAN_CV_COVERAGE)
colnames(design_mat5) <- c("IgD","IgG","IgM","MedCV")
contrast.matrix5 <- makeContrasts(IgD-IgM, IgD-IgG, IgM-IgG, levels=design_mat5)
vwts5 <- voomWithQualityWeights(DGEcounts, design=design_mat5, plot=TRUE, span=0.1)
#Get model fits
vfit5a <- lmFit(vwts5, design = design_mat5)
vfit5 <- contrasts.fit(vfit5a, contrast.matrix5)
vfit_eb5 <- eBayes(vfit5)
volcanoplot(vfit_eb5, coef = 1)
volcanoplot(vfit_eb5, coef = 2)
volcanoplot(vfit_eb5, coef = 3)

design_batches <- model.matrix(~0+design$Ig)
counts_batch_corr5 <- removeBatchEffect(vwts5$E, covariates = design$MEDIAN_CV_COVERAGE,
                      design = design_batches)
top_genes_DvM5 <-topTable (vfit_eb5, coef = 1, number=Inf, sort.by="P")
top_genes_DvG5 <-topTable (vfit_eb5, coef = 2, number=Inf, sort.by="P")
top_genes_MvG5 <-topTable (vfit_eb5, coef = 3, number=Inf, sort.by="P")
sig_genes_DvM5 <- get_sig_genes(top_genes_DvM5) #0 genes
sig_genes_DvG5 <- get_sig_genes(top_genes_DvG5) #2 genes
sig_genes_MvG5 <- get_sig_genes(top_genes_MvG5) #0 genes
less_sig_genes_DvM5 <- get_less_sig_genes(top_genes_DvM5) #84 genes (0 if you exclude mouse 1)
less_sig_genes_DvG5 <- get_less_sig_genes(top_genes_DvG5) #114 genes (5 if you exclude mouse 1 - Zdhhc2, Tmem177, Ctnnal1, Kctd17, 3300002I08Rik)
less_sig_genes_MvG5 <- get_less_sig_genes(top_genes_MvG5) #2 genes (only 1 if you exclude mouse 1 - Zdhhc2)
results <- decideTests(vfit5)
vennDiagram(results)

plotMD(vfit_eb5, coef=1, status=results[,1], main=colnames(vfit_eb5)[1], xlim=c(-2,11), hl.col = c("royalblue3","indianred"), bg.col="grey50")
abline(0,0,col="grey")
plotMD(vfit_eb5, coef=2, status=results[,2], main=colnames(vfit_eb5)[2], xlim=c(-2,11), hl.col = c("palegreen4","indianred"), bg.col="grey50")
abline(0,0,col="grey")
plotMD(vfit_eb5, coef=3, status=results[,3], main=colnames(vfit_eb5)[3], xlim=c(-2,11), hl.col = c("royalblue3","palegreen4"), bg.col="grey50")
abline(0,0,col="grey")


#Write DE gene lists
write.csv(sig_genes_DvM5,"Limma/With MedCV/sig_genes_DvM5.csv",
          quote=FALSE, row.names = TRUE)
write.csv(less_sig_genes_DvM5, "Limma/With MedCV/less_sig_genes_DvM5.csv",
          quote = F, row.names = T)
write.csv(top_genes_DvM5, "Limma/With MedCV/top_genes_DvM5.csv",
          quote = F, row.names = T)

write.csv(sig_genes_DvG5,"Limma/With MedCV/sig_genes_DvG5.csv",
          quote=FALSE, row.names = TRUE)
write.csv(less_sig_genes_DvG5, "Limma/With MedCV/less_sig_genes_DvG5.csv",
          quote = F, row.names = T)
write.csv(top_genes_DvG5, "Limma/With MedCV/top_genes_DvG5.csv",
          quote = F, row.names = T)

write.csv(sig_genes_MvG5,"Limma/With MedCV/sig_genes_MvG5.csv",
          quote=FALSE, row.names = TRUE)
write.csv(less_sig_genes_MvG5, "Limma/With MedCV/less_sig_genes_MvG5.csv",
          quote = F, row.names = T)
write.csv(top_genes_MvG5, "Limma/With MedCV/top_genes_MvG5.csv",
          quote = F, row.names = T)

###
#Add mouse into the model (since the mice faintly cluster on PCA) -> loses statistical power
design_mat7 <- model.matrix(~0+design$Mouse+design$Ig+design$MEDIAN_CV_COVERAGE)
colnames(design_mat7) <- c("Mouse1","Mouse2","Mouse3","Mouse4","Mouse5","IgG","IgM","MedCV")
contrast.matrix7 <- makeContrasts(IgM, IgG, IgM-IgG, levels=design_mat7)
vwts7 <- voomWithQualityWeights(DGEcounts, design=design_mat7, plot=TRUE, span=0.1)

vfit7a <- lmFit(vwts7, design = design_mat7)
vfit7 <- contrasts.fit(vfit7a, contrast.matrix7)
vfit_eb7 <- eBayes(vfit7)
design_batches <- model.matrix(~design$Ig)
counts_batch_corr7 <- removeBatchEffect(vwts7$E, batch = design$Mouse,
                                        covariates = design$MEDIAN_CV_COVERAGE,
                                        design = design_batches)
top_genes_DvM7 <-topTable (vfit_eb7, coef = 1, number=Inf, sort.by="P")
top_genes_DvG7 <-topTable (vfit_eb7, coef = 2, number=Inf, sort.by="P")
top_genes_MvG7 <-topTable (vfit_eb7, coef = 3, number=Inf, sort.by="P")

sig_genes_DvM7 <- get_sig_genes(top_genes_DvM7) #0 genes
sig_genes_DvG7 <- get_sig_genes(top_genes_DvG7) #0 genes
sig_genes_MvG7 <- get_sig_genes(top_genes_MvG7) #0 genes
less_sig_genes_DvM7 <- get_less_sig_genes(top_genes_DvM7) #0 genes
less_sig_genes_DvG7 <- get_less_sig_genes(top_genes_DvG7) #0 genes
less_sig_genes_MvG7 <- get_less_sig_genes(top_genes_MvG7) #0 genes
###
design_mat7b <- model.matrix(~0+design$Ig+design$MEDIAN_CV_COVERAGE)
colnames(design_mat7b) <- c("IgD","IgG","IgM","MedCV")
contrast.matrix7b <- makeContrasts(IgD-IgM, IgD-IgG, IgM-IgG, levels=design_mat7b)

vwts7b <-  voomWithQualityWeights(DGEcounts, design=design_mat7b)
corfit <- duplicateCorrelation(vwts7b, design_mat7b, block=design$Mouse)
vwts7b <-  voomWithQualityWeights(DGEcounts, design=design_mat7b, block=design$Mouse, correlation=corfit$consensus)
fit7 <- lmFit(vwts7b, design_mat7b, block=design$Mouse, correlation=corfit$consensus)
vfit7b <- contrasts.fit(fit7, contrast.matrix7b)
vfit_eb7b <- eBayes(vfit7b)
design_batches <- model.matrix(~design$Ig)
counts_batch_corr7b <- removeBatchEffect(vwts7b$E, batch = design$Mouse,
                                        covariates = design$MEDIAN_CV_COVERAGE,
                                        design = design_batches)
top_genes_DvM7 <-topTable (vfit_eb7b, coef = 1, number=Inf, sort.by="P")
top_genes_DvG7 <-topTable (vfit_eb7b, coef = 2, number=Inf, sort.by="P")
top_genes_MvG7 <-topTable (vfit_eb7b, coef = 3, number=Inf, sort.by="P")

sig_genes_DvM7 <- get_sig_genes(top_genes_DvM7) #0 genes
sig_genes_DvG7 <- get_sig_genes(top_genes_DvG7) #0 genes
sig_genes_MvG7 <- get_sig_genes(top_genes_MvG7) #0 genes
less_sig_genes_DvM7 <- get_less_sig_genes(top_genes_DvM7) #0 genes
less_sig_genes_DvG7 <- get_less_sig_genes(top_genes_DvG7) #4 genes - Ctnnal1, Tmem177, Zdhhc2, Gm15991
less_sig_genes_MvG7 <- get_less_sig_genes(top_genes_MvG7) #1 genes - Ighg3

