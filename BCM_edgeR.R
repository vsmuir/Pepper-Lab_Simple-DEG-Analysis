setwd("~/Box Sync/Memory_B_Pepper")

library(edgeR)

options(stringsAsFactors = FALSE);
#Set seed for reproducability 
set.seed(123456)

#Load data from QC.R
load(file = "R_Analysis/counts_design_with_BCMinfo.RData")

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

#Load function to show DE genes with FC > 1 and p value <0.01
get_sig_genes <- function(top_genes){
  #Threshold 
  pCut <- 0.01
  sig_genes = top_genes[top_genes$FDR<pCut, ]
  return(sig_genes)
} # end get_genes function

#Load function that eases up on significance calling
get_less_sig_genes <- function(top_genes){
  #Threshold 
  pCut <- 0.05
  sig_genes = top_genes[top_genes$FDR<pCut, ]
  return(sig_genes)
} # end get_genes function

###Classic edgeR analysis (Basic)
basicBCM <- DGEcounts
basicBCM$samples$group <- design$Ig
basicBCM <- estimateDisp(basicBCM)
et_basic_DvG <- exactTest(basicBCM, pair=c("G","D"))
et_basic_DvM <- exactTest(basicBCM, pair=c("M","D"))
et_basic_MvG <- exactTest(basicBCM, pair=c("G","M"))

top_genes_DvG <- topTags(et_basic_DvG, n = Inf, sort.by="P")
top_genes_DvG <- top_genes_DvG$table
top_genes_DvM <- topTags(et_basic_DvM, n = Inf, sort.by="P")
top_genes_DvM <- top_genes_DvM$table
top_genes_MvG <- topTags(et_basic_MvG, n = Inf, sort.by="P")
top_genes_MvG <- top_genes_MvG$table

sig_genes_DvM <- get_sig_genes(top_genes_DvM) #3 genes
sig_genes_DvG <- get_sig_genes(top_genes_DvG) #7 genes
sig_genes_MvG <- get_sig_genes(top_genes_MvG) #0 genes 
less_sig_genes_DvM <- get_less_sig_genes(top_genes_DvM) #28 genes 
less_sig_genes_DvG <- get_less_sig_genes(top_genes_DvG) #28 genes 
less_sig_genes_MvG <- get_less_sig_genes(top_genes_MvG) #0 genes 

write.csv(sig_genes_DvM,"edgeR/Basic/sig_genes_DvM.csv",
          quote=FALSE, row.names = TRUE)
write.csv(less_sig_genes_DvM, "edgeR/Basic/less_sig_genes_DvM.csv",
          quote = F, row.names = T)
write.csv(top_genes_DvM, "edgeR/Basic/top_genes_DvM.csv",
          quote = F, row.names = T)

write.csv(sig_genes_DvG,"edgeR/Basic/sig_genes_DvG.csv",
          quote=FALSE, row.names = TRUE)
write.csv(less_sig_genes_DvG, "edgeR/Basic/less_sig_genes_DvG.csv",
          quote = F, row.names = T)
write.csv(top_genes_DvG, "edgeR/Basic/top_genes_DvG.csv",
          quote = F, row.names = T)

write.csv(sig_genes_MvG,"edgeR/Basic/sig_genes_MvG.csv",
          quote=FALSE, row.names = TRUE)
write.csv(less_sig_genes_MvG, "edgeR/Basic/less_sig_genes_MvG.csv",
          quote = F, row.names = T)
write.csv(top_genes_MvG, "edgeR/Basic/top_genes_MvG.csv",
          quote = F, row.names = T)

###GLM edgeR analysis (with median CV)
#Build Model Matrix for analysis
cvBCM <- DGEcounts
design_mat5 <- model.matrix(~0+design$Ig+design$MEDIAN_CV_COVERAGE)
colnames(design_mat5) <- c("IgD","IgG","IgM","MedCV")
#Find dispersions and test with model
cvBCM <- estimateDisp(cvBCM, design_mat5)
anova.con <- makeContrasts(DvG=IgD-IgG, DvM=IgD-IgM, MvG=IgM-IgG, levels=design_mat5)
fit <- glmFit(cvBCM, design_mat5)
lrt.anova <- glmLRT(fit, contrast=anova.con)
lrt.DvM <- glmLRT(fit, contrast = anova.con[,"DvM"])
lrt.DvG <- glmLRT(fit, contrast = anova.con[,"DvG"])
lrt.MvG <- glmLRT(fit, contrast = anova.con[,"MvG"])
#Output list of genes with fold changes and p-values
top_genes_DvG1.lrt <- topTags(lrt.DvG, n = Inf, sort.by="P")
top_genes_DvG1.lrt <- top_genes_DvG1.lrt$table
top_genes_DvM1.lrt <- topTags(lrt.DvM, n = Inf, sort.by="P")
top_genes_DvM1.lrt <- top_genes_DvM1.lrt$table
top_genes_MvG1.lrt <- topTags(lrt.MvG, n = Inf, sort.by="P")
top_genes_MvG1.lrt <- top_genes_MvG1.lrt$table
top_genes_anova.lrt <- topTags(lrt.anova, n = Inf, sort.by="P")
top_genes_anova.lrt <- top_genes_anova.lrt$table
#Select significantly DE genes
sig_genes_DvM1.lrt <- get_sig_genes(top_genes_DvM1.lrt) #27 genes
sig_genes_DvG1.lrt <- get_sig_genes(top_genes_DvG1.lrt) #16 genes
sig_genes_MvG1.lrt <- get_sig_genes(top_genes_MvG1.lrt) #0 genes 
sig_genes_anova.lrt <- get_sig_genes(top_genes_anova.lrt) #24 genes 
less_sig_genes_DvM1.lrt <- get_less_sig_genes(top_genes_DvM1.lrt) #110 genes 
less_sig_genes_DvG1.lrt <- get_less_sig_genes(top_genes_DvG1.lrt) #104 genes 
less_sig_genes_MvG1.lrt <- get_less_sig_genes(top_genes_MvG1.lrt) #0 genes 
less_sig_genes_anova.lrt <- get_less_sig_genes(top_genes_anova.lrt) #83 genes 
## The QL F-test is preferred as it reflects the uncertainty in estimating the dispersion for each gene. 
## It provides more robust and reliable error rate control when the number of replicates is small.
## It also was unable to find any DE genes when I applied it to this dataset.

write.csv(sig_genes_DvM1.lrt,"edgeR/With MedCV/sig_genes_DvM.csv",
          quote=FALSE, row.names = TRUE)
write.csv(less_sig_genes_DvM1.lrt, "edgeR/With MedCV/less_sig_genes_DvM.csv",
          quote = F, row.names = T)
write.csv(top_genes_DvM1.lrt, "edgeR/With MedCV/top_genes_DvM.csv",
          quote = F, row.names = T)

write.csv(sig_genes_DvG1.lrt,"edgeR/With MedCV/sig_genes_DvG.csv",
          quote=FALSE, row.names = TRUE)
write.csv(less_sig_genes_DvG1.lrt, "edgeR/With MedCV/less_sig_genes_DvG.csv",
          quote = F, row.names = T)
write.csv(top_genes_DvG1.lrt, "edgeR/With MedCV/top_genes_DvG.csv",
          quote = F, row.names = T)

write.csv(sig_genes_MvG1.lrt,"edgeR/With MedCV/sig_genes_MvG.csv",
          quote=FALSE, row.names = TRUE)
write.csv(less_sig_genes_MvG1.lrt, "edgeR/With MedCV/less_sig_genes_MvG.csv",
          quote = F, row.names = T)
write.csv(top_genes_MvG1.lrt, "edgeR/With MedCV/top_genes_MvG.csv",
          quote = F, row.names = T)

write.csv(less_sig_genes_anova.lrt,"edgeR/With MedCV/less_sig_genes_anova.csv",
          quote = F, row.names = T)

###GLM edgeR analysis (with mouse and median CV) <- This model was used for Marion and Akshay
cvmBCM <- DGEcounts
#Build Model Matrix for analysis
design_mat7 <- model.matrix(~0+design$Mouse+design$MEDIAN_CV_COVERAGE+design$Ig)
colnames(design_mat7) <- c("Mouse1","Mouse2","Mouse3","Mouse4","Mouse5","MedCV","IgG","IgM")
#For the heatmap
vwts7 <- voomWithQualityWeights(DGEcounts, design=design_mat7, plot=TRUE, span=0.1)
design_batches <- model.matrix(~design$Ig)
counts_batch_corr7 <- removeBatchEffect(vwts7$E, batch = design$Mouse,
                                        covariates = design$MEDIAN_CV_COVERAGE,
                                        design = design_batches)


cvmBCM <- estimateDisp(cvmBCM, design_mat7)
anova.con.m <- makeContrasts(DvG=c(0,0,0,0,0,0,1,0), DvM=c(0,0,0,0,0,0,0,1), MvG=IgM-IgG, levels=design_mat7)
fit2 <- glmFit(cvmBCM, design_mat7)
lrt.anova2 <- glmLRT(fit2, contrast=anova.con.m)
lrt.DvM2 <- glmLRT(fit2, contrast = anova.con.m[,"DvM"])
lrt.DvG2 <- glmLRT(fit2, contrast = anova.con.m[,"DvG"])
lrt.MvG2 <- glmLRT(fit2, contrast = anova.con.m[,"MvG"])

top_genes_DvG2.lrt <- topTags(lrt.DvG2, n = Inf, sort.by="P")
top_genes_DvG2.lrt <- top_genes_DvG2.lrt$table
top_genes_DvM2.lrt <- topTags(lrt.DvM2, n = Inf, sort.by="P")
top_genes_DvM2.lrt <- top_genes_DvM2.lrt$table
top_genes_MvG2.lrt <- topTags(lrt.MvG2, n = Inf, sort.by="P")
top_genes_MvG2.lrt <- top_genes_MvG2.lrt$table
top_genes_anova2.lrt <- topTags(lrt.anova2, n = Inf, sort.by="P")
top_genes_anova2.lrt <- top_genes_anova2.lrt$table

sig_genes_DvM2.lrt <- get_sig_genes(top_genes_DvM2.lrt) #37 genes
sig_genes_DvG2.lrt <- get_sig_genes(top_genes_DvG2.lrt) #15 genes
sig_genes_MvG2.lrt <- get_sig_genes(top_genes_MvG2.lrt) #2 genes 
sig_genes_anova2.lrt <- get_sig_genes(top_genes_anova2.lrt) #33 genes 

less_sig_genes_DvM2.lrt <- get_less_sig_genes(top_genes_DvM2.lrt) #115 genes 
less_sig_genes_DvG2.lrt <- get_less_sig_genes(top_genes_DvG2.lrt) #92 genes 
less_sig_genes_MvG2.lrt <- get_less_sig_genes(top_genes_MvG2.lrt) #4 genes 
less_sig_genes_anova2.lrt <- get_less_sig_genes(top_genes_anova2.lrt) #130 genes 

DvM2up <- rownames(cvmBCM) %in% rownames(less_sig_genes_DvM2.lrt[(less_sig_genes_DvM2.lrt$logFC > 0),])
DvM2dn <- rownames(cvmBCM) %in% rownames(less_sig_genes_DvM2.lrt[(less_sig_genes_DvM2.lrt$logFC < 0),])
DvG2up <- rownames(cvmBCM) %in% rownames(less_sig_genes_DvG2.lrt[(less_sig_genes_DvG2.lrt$logFC > 0),])
DvG2dn <- rownames(cvmBCM) %in% rownames(less_sig_genes_DvG2.lrt[(less_sig_genes_DvG2.lrt$logFC < 0),])
MvG2up <- rownames(cvmBCM) %in% rownames(less_sig_genes_MvG2.lrt[(less_sig_genes_MvG2.lrt$logFC > 0),])
MvG2dn <- rownames(cvmBCM) %in% rownames(less_sig_genes_MvG2.lrt[(less_sig_genes_MvG2.lrt$logFC < 0),])

with(lrt.DvM2$table, plot(logCPM,logFC,pch=16,cex=0.2, col="grey50", main = "IgD - IgM"))
with(lrt.DvM2$table, (points(logCPM[DvM2up],logFC[DvM2up],pch=16,col="royalblue3")))
with(lrt.DvM2$table, (points(logCPM[DvM2dn],logFC[DvM2dn],pch=16,col="indianred"))) 
abline(0,0,col="grey")

with(lrt.DvG2$table, plot(logCPM,logFC,pch=16,cex=0.2, col="grey50", main = "IgD - IgG"))
with(lrt.DvG2$table, (points(logCPM[DvG2up],logFC[DvG2up],pch=16,col="palegreen4")))
with(lrt.DvG2$table, (points(logCPM[DvG2dn],logFC[DvG2dn],pch=16,col="indianred"))) 
abline(0,0,col="grey")

with(lrt.MvG2$table, plot(logCPM,logFC,pch=16,cex=0.2, col="grey50", main = "IgM - IgG"))
with(lrt.MvG2$table, (points(logCPM[MvG2up],logFC[MvG2up],pch=16,col="royalblue3")))
with(lrt.MvG2$table, (points(logCPM[MvG2dn],logFC[MvG2dn],pch=16,col="palegreen4"))) 
abline(0,0,col="grey")

write.csv(sig_genes_DvM2.lrt,"edgeR/With Mouse and MedCV/sig_genes_DvM.csv",
          quote=FALSE, row.names = TRUE)
write.csv(less_sig_genes_DvM2.lrt, "edgeR/With Mouse and MedCV/less_sig_genes_DvM.csv",
          quote = F, row.names = T)
write.csv(top_genes_DvM2.lrt, "edgeR/With Mouse and MedCV/top_genes_DvM.csv",
          quote = F, row.names = T)

write.csv(sig_genes_DvG2.lrt,"edgeR/With Mouse and MedCV/sig_genes_DvG.csv",
          quote=FALSE, row.names = TRUE)
write.csv(less_sig_genes_DvG2.lrt, "edgeR/With Mouse and MedCV/less_sig_genes_DvG.csv",
          quote = F, row.names = T)
write.csv(top_genes_DvG2.lrt, "edgeR/With Mouse and MedCV/top_genes_DvG.csv",
          quote = F, row.names = T)

write.csv(sig_genes_MvG2.lrt,"edgeR/With Mouse and MedCV/sig_genes_MvG.csv",
          quote=FALSE, row.names = TRUE)
write.csv(less_sig_genes_MvG2.lrt, "edgeR/With Mouse and MedCV/less_sig_genes_MvG.csv",
          quote = F, row.names = T)
write.csv(top_genes_MvG2.lrt, "edgeR/With Mouse and MedCV/top_genes_MvG.csv",
          quote = F, row.names = T)

write.csv(less_sig_genes_anova2.lrt,"edgeR/With Mouse and MedCV/less_sig_genes_anova.csv",
          quote = F, row.names = T)


###GLM edgeR analysis (with mouse)
mBCM <- DGEcounts
#Build Model Matrix for analysis
design_mat6 <- model.matrix(~0+design$Mouse+design$Ig)
colnames(design_mat6) <- c("Mouse1","Mouse2","Mouse3","Mouse4","Mouse5","IgG","IgM")

mBCM <- estimateDisp(mBCM, design_mat6)
anova.conm <- makeContrasts(DvG=c(0,0,0,0,0,1,0), DvM=c(0,0,0,0,0,0,1), MvG=IgM-IgG, levels=design_mat6)
fit3 <- glmFit(mBCM, design_mat6)
lrt.anova3 <- glmLRT(fit3, contrast=anova.conm)
lrt.DvM3 <- glmLRT(fit3, contrast = anova.conm[,"DvM"])
lrt.DvG3 <- glmLRT(fit3, contrast = anova.conm[,"DvG"])
lrt.MvG3 <- glmLRT(fit3, contrast = anova.conm[,"MvG"])

top_genes_DvG3.lrt <- topTags(lrt.DvG3, n = Inf, sort.by="P")
top_genes_DvG3.lrt <- top_genes_DvG3.lrt$table
top_genes_DvM3.lrt <- topTags(lrt.DvM3, n = Inf, sort.by="P")
top_genes_DvM3.lrt <- top_genes_DvM3.lrt$table
top_genes_MvG3.lrt <- topTags(lrt.MvG3, n = Inf, sort.by="P")
top_genes_MvG3.lrt <- top_genes_MvG3.lrt$table
top_genes_anova3.lrt <- topTags(lrt.anova3, n = Inf, sort.by="P")
top_genes_anova3.lrt <- top_genes_anova3.lrt$table

sig_genes_DvM3.lrt <- get_sig_genes(top_genes_DvM3.lrt) #16 genes
sig_genes_DvG3.lrt <- get_sig_genes(top_genes_DvG3.lrt) #13 genes
sig_genes_MvG3.lrt <- get_sig_genes(top_genes_MvG3.lrt) #2 genes 
sig_genes_anova3.lrt <- get_sig_genes(top_genes_anova3.lrt) #24 genes 

less_sig_genes_DvM3.lrt <- get_less_sig_genes(top_genes_DvM3.lrt) #66 genes 
less_sig_genes_DvG3.lrt <- get_less_sig_genes(top_genes_DvG3.lrt) #76 genes 
less_sig_genes_MvG3.lrt <- get_less_sig_genes(top_genes_MvG3.lrt) #3 genes 
less_sig_genes_anova3.lrt <- get_less_sig_genes(top_genes_anova3.lrt) #94 genes 

write.csv(sig_genes_DvM3.lrt,"edgeR/With Mouse/sig_genes_DvM.csv",
          quote=FALSE, row.names = TRUE)
write.csv(less_sig_genes_DvM3.lrt, "edgeR/With Mouse/less_sig_genes_DvM.csv",
          quote = F, row.names = T)
write.csv(top_genes_DvM3.lrt, "edgeR/With Mouse/top_genes_DvM.csv",
          quote = F, row.names = T)

write.csv(sig_genes_DvG3.lrt,"edgeR/With Mouse/sig_genes_DvG.csv",
          quote=FALSE, row.names = TRUE)
write.csv(less_sig_genes_DvG3.lrt, "edgeR/With Mouse/less_sig_genes_DvG.csv",
          quote = F, row.names = T)
write.csv(top_genes_DvG3.lrt, "edgeR/With Mouse/top_genes_DvG.csv",
          quote = F, row.names = T)

write.csv(sig_genes_MvG3.lrt,"edgeR/With Mouse/sig_genes_MvG.csv",
          quote=FALSE, row.names = TRUE)
write.csv(less_sig_genes_MvG3.lrt, "edgeR/With Mouse/less_sig_genes_MvG.csv",
          quote = F, row.names = T)
write.csv(top_genes_MvG3.lrt, "edgeR/With Mouse/top_genes_MvG.csv",
          quote = F, row.names = T)

write.csv(less_sig_genes_anova3.lrt,"edgeR/With Mouse/less_sig_genes_anova.csv",
          quote = F, row.names = T)
