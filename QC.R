library(plyr)
library(ggplot2); theme_set(theme_bw(20) + theme(panel.grid.major = element_blank(),
                                                 panel.grid.minor = element_blank()))
library(edgeR)
library(xlsx)
library(limma)
library(RColorBrewer)
library(gplots)
library(ggthemes)
library(stringr)
library(dplyr)
options(stringsAsFactors = FALSE)

# Get colorblind palette
cb_pal <- colorblind_pal()(8)
my_cb_pal <- cb_pal[c(1,2,4,6,7,8)]

#Read in counts
counts <- read.csv("counts/P141-1_C9E23ANXX_161017_combined_counts.csv", sep = ",")
rownames(counts) <- counts$geneName
counts <- counts[,-1]

#Read in metrics
Metrics <- read.csv("metrics/P141-1_C9E23ANXX_161017_combined_metrics.csv")

#Trim Lib ID's to lib##### (lib plus 5 digit)
colnames(counts) <- substr(colnames(counts), 1, 8)
names(Metrics)[1] <- "Lib.ID" # Reformat libID variable name
Metrics$Lib.ID <- substr(Metrics$Lib.ID, 1, 8)
Metrics.Sort <- arrange(Metrics, fastq_total_reads)

pdf(file="QC.plots/QC_FastqTotalReads.pdf")
barplot(Metrics.Sort$fastq_total_reads/10^6,
        xlab="15 libraries", ylab = "total reads (in millions)", main = "RA_DoD")
hist(Metrics$fastq_total_reads/10^6,
     xlab = "total reads (in millions)",ylab="frequency",main="B Cell Memory Project", col=24)
dev.off()

pdf(file="QC.plots/LibraryCounts.pdf")
barplot(Metrics.Sort$fastq_total_reads[1:20]/10^6, names.arg = Metrics.Sort$Lib.ID[1:20],
        ylab = "total reads (in millions)", main = "B Cell Memory Project",las=2)
dev.off()

ggplot(Metrics, aes(x=MEDIAN_CV_COVERAGE, y=mapped_reads_w_dups)) + geom_point(color="dark green") +
  geom_text(data=subset(Metrics, MEDIAN_CV_COVERAGE > 0.9), aes(y=mapped_reads_w_dups-.01,label=Lib.ID),size=4, vjust=1,hjust=0.5)+
  labs(x = "Median CV Coverage", y = "Percent Alignment") + xlim(0,1.05) + ylim(.2,1) +
  geom_hline(yintercept=0.8) + geom_vline(xintercept=0.9)
ggsave("QC.Plots/PctAlign_v_MedCV.pdf")

ggplot(Metrics, aes(x=fastq_total_reads/10^6, y=mapped_reads_w_dups)) + geom_point(color="red") +
  geom_text(data=subset(Metrics, fastq_total_reads/10^6 < 5), aes(y=mapped_reads_w_dups-.01,label=Lib.ID),size=4, vjust=1,hjust=0.5)+
  #geom_text(aes(y=mapped_reads_w_dups-.01,label=Lib.ID), size=4, vjust=1,hjust=0.5) +
  labs(x = "Total Counts (in millions)", y = "Percent Alignment") + xlim(0,11) + ylim(.2,1) +
  geom_hline(yintercept=0.8) + geom_vline(xintercept=5)
ggsave("QC.Plots/PctAlign_v_TotCounts.pdf")

ggplot(Metrics, aes(x=fastq_total_reads/10^6, y=MEDIAN_CV_COVERAGE)) + geom_point(color="blue") +
  geom_text(data=subset(Metrics, MEDIAN_CV_COVERAGE > 0.9 | fastq_total_reads/10^6 < 5), aes(y=MEDIAN_CV_COVERAGE-.01,label=Lib.ID),size=4, vjust=1,hjust=0.5)+
  labs(x = "Total Counts (in millions)", y = "Median CV Coverage") + xlim(0,11) + ylim(.2,1.05) +
  geom_hline(yintercept=0.9) + geom_vline(xintercept=5)
ggsave("QC.Plots/MedCV_v_TotCounts.pdf")

#Get design info
design <- read.xlsx2("pepper_design.xlsx", sheetIndex = 1)

#Order design according to counts lib 
IndDesign <- match(colnames(counts), design$Lib.ID)
IndDesign <- IndDesign[!is.na(IndDesign)]
design <- design[IndDesign,]
metdes <- merge(design, Metrics, by.x = "Lib.ID", by.y="Lib.ID")
metdes <- metdes[1:15,c(1,64,70,20,2,3)]

ggplot(metdes, aes(x=fastq_total_reads/10^6, y=mapped_reads_w_dups, shape=Ig, color=Mouse)) + geom_point(size=2.5) +
  labs(x = "Total Counts (in millions)", y = "Percent Alignment", shape="Ig", color="Mouse") + xlim(0,11) + ylim(.2,1) +
  geom_hline(yintercept=0.8) + geom_vline(xintercept=5) +  scale_shape_manual(values=c(8,15,17,18,19))+
  scale_color_manual(values=my_cb_pal)
ggsave("QC.Plots/PctAlign_v_TotCounts_byMouseIg.pdf")

ggplot(metdes, aes(x=fastq_total_reads/10^6, y=MEDIAN_CV_COVERAGE, shape=Ig, color=Mouse)) + geom_point(size=2.5) +
  labs(x = "Total Counts (in millions)", y = "Median CV Coverage", shape="Ig", color="Mouse") + xlim(0,11) + ylim(.2,1.05) +
  geom_hline(yintercept=0.9) + geom_vline(xintercept=5) +  scale_shape_manual(values=c(8,15,17,18,19))+
  scale_color_manual(values=my_cb_pal)
ggsave("QC.Plots/MedCV_v_TotCounts_byMouseIg.pdf")

Metrics = subset(Metrics, Metrics$fastq_total_reads > 5*10^6 &
                       Metrics$MEDIAN_CV_COVERAGE <0.9 &
                       Metrics$mapped_reads_w_dups > 0.8) 

counts <- counts[, which(colnames(counts) %in% Metrics$Lib.ID)]

#Reading in HGNC symbols to replace Ensembl gene IDs
geneKey <- read.table("EnsembltoMGI_GRCm38.txt", header = TRUE,sep = "\t",na.strings = "") 
genesWithID <- geneKey[!is.na(geneKey$Associated.Gene.Name),]
#pcGenes <- subset(genesWithHGNC, genesWithHGNC$Gene.type == "protein_coding")
dpGenes <- genesWithID[duplicated(genesWithID$Ensembl.Gene.ID),] #check duplicated ensembl genes
genesWithID <- genesWithID[!duplicated(genesWithID$Ensembl.Gene.ID),] #discard duplicated ensembl genes
#test <- merge(genesWithID, counts, by.x="Ensembl.Gene.ID", by.y ="row.names")
#missing <- subset(counts, !(rownames(counts) %in% test$Ensembl.Gene.ID))
#write.csv(rownames(missing), file = "QC.plots/data_lost_no_geneID.csv")
counts <- merge(genesWithID, counts, by.x="Ensembl.Gene.ID", by.y ="row.names")

#Average counts for duplicated HGNC symbols
rownames(counts) <- counts$Ensembl.Gene.ID
counts <- aggregate( counts, list(counts$Associated.Gene.Name), mean) #will produce warnings when it tries to average names, ignore these
counts$Gene.type = NULL
counts$Ensembl.Gene.ID = NULL
counts$HGNC.symbol = NULL
rownames(counts) <- counts$Group.1
counts <- counts[,-2]

#Find lib columns that corresponded to desired samples
#design <- design[(design$Mouse == "2" | design$Mouse == "3" | design$Mouse == "4" | design$Mouse == "5"),] ##remove mouse 1 data
keepCounts <- colnames(counts) %in% design$Lib.ID
counts <- counts[,keepCounts]
keepDesign <- design$Lib.ID %in% colnames(counts)
design <- design[keepDesign,]

#Order design according to counts lib 
IndDesign <- match(colnames(counts), design$Lib.ID)
IndDesign <- IndDesign[!is.na(IndDesign)]
design <- design[IndDesign,]

#Normalize counts data
libsize_norm_factors <- colSums(counts)/(10^6)
counts_libsize_norm <- as.data.frame(t(t(counts)/libsize_norm_factors))

#Filter to keep genes that have a count of at least one in 10% of libraries
keepRows <- rowSums((counts_libsize_norm) >= 1) >= 0.10*ncol(counts_libsize_norm)
counts <- counts[keepRows,]

#Make DGEList and apply TMM normalization 
DGEcounts <- DGEList(counts = counts)
DGEcounts <- calcNormFactors(DGEcounts) #TMM default

#Apply log2 cpm normalization
v <- voom(DGEcounts)

#Remove outlier libs from ???
#outlier_libs <- as.factor(c("libX", "libY"))
#design <- design[!(design$Lib.ID %in% outlier_libs),]

save(design, DGEcounts, Metrics, v, file = "R_Analysis/counts_design_with_BCMinfo.RData")

#set seed for reproducability
set.seed(123456)

#Run PCA on the v$E log normalized counts
pca = prcomp(t(v$E), center=TRUE, scale=FALSE) #center to change mean to 0, scale to change SD to 1
#To test PCA using MedCV-corrected counts (from Limma analysis)
#pca = prcomp(t(counts_batch_corr5), center = T, scale = F)
#For PCAs with Mouse and MedCV corrected counts,  
#pca = prcomp(t(counts_batch_corr7), center = T, scale = F)
sum_pca = summary(pca)
scores_pca= as.data.frame(pca$x)

# attach sample info to PC scores  
pdatscores <- merge(design, scores_pca, by.x = "Lib.ID", by.y="row.names")

#scree plot, the number of informative PCs = elbow
pdf(file = "QC.plots/ScreePlot.pdf")
plot(pca, type="l") # l = small L
dev.off()

#separate out the number of PCs that explain most of the variation
pc1All = paste("PC1 (", round(100*sum_pca$importance[2, 1], 1),  "%)", sep="")
pc2All = paste("PC2 (", round(100*sum_pca$importance[2, 2], 1),  "%)", sep="")
pc3All = paste("PC3 (", round(100*sum_pca$importance[2, 3], 1),  "%)", sep="")
pc4All = paste("PC4 (", round(100*sum_pca$importance[2, 4], 1),  "%)", sep="")
pc5All = paste("PC5 (", round(100*sum_pca$importance[2, 5], 1),  "%)", sep="")
pc6All = paste("PC6 (", round(100*sum_pca$importance[2, 6], 1),  "%)", sep="")
pc7All = paste("PC7 (", round(100*sum_pca$importance[2, 7], 1),  "%)", sep="")

g <- ggplot(pdatscores, aes(x=PC1, y=PC2, color=factor(Mouse),shape=factor(Ig))) + geom_point(size=2.5) +
  #geom_text(aes(label=Lib.ID), size=4, vjust=1,hjust=0.5) +
  labs(x = pc1All, y = pc2All, color = "Mouse", shape = "Ig")+
  #theme(panel.background = element_rect(fill='gray'))+
  scale_shape_manual(values=c(8,15,17,18,19))+
  scale_color_manual(values=my_cb_pal)+
  theme(legend.key = element_blank())
filename <- "QC.Plots/PCA/PC1_PC2_by_mouse_and_Ig.pdf"
pdf(filename,width=6, height=4, useDingbats = FALSE)
print(g)
dev.off()

h <- ggplot(pdatscores, aes(x=PC1, y=PC3, color=factor(Mouse),shape=factor(Ig))) + geom_point(size=2.5) +
  #geom_text(aes(label=Lib.ID), size=4, vjust=1,hjust=0.5) +
  labs(x = pc1All, y = pc3All, color = "Mouse", shape = "Ig")+
  #theme(panel.background = element_rect(fill='gray'))+
  scale_shape_manual(values=c(8,15,17,18,19))+
  scale_color_manual(values=my_cb_pal)+
  theme(legend.key = element_blank())
filename <- "QC.Plots/PCA/PC1_PC3_by_mouse_and_Ig.pdf"
pdf(filename,width=6, height=4, useDingbats = FALSE)
print(h)
dev.off()

i <- ggplot(pdatscores, aes(x=PC2, y=PC3, color=factor(Mouse),shape=factor(Ig))) + geom_point(size=2.5) +
  #geom_text(aes(label=Lib.ID), size=4, vjust=1,hjust=0.5) +
  labs(x = pc2All, y = pc3All, color = "Mouse", shape = "Ig")+
  #theme(panel.background = element_rect(fill='gray'))+
  scale_shape_manual(values=c(8,15,17,18,19))+
  scale_color_manual(values=my_cb_pal)+
  theme(legend.key = element_blank())
filename <- "QC.Plots/PCA/PC2_PC3_by_mouse_and_Ig.pdf"
pdf(filename,width=6, height=4, useDingbats = FALSE)
print(i)
dev.off()

#Make sure metrics libs are in same order as in counts
lib_order <- match(colnames(v$E), KeepMetrics$Lib.ID)
metrics <- KeepMetrics[lib_order,]

pdatscores<- merge(metrics, pdatscores, by.x = "Lib.ID", by.y = "Lib.ID")

pdatscores2 <- pdatscores
pdatscores2$Mouse <- as.numeric(pdatscores2$Mouse)
pdatscores2$Ig <- as.numeric(as.factor(pdatscores2$Ig))
pdatscores2$Sample.ID <- NULL

pdatscores_numeric <- as.data.frame(apply(pdatscores2,2,as.numeric))
pdatscores_numeric <- pdatscores_numeric[, colSums(is.na(pdatscores_numeric)) != nrow(pdatscores_numeric)]

correlations <- cor(pdatscores_numeric, use = "pairwise.complete.obs")

#Plot a heatmap of correlation of interest
correlations_sub <- correlations[61:70,1:60]
heatmap.col <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)

pdf(file = "QC.Plots/PCA/Top10_PCA_Heatmap.pdf")
heatmap.2(correlations_sub, scale="none",Rowv = FALSE, Colv=FALSE, col=heatmap.col, 
          trace="none",notecol = "black", key=TRUE,
          density.info = "none", keysize = 1, key.title = "correlation", margins = c(11,4),
          srtCol = 60)
dev.off()
##Or export as 8x17 landscape pdf

ggplot(pdatscores, aes(x=PC1, y=PC2, color=PF_INDEL_RATE,shape=Ig)) + geom_point(size=2.5) +
  labs(x = pc1All, y = pc2All, color = "PF_INDEL_RATE", shape = "Ig")+
  scale_shape_manual(values=c(8,15,17,18,19))+
  scale_colour_gradientn(colours = rev(rainbow(n=7,start=0, end=0.70)))+
  theme(legend.key = element_blank())

ggplot(pdatscores, aes(x=PC1, y=PC2, color=MEDIAN_CV_COVERAGE,shape=Ig)) + geom_point(size=2.5) +
  geom_text(data=subset(pdatscores, (PC1 > 100  | PC2 < -100)), aes(y=PC2,label=Lib.ID),size=4, vjust=1,hjust=0.5)+
  labs(x = pc1All, y = pc2All, color = "Median CV Coverage", shape = "Ig")+
  scale_shape_manual(values=c(8,15,17,18,19))+
  scale_colour_gradientn(colours = rev(rainbow(n=7,start=0, end=0.70)))+
  theme(legend.key = element_blank())

ggplot(pdatscores, aes(x=PC1, y=PC2, color=MEDIAN_5PRIME_TO_3PRIME_BIAS,shape=Ig)) + geom_point(size=2.5) +
  geom_text(data=subset(pdatscores, (PC1 > 100  | PC2 < -100)), aes(y=PC2,label=Lib.ID),size=4, vjust=1,hjust=0.5)+
  labs(x = pc1All, y = pc2All, color = "Median 5' to 3' bias", shape = "Ig")+
  scale_shape_manual(values=c(8,15,17,18,19))+
  scale_colour_gradientn(colours = rev(rainbow(n=7,start=0, end=0.70)))+
  theme(legend.key = element_blank())


pdf(file = "QC.Plots/PCA/PC1vPF_INDEL.pdf")
ggplot(pdatscores, aes(x=PC1, y=PF_INDEL_RATE, shape=Ig, color=Mouse)) + geom_point(size=2.5) +
  labs(x = pc1All, y = "PF Indel Rate", color = "Mouse", Shape = "Ig")+
  scale_shape_manual(values=c(8,15,17,18,19))+
  scale_color_manual(values=my_cb_pal)+
  theme(legend.key = element_blank())
dev.off()

ggplot(pdatscores, aes(x=PC1, y=Mouse, shape=Ig, color=Mouse)) + geom_point(size=2.5) +
  labs(x = pc1All, y = "Mouse", color = "Mouse", Shape = "Ig")+
  scale_shape_manual(values=c(8,15,17,18,19))+
  scale_color_manual(values=my_cb_pal)+
  theme(legend.key = element_blank())

ggplot(pdatscores, aes(x=PC1, y=Ig, shape=Ig, color=Mouse)) + geom_point(size=2.5) +
  labs(x = pc1All, y = "Ig", color = "Mouse", Shape = "Ig")+
  scale_shape_manual(values=c(8,15,17,18,19))+
  scale_color_manual(values=my_cb_pal)+
  theme(legend.key = element_blank())

ggplot(pdatscores, aes(x=PC2, y=MEDIAN_5PRIME_TO_3PRIME_BIAS, shape=Ig, color=Mouse)) + geom_point(size=2.5) +
  labs(x = pc2All, y = "Median 5' to 3' Bias", color = "Mouse", Shape = "Ig")+
  scale_shape_manual(values=c(8,15,17,18,19))+
  scale_color_manual(values=my_cb_pal)+
  theme(legend.key = element_blank())
ggsave(filename = "QC.plots/PCA/PC2v5to3Bias.pdf")

ggplot(pdatscores, aes(x=PC6, y=Ig, shape=Ig, color=Mouse)) + geom_point(size=2.5) +
  labs(x = pc6All, y = "Ig", color = "Mouse", Shape = "Ig")+
  scale_shape_manual(values=c(8,15,17,18,19))+
  scale_color_manual(values=my_cb_pal)+
  theme(legend.key = element_blank())

ggplot(pdatscores, aes(x=PC2, y=Ig, shape=Ig, color=Mouse)) + geom_point(size=2.5) +
  labs(x = pc2All, y = "Ig", color = "Mouse", Shape = "Ig")+
  scale_shape_manual(values=c(8,15,17,18,19))+
  scale_color_manual(values=my_cb_pal)+
  theme(legend.key = element_blank())

pdf(file = "QC.Plots/PCA/PC1vPF_INDEL.pdf")
ggplot(pdatscores, aes(x=PC1, y=PF_INDEL_RATE, shape=Ig, color=Mouse)) + geom_point(size=2.5) +
  labs(x = pc1All, y = "PF Indel Rate", color = "Mouse", Shape = "Ig")+
  scale_shape_manual(values=c(8,15,17,18,19))+
  scale_color_manual(values=my_cb_pal)+
  theme(legend.key = element_blank())
dev.off()

#Sex Check (all shoucd be females)
count_matrix <- v$E
XYgenes <- read.delim("XYgenes.txt")
XYgenes <- data.frame(XYgenes)
Ygenes <- XYgenes[XYgenes$Gene.type == "protein_coding" & XYgenes$Chromosome.Name == "Y",]
Xgenes <- XYgenes[XYgenes$Gene.type == "protein_coding" & XYgenes$Chromosome.Name == "X",]

keepCounts <- rownames(count_matrix) %in% Ygenes$Associated.Gene.Name
Ycounts <- count_matrix[keepCounts,]
Ycounts <- t(Ycounts)
Ycounts <- data.frame(Ycounts)
Ycounts$Lib.ID <- rownames(Ycounts)
Ycounts$sex <- rep("female",14)
Ycounts$sample <- design[match(Ycounts$Lib.ID, design$Lib.ID),"Sample.ID"]
Ycounts$MeanY <- rep(0, 14)
ggplot(Ycounts, aes(x=Lib.ID, y=MeanY, color=factor(sex))) + geom_bar(stat = "identity") + facet_grid(~sex) +
  scale_color_brewer()

keepCountsX <- rownames(count_matrix) %in% Xgenes$Associated.Gene.Name
Xcounts <- count_matrix[keepCountsX,]
Xcounts <- t(Xcounts)
Xcounts <- data.frame(Xcounts)
Xcounts$Lib.ID <- rownames(Xcounts)
Xcounts$sex <- rep("female",14)
Xcounts$sample <- design[match(Xcounts$Lib.ID, design$Lib.ID),"Sample.ID"]
Xcounts$MeanX <- rowMeans(Xcounts[,c(1:374)], na.rm = T)
ggplot(Xcounts, aes(x=Lib.ID, y=MeanX, color=factor(sex))) + geom_bar(stat = "identity") + facet_grid(~sex) +
  scale_color_brewer()

#Ig Check
IgG1 <- as.data.frame(v$E["Ighg1",])
IgG2b <- as.data.frame(v$E["Ighg2b",])
IgG2c <- as.data.frame(v$E["Ighg2c",])
IgG3 <- as.data.frame(v$E["Ighg3",])
IgM <- as.data.frame(v$E["Ighm",])
IgD <- as.data.frame(v$E["Zfp318",])
Ig <- data.frame(v$E["Ighg1",], v$E["Ighg2b",], v$E["Ighg2c",], v$E["Ighg3",], v$E["Ighm",], v$E["Cd80",], v$E["Nt5e",], 
                 v$E["Cd38",], v$E["Sdc1",], design$Ig, design$Mouse)
colnames(Ig) <- c("IgG1", "IgG2b", "IgG2c", "IgG3", "IgM", "CD80", "CD73", "CD38", "CD138", "SampleIg", "Mouse")

design <- merge(design, Metrics, by.x = "Lib.ID", by.y = "Lib.ID")
#Figure out if there's a biological skew with indels between Ig-selected categories
ggplot(design, aes(x=Ig, y=PF_INDEL_RATE, color=Mouse)) + geom_point(size=2.5) + labs(x="Ig", y="PF Indel Rate", color = "Mouse") + 
  scale_color_manual(values=my_cb_pal) + theme(legend.key = element_blank())
#Sort genes by correlation with PF_INDEL_RATE
test <- v$E
set <- data.frame(design$Lib.ID, design$PF_INDEL_RATE)
out <- apply(test,1, function(row)cor(row, set$design.PF_INDEL_RATE))
genecor <- data.frame(out)
genecor2 <- genecor[order(-out),,drop = F]
write.csv(genecor2, file="QC.Plots/GenesByIndelCorr.csv")

##Why does median CV help me pull out so many more DE genes?  Can I trust this or is correcting for Median CV bringing back a library size effect?
#Plot Median CV, my samples, and the # of total reads
ggplot(pdatscores, aes(x=PC2, y=MEDIAN_CV_COVERAGE, shape=Ig, color=Mouse)) + geom_point(size=2.5) +
  labs(x = pc2All, y = "Median CV Coverage", color = "Mouse", Shape = "Ig")+
  scale_shape_manual(values=c(8,15,17,18,19))+
  scale_color_manual(values=my_cb_pal)+
  theme(legend.key = element_blank())
ggplot(design, aes(x=MEDIAN_CV_COVERAGE, y=fastq_total_reads)) + geom_point(size=2.5)+ stat_smooth(method=lm)
MCVvTR.lm = lm(MEDIAN_CV_COVERAGE ~ fastq_total_reads, data=design)
summary(MCVvTR.lm) #r^2 =
plot(MCVvTR.lm) #data doesn't look super normally distributed
cor.test(design$MEDIAN_CV_COVERAGE, design$fastq_total_reads) #Pearson corr = 

ggplot(design, aes(x=MEDIAN_5PRIME_TO_3PRIME_BIAS, y=fastq_total_reads)) + geom_point(size=2.5)+ stat_smooth(method=lm)
M53vTR.lm = lm(MEDIAN_5PRIME_TO_3PRIME_BIAS ~ fastq_total_reads, data=design)
summary(M53vTR.lm) #r^2 = 0.4, p = 0.009
plot(M53vTR.lm) #data look normally distributed
cor.test(design$MEDIAN_5PRIME_TO_3PRIME_BIAS, design$fastq_total_reads) #Pearson corr = 0.668, p = 0.009016


ggplot(design, aes(x=Ig, y=fastq_total_reads)) + geom_point(size=2.5)
ggplot(design, aes(x=Ig, y=fastq_total_reads)) + geom_boxplot()
ggplot(design, aes(x=Ig, y=MEDIAN_CV_COVERAGE)) + geom_point(size=2.5)
ggplot(design, aes(x=Ig, y=MEDIAN_CV_COVERAGE)) + geom_boxplot()
#Sort genes by correlation with MEDIAN_CV_COVERAGE
test <- v$E
set <- data.frame(design$Lib.ID, design$MEDIAN_CV_COVERAGE)
MedCVcorr <- apply(test,1, function(row)cor(row, set$design.MEDIAN_CV_COVERAGE))
exprQuintile <- apply(test,2,function(col)ntile(col,5))
genecor <- data.frame(MedCVcorr)
genecor2 <- merge(genecor, exprQuintile, by="row.names")
rownames(genecor2) <- genecor2$Row.names
genecor2 <- genecor2[,-1]
genecor2$avgQuint <- rowMeans(genecor2[,2:15], na.rm=TRUE)
genecor2$G_Quint <- rowMeans(genecor2[,c(2:4,8)], na.rm=TRUE)
genecor2$M_Quint <- rowMeans(genecor2[,c(5,6,9,10,12)], na.rm=TRUE)
genecor2$D_Quint <- rowMeans(genecor2[,c(7,11,13:15)], na.rm=TRUE)
genecor2 <- genecor2[with(genecor2, order(-MedCVcorr)),]
genecor2 <- genecor2[,c(1,16:19)]
write.csv(genecor2, file="QC.Plots/GenesByMedCVcorr.csv")

