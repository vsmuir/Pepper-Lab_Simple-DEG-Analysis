########################################
#Head toward heatmaps for less sig genes
#(With Batch Correction)
########################################
normCounts_all <- counts_batch_corr7
normCounts_IgD <- counts_batch_corr7[,design$Ig == "D"]
normCounts_IgM <- counts_batch_corr7[,design$Ig == "M"]
normCounts_IgG <- counts_batch_corr7[,design$Ig == "G"]
design_all <- design
design_IgD <- design[design$Ig == "D",] 
design_IgG <- design[design$Ig == "G",] 
design_IgM <- design[design$Ig == "M",] 

#Build matrix with all DE genes
DE_IgD <- union(rownames(less_sig_genes_DvM2.lrt), rownames(less_sig_genes_DvG2.lrt))
DE_IgM <- union(rownames(less_sig_genes_DvM2.lrt), rownames(less_sig_genes_MvG2.lrt))  
DE_IgG <- union(rownames(less_sig_genes_DvG2.lrt), rownames(less_sig_genes_MvG2.lrt)) 
DE_all <- union(rownames(less_sig_genes_DvG2.lrt),union(rownames(less_sig_genes_DvM2.lrt), rownames(less_sig_genes_MvG2.lrt)))
DE_counts_IgD <- normCounts_IgD[rownames(normCounts_IgD) %in% DE_IgD,]
DE_counts_IgM <- normCounts_IgM[rownames(normCounts_IgM) %in% DE_IgM,]
DE_counts_IgG <- normCounts_IgG[rownames(normCounts_IgG) %in% DE_IgG,]
DE_counts_all <- normCounts_all[rownames(normCounts_all) %in% DE_all,]

#Heatmap of counts normalized 0 to 1 in each row
norm_DE_counts_IgD <- t(apply(DE_counts_IgD,1,function(x)(x-min(x))/(max(x)-min(x))))
norm_DE_counts_IgM <- t(apply(DE_counts_IgM,1,function(x)(x-min(x))/(max(x)-min(x))))
norm_DE_counts_IgG <- t(apply(DE_counts_IgG,1,function(x)(x-min(x))/(max(x)-min(x))))
norm_DE_counts_all <- t(apply(DE_counts_all,1,function(x)(x-min(x))/(max(x)-min(x))))

#Set up colors for heatmaps
my.col <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)

#Subset data by Ig, order genes according to contrast
M_up_vs_DG <- intersect(rownames(less_sig_genes_MvG2.lrt[less_sig_genes_MvG2.lrt$logFC > 0,]),
                        rownames(less_sig_genes_DvM2.lrt[(less_sig_genes_DvM2.lrt$logFC > 0),]))
MG_up_vs_D <- intersect(rownames(less_sig_genes_DvG2.lrt[(less_sig_genes_DvG2.lrt$logFC > 0),]),
                        rownames(less_sig_genes_DvM2.lrt[(less_sig_genes_DvM2.lrt$logFC > 0),]))
M_up_vs_D <- setdiff(rownames(less_sig_genes_DvM2.lrt[(less_sig_genes_DvM2.lrt$logFC > 0),]), union(MG_up_vs_D, M_up_vs_DG))
G_up_vs_D <- setdiff(rownames(less_sig_genes_DvG2.lrt[(less_sig_genes_DvG2.lrt$logFC > 0),]), MG_up_vs_D)
M_down_vs_D <- rownames(less_sig_genes_DvM2.lrt[(less_sig_genes_DvM2.lrt$logFC < 0),])
G_down_vs_D <- rownames(less_sig_genes_DvG2.lrt[(less_sig_genes_DvG2.lrt$logFC < 0),])

#Column coloring for heatmap
Ig_labs <- ifelse (design$Ig == "D", "indianred",
                  ifelse(design$Ig == "M", "royalblue3",
                        ifelse(design$Ig == "G", "palegreen4", 0)))

mouse_labs <- ifelse(design$Mouse == "1","gray8",
                     ifelse(design$Mouse == "2","gray28",
                            ifelse(design$Mouse == "3","gray48",
                                   ifelse(design$Mouse == "4","gray68",
                                          ifelse(design$Mouse == "5","gray88", 0)))))
col_lab_matrix <- cbind(Ig_labs, mouse_labs)
#Row coloring for heatmap
row_labs_updown <- ifelse(rownames(norm_DE_counts_all) %in% M_up_vs_DG, "mediumblue",
                          ifelse(rownames(norm_DE_counts_all) %in% MG_up_vs_D, "lightseagreen", 
                                 ifelse(rownames(norm_DE_counts_all) %in% M_up_vs_D, "royalblue4",
                                        ifelse(rownames(norm_DE_counts_all) %in% G_up_vs_D, "forestgreen",
                                               ifelse(rownames(norm_DE_counts_all) %in% M_down_vs_D, "coral3",
                                                      ifelse(rownames(norm_DE_counts_all) %in% G_down_vs_D, "palevioletred3",0))))))
up_in_D <- ifelse(rownames(norm_DE_counts_all) %in% union(rownames(less_sig_genes_DvG2.lrt[(less_sig_genes_DvG2.lrt$logFC < 0),]),
                                                         rownames(less_sig_genes_DvM2.lrt[(less_sig_genes_DvM2.lrt$logFC < 0),])), 
                                                        "lightpink1", 0)
up_in_M <-ifelse(rownames(norm_DE_counts_all) %in% union(rownames(less_sig_genes_MvG2.lrt[less_sig_genes_MvG2.lrt$logFC > 0,]),
                                                         rownames(less_sig_genes_DvM2.lrt[(less_sig_genes_DvM2.lrt$logFC > 0),])), 
                                                        "skyblue3", 0)
up_in_G <-ifelse(rownames(norm_DE_counts_all) %in% union(rownames(less_sig_genes_MvG2.lrt[less_sig_genes_MvG2.lrt$logFC < 0,]),
                                                         rownames(less_sig_genes_DvG2.lrt[(less_sig_genes_DvG2.lrt$logFC > 0),])), 
                                                        "palegreen", 0)
row_lab_matrix <- cbind(up_in_D, up_in_M, up_in_G)
#Ordering components of the heatmap so that labels will be correct
normCounts_byIg <- norm_DE_counts_all[,order(Ig_labs)]
normCounts_byMouse <- norm_DE_counts_all[,order(mouse_labs)]
Ig_labs_byIg <- Ig_labs[order(Ig_labs)]
mouse_labs_byIg <- mouse_labs[order(Ig_labs)]
Ig_labs_byMouse <- Ig_labs[order(mouse_labs)]

row_labs_updown_roworder <- row_labs_updown[order(row_labs_updown, decreasing = TRUE)]
normCounts_byIg_roworder <- normCounts_byIg[order(row_labs_updown, decreasing = TRUE),]
normCounts_byMouse_roworder <- normCounts_byMouse[order(row_labs_updown, decreasing = TRUE),]

#Drawing heatmaps for the simple model
pdf(file="edgeR/With Mouse and MedCV/Heatmaps/IgComp_1.pdf")
heatmap.2(norm_DE_counts_all, dendrogram = "none", scale="none", col=my.col,
          ColSideColors = Ig_labs,RowSideColors = row_labs_updown,trace = "none", keysize = 1,
          key=TRUE,density.info="none", margins = c(6,9), main = "Plasmodium-specific Memory B Cells", xlab = "Ig")
dev.off()

heatmap.plus(norm_DE_counts_all, keep.dendro = TRUE, scale="none", col=my.col,
          ColSideColors = col_lab_matrix,RowSideColors = row_lab_matrix,
          margins = c(7,9), main = "Plasmodium-specific Memory B Cells")

pdf(file="edgeR/With Mouse and MedCV/Heatmaps/IgComp_byIg.pdf")
heatmap.2(normCounts_byIg, Colv = FALSE, dendrogram = "none", scale="none", col=my.col, keysize = 1,
          ColSideColors = Ig_labs_byIg, RowSideColors = row_labs_updown, trace = "none",
          key=TRUE,density.info="none", margins = c(6,9), main = "Plasmodium-specific Memory B Cells", xlab = "Ig")
dev.off()


pdf(file="edgeR/With Mouse and MedCV/Heatmaps/IgComp_byMouse.pdf")
heatmap.2(normCounts_byIg, dendrogram = "none", scale="none", col=my.col, keysize = 1,
          ColSideColors = mouse_labs_byIg,RowSideColors = row_labs_updown,trace = "none",
          key=TRUE,density.info="none", margins = c(6,9), main = "Plasmodium-specific Memory B Cells", xlab = "Mouse")
dev.off()

