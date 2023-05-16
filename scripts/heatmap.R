library(DESeq2)
library(gplots)
library(RColorBrewer) #colour scheme for hierarchical clustering heatmap
library(ggplot2) 
library(ggrepel) #labels for ggplot
library(reshape)
library(tidyr)
library(dplyr)
library(ggdendro)
library(reshape2)
library(apeglm) #for MAplot function in DESEq2
library(genefilter) #top var genes
library(pheatmap)
library(readxl)
library(ggrepel)

### PCA biplot from rlog transformation
my_PCAbiplot <- function(object, ...) {
  .local <- function (object,  X = "PC1", Y = "PC2", intgroup = "condition", ntop = 500, nplot = 10, 
                      returnData = FALSE) {
    
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                       length(rv)))]
    genes <- mcols(object)$GeneSymbol[select]
    
    pca <- prcomp(t(assay(object)[select, ]))
    row.names(pca$rotation) <- genes
    
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    
    if (!all(intgroup %in% names(colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                                 drop = FALSE])
    group <- if (length(intgroup) > 1) {
      factor(apply(intgroup.df, 1, paste, collapse = ":"))
    }
    else {
      colData(object)[[intgroup]]
    }
    d <- data.frame(pca$x, group = group, 
                    intgroup.df, name = colnames(object))
    
    if (returnData) {
      attr(d, "percentVar") <- percentVar
      return(d)
    }
    
    plot <- ggplot(data = d, aes_string(x = X, y = Y, color = "group")) + 
      geom_point(size = 3) + xlab(paste0(X,": ", round(percentVar[which(colnames(pca$x)==X)] * 100), "% variance")) + 
      ylab(paste0("Y",": ", round(percentVar[which(colnames(pca$x)==Y)] * 100), "% variance")) + theme_bw()
    
    plot <- plot + geom_hline(yintercept = 0, size=.2) + geom_vline(xintercept = 0, size=.2)
    
    datapc <- data.frame(varnames=rownames(pca$rotation), pca$rotation)
    mult <- min(
      (max(d[,Y]) - min(d[,Y])/(max(datapc[,Y])-min(datapc[,Y]))),
      (max(d[,X]) - min(d[,X])/(max(datapc[,X])-min(datapc[,X])))
    )
    datapc <- transform(datapc,
                        v1 = .7 * mult * (get(X)),
                        v2 = .7 * mult * (get(Y))
    )
    datapc$distance <- sqrt(datapc$v1^2 + datapc$v2^2)
    #datapc <- datapc[select,]
    
    plot <- plot + coord_equal() + geom_text_repel(data=datapc[1:nplot,], aes(x=v1, y=v2, label=varnames), size = 3, vjust=1, color="black")
    plot <- plot + geom_segment(data=datapc[1:nplot,], aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="blue")
    plot
  }
  .local(object, ...)
}

## MA and volcano plot with LFC shrinkage
DE_plot_func <- function(dds, contrast ,LFC, sign, out) {
  
  res_LFC <- lfcShrink(dds, coef = (contrast), type = "apeglm")
  
  res_LFC <- as.data.frame(res_LFC)
  res_LFC$gene_name <- mcols(dds)$GeneSymbol
  res_LFC[["padj"]][is.na(res_LFC[["padj"]])]<-0.99
  res_LFC[["log2FoldChange"]][is.na(res_LFC[["log2FoldChange"]])]<-0
  res_LFC$sig <- "No"
  res_LFC$sig[abs(res_LFC$log2FoldChange) > LFC & res_LFC$padj < sign] <- "Yes"
  res_LFC$sig <- factor(res_LFC$sig, levels <- c("Yes","No"))
  
  MA_plot<-ggplot(res_LFC,aes(x=baseMean, y=log2FoldChange,col=sig))+
    geom_point(alpha=0.5,size=0.8)+
    geom_hline(yintercept=0,alpha=0.75,col="red")+
    labs(col="")+
    xlab("Base Mean")+
    ylab("Log2 Fold Change") +
    scale_color_manual(labels=c("DE","Not DE"),values=c("red","black"))+ xlim(0,1000)
  theme_bw()
  
  ggsave(plot = MA_plot, filename = paste0(dir,out,"_MAplot.pdf"), 
         height = 7, width = 7)
  
  volcano_plot <- ggplot(res_LFC, aes(x=log2FoldChange, y = -log10(padj))) +
    geom_jitter(aes(color = sig)) + theme_bw() + 
    geom_hline(yintercept=-log10(sign), linetype="dashed", color = "black") + 
    geom_vline (xintercept = c(-LFC, LFC), linetype = "dashed", color = "black") + 
    guides(color = FALSE) + ylab("-log10(Adjusted p-value)") + xlab("Log2 FC")
  
  ggsave(plot = volcano_plot, filename = paste0(dir,out,"_volcano.pdf"), 
         height = 7, width = 7)
}

##### Loading initial data #####
dir <- "/Volumes/groupdir/SUN-DAN-Brickman/Jose/Elena/RNAseq_analysis/nascent_RNAseq_20200526/"
data_path <- paste0(dir, "data/counts_gene-id_gene.txt")
data <- read.delim(data_path,as.is=T)
## Annotations of the data
annot <- data[,1:6]
## Counts. "-11" is to pick the 8h samples (which are the 12 last samples of the file)
counts <- as.matrix(data[,(ncol(data)-11):ncol(data)])

#### Getting gene names from the mouse ensembl data base ####
genes <- annot$Geneid
require(org.Mm.eg.db)
db <- org.Mm.eg.db
ensembl <- suppressWarnings(mapIds(db, keys=genes, keytype="ENSEMBL", column="SYMBOL",multiVals="first"))
genes <- data.frame(SYMBOL = ensembl, ENSEMBL = genes, stringsAsFactors = F)
genes$SYMBOL[is.na(genes$SYMBOL)] = genes$ENSEMBL[is.na(genes$SYMBOL)]

annot$GeneSymbol = genes$SYMBOL

#### Getting sample annotations for DEseq analysis ####

require(stringr)

# extract information (metadata) from sample column names
names = colnames(counts)
names = sapply(names,function(x) gsub("X8h","X_8h",x))
p1 = "(.*)_(.*)_(.*)_(.*)_(.*)_uniq"
colData <- data.frame(str_match(names,p1),stringsAsFactors = F)
colnames(colData) = c("Sample","Dox","Time","Tam","Rep","Rep2")
colData$SampleName = colnames(counts)
colData$TimePt = as.numeric(str_match(colData$Time,"^(.*)h")[,2])

# samples to fix: the Tam field should be empty
# and the values there should move to Dox
## Jose: This does not seem to apply to the 8h time point, i believe it is for the ones before the 2nd treatment
idx = colData$Dox=="X"
colData[idx,"Dox"] = colData[idx,"Tam"]
colData[idx,"Tam"] = NA

colData$SampleNameShort = with(colData,paste(Dox,Tam,Time,Rep2, sep="_"))
colData$ShortLabel = with(colData,paste(Dox,Tam,Time, sep="_"))
colData$ShortLabelTime = with(colData,paste(Dox,Tam, sep="_"))

##########################These are the conditions of the experiments#############################################

################# DOX_4OHT_vs_H2O_4OHT #######################

colData$condition = factor(colData$ShortLabel, levels = c("H2O_4OHT_8h","H2O_ETOH_8h","DOX_ETOH_8h","DOX_4OHT_8h"))
colData$Dox = factor(colData$Dox, levels = c("H2O","DOX"))
colData$Tam = factor(colData$Tam, levels = c("ETOH", "4OHT"))

### DESeq analysis ###
#full condition
dds <- DESeqDataSetFromMatrix(counts, colData, design = ~ condition)

#Adding the extrated annotation that we did at the beginning and add it to the information that we just retrieved
mcols(dds) = cbind(annot,mcols(dds))

## Pre-filtering the dataset
nrow(dds)
dds <- dds[rowSums(counts(dds)) > 10, ]
nrow(dds)

## Results 
dds <- DESeq(dds)
res <- results(dds, name = "condition_DOX_4OHT_8h_vs_H2O_4OHT_8h")
res <- as.data.frame(res)
res$gene_name <- mcols(dds)$Geneid
res[["padj"]][is.na(res[["padj"]])]<-0.99
res[["log2FoldChange"]][is.na(res[["log2FoldChange"]])]<-0

## Up and down regulated genes
up_dox_tam_vs_h2o_tam<-filter(res,log2FoldChange>1, padj<0.05)$gene_name
down_dox_tam_vs_h2o_tam<-filter(res,log2FoldChange<(-1),padj<0.05)$gene_name

up_and_down_dox_tam_vs_h2o_tam <- filter(res, abs(log2FoldChange)> 1, padj < 0.05)

### H2O_ETOH as reference ###
colData$condition = factor(colData$ShortLabel, levels = c("H2O_ETOH_8h","H2O_4OHT_8h","DOX_ETOH_8h","DOX_4OHT_8h"))
colData$Dox = factor(colData$Dox, levels = c("H2O","DOX"))
colData$Tam = factor(colData$Tam, levels = c("ETOH", "4OHT"))

## DESeq analysis ##
#full condition
dds <- DESeqDataSetFromMatrix(counts, colData, design = ~ condition)

#Adding the extrated annotation that we did at the beginning and add it to the information that we just retrieved
mcols(dds) = cbind(annot,mcols(dds))

## Pre-filtering the dataset
nrow(dds)
dds <- dds[rowSums(counts(dds)) > 10, ]
nrow(dds)

### DE analysis 4OHT vs Control

## Up and down regulated pErk regulated genes
dds <- DESeq(dds)
res <- results(dds, name = "condition_H2O_4OHT_8h_vs_H2O_ETOH_8h")
res <- as.data.frame(res)
res$gene_name <- mcols(dds)$Geneid
res[["padj"]][is.na(res[["padj"]])]<-0.99
res[["log2FoldChange"]][is.na(res[["log2FoldChange"]])]<-0

up<-filter(res,log2FoldChange>1, padj<0.05)
down<-filter(res,log2FoldChange<(-1),padj<0.05)

up_and_down_tam_vs_etoh <- filter(res, abs(log2FoldChange)> 1, padj < 0.05)
up_and_down_tam_vs_etoh$DE <- "Up"
up_and_down_tam_vs_etoh$DE[up_and_down_tam_vs_etoh$log2FoldChange < (-1)] <- "Down"

#tam_regulated <- data.frame (genes = rbind(up,down), DE = c(rep("up",nrow(up)), rep("down", nrow(down))))
#write.table(tam_regulated, file = paste0(dir,"results/4OHT_vs_Control_DE_genes.tsv"), sep = "\t",
#            quote = F, col.names = T, row.names = F)

## Hetmap of genes regulated by pErk
#mat <- assay(normTransform(dds))
mat <- assay(rlog(dds))
rownames(mat) <- mcols(dds)$Geneid
mat  <- mat[up_and_down_tam_vs_etoh$gene_name,]
#mat  <- mat - rowMeans(mat)
mat <- mat[,c("DOX_8h_4OHT_a_S28_uniq", "DOX_8h_4OHT_b_S29_uniq", "DOX_8h_4OHT_d_S30_uniq", 
              "H2O_8h_4OHT_a_S25_uniq", "H2O_8h_4OHT_b_S26_uniq", "H2O_8h_4OHT_d_S27_uniq", 
              "DOX_8h_ETOH_a_S22_uniq", "DOX_8h_ETOH_b_S23_uniq", "DOX_8h_ETOH_d_S24_uniq",
              "H2O_8h_ETOH_a_S19_uniq", "H2O_8h_ETOH_b_S20_uniq", "H2O_8h_ETOH_d_S21_uniq")]


#row_clustering <- hclust(dist(mat[,c("H2O_8h_4OHT_a_S25_uniq", "H2O_8h_4OHT_b_S26_uniq", "H2O_8h_4OHT_d_S27_uniq",
#                                     "H2O_8h_ETOH_a_S19_uniq", "H2O_8h_ETOH_b_S20_uniq", "H2O_8h_ETOH_d_S21_uniq",
#                                     "DOX_8h_4OHT_a_S28_uniq", "DOX_8h_4OHT_b_S29_uniq", "DOX_8h_4OHT_d_S30_uniq")]))

row_anno <- data.frame(row.names = rownames(mat), pErk_DE = up_and_down_tam_vs_etoh$DE)

row_anno$pErk_mutant_DE <- "No change"
row_anno$pErk_mutant_DE[row.names(row_anno) %in% intersect(row.names(row_anno), up_dox_tam_vs_h2o_tam)] <- "Up"
row_anno$pErk_mutant_DE[row.names(row_anno) %in% intersect(row.names(row_anno), down_dox_tam_vs_h2o_tam) ] <- "Down"

anno <- as.data.frame(colData(dds)[c("DOX_8h_4OHT_a_S28_uniq", "DOX_8h_4OHT_b_S29_uniq", "DOX_8h_4OHT_d_S30_uniq", 
                                     "H2O_8h_4OHT_a_S25_uniq", "H2O_8h_4OHT_b_S26_uniq", "H2O_8h_4OHT_d_S27_uniq", 
                                     "DOX_8h_ETOH_a_S22_uniq", "DOX_8h_ETOH_b_S23_uniq", "DOX_8h_ETOH_d_S24_uniq",
                                     "H2O_8h_ETOH_a_S19_uniq", "H2O_8h_ETOH_b_S20_uniq", "H2O_8h_ETOH_d_S21_uniq"), c("Dox","Tam")])

colnames(anno) <- c("Oct4", "pErk")
anno$Oct4 <- as.character(anno$Oct4)
anno$Oct4[which(anno$Oct4 == "DOX")] <- "-"
anno$Oct4[which(anno$Oct4 == "H2O")] <- "+"
anno$Oct4 <- as.factor(anno$Oct4)

anno$pErk <- as.character(anno$pErk)
anno$pErk[which(anno$pErk == "4OHT")] <- "+"
anno$pErk[which(anno$pErk == "ETOH")] <- "-"
anno$pErk <- as.factor(anno$pErk)

ann_colors <- list(Oct4 = c("+"="grey45", "-"= "grey75"), pErk = c("+"="grey45", "-"= "grey75"),
                   pErk_DE = c("Up" = "firebrick2", "Down" = "skyblue2"),
                   pErk_mutant_DE = c("Up" = "firebrick2", "Down" = "skyblue2", "No change" = "grey"))

hm<- pheatmap(#filename = paste0(dir,"results/tam_vs_etoh_annotated_heatmap.png"),  
              mat, annotation_col = anno, annotation_colors = ann_colors, 
              annotation_row = row_anno, 
              show_rownames = F, show_colnames = F, 
              cluster_cols = F, cluster_rows = T, 
              main = "Dysregulation of Fgf-Erk dependent genes", scale = "row",
              treeheight_row = 0, treeheight_col = 0, width = 7, height = 7)

## Row annotations
row_anno$Geneid <- rownames(row_anno)

row_anno <- merge(row_anno, mcols(dds)[,c("Geneid","GeneSymbol")], by = "Geneid")
write.table(row_anno, file = paste0(dir,"results/tam_vs_etoh_annotated_heatmap_row_anno.tsv"), row.names = F, col.names = T, sep = "\t", quote = F) 

### Bed files for GREAT

# Tam vs EtOH up
bed_file <- annot[annot$Geneid %in% up_and_down_tam_vs_etoh$gene_name[up_and_down_tam_vs_etoh$DE =="Up"],]
bed_file <- bed_file[,c(2,3,4,7)]
bed_file$Chr <- paste0("chr",bed_file$Chr)
write.table(bed_file, row.names = F, col.names = F, sep = "\t", na = "", quote = F,
            file = paste0(dir,"results/GREAT_Tam_vs_EtOH_Up.bed"))

# Tam vs EtOH down
bed_file <- annot[annot$Geneid %in% up_and_down_tam_vs_etoh$gene_name[up_and_down_tam_vs_etoh$DE =="Down"],]
bed_file <- bed_file[,c(2,3,4,7)]
bed_file$Chr <- paste0("chr",bed_file$Chr)
write.table(bed_file, row.names = F, col.names = F, sep = "\t", na = "", quote = F,
            file = paste0(dir,"results/GREAT_Tam_vs_EtOH_Down.bed"))

# Tam vs EtOH up and Dox-Tam vs Tam up
intersection <- rownames(row_anno)[row_anno$pErk_DE == "Up" & row_anno$pErk_mutant_DE == "Up"] 
bed_file <- annot[annot$Geneid %in% intersection,]
bed_file <- bed_file[,c(2,3,4,7)]
bed_file$Chr <- paste0("chr",bed_file$Chr)
write.table(bed_file, row.names = F, col.names = F, sep = "\t", na = "", quote = F,
            file = paste0(dir,"results/GREAT_Tam_vs_EtOH_up_Dox-Tam_vs_Tam_up.bed"))

# Tam vs EtOH up and Dox-Tam vs Tam down
intersection <- rownames(row_anno)[row_anno$pErk_DE == "Up" & row_anno$pErk_mutant_DE == "Down"] 
bed_file <- annot[annot$Geneid %in% intersection,]
bed_file <- bed_file[,c(2,3,4,7)]
bed_file$Chr <- paste0("chr",bed_file$Chr)
write.table(bed_file, row.names = F, col.names = F, sep = "\t", na = "", quote = F,
            file = paste0(dir,"results/GREAT_Tam_vs_EtOH_up_Dox-Tam_vs_Tam_down.bed"))

# Tam vs EtOH down and Dox-Tam vs Tam up
intersection <- rownames(row_anno)[row_anno$pErk_DE == "Down" & row_anno$pErk_mutant_DE == "Up"] 
bed_file <- annot[annot$Geneid %in% intersection,]
bed_file <- bed_file[,c(2,3,4,7)]
bed_file$Chr <- paste0("chr",bed_file$Chr)
write.table(bed_file, row.names = F, col.names = F, sep = "\t", na = "", quote = F,
            file = paste0(dir,"results/GREAT_Tam_vs_EtOH_down_Dox-Tam_vs_Tam_up.bed"))

# Tam vs EtOH down and Dox-Tam vs Tam down
intersection <- rownames(row_anno)[row_anno$pErk_DE == "Down" & row_anno$pErk_mutant_DE == "Down"] 
bed_file <- annot[annot$Geneid %in% intersection,]
bed_file <- bed_file[,c(2,3,4,7)]
bed_file$Chr <- paste0("chr",bed_file$Chr)
write.table(bed_file, row.names = F, col.names = F, sep = "\t", na = "", quote = F,
            file = paste0(dir,"results/GREAT_Tam_vs_EtOH_down_Dox-Tam_vs_Tam_down.bed"))

