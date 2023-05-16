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
dir <- "/Volumes/danstem/Brickman/Jose/Elena/RNAseq_analysis/nascent_RNAseq_20200526/"
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

### H2O_ETOH as reference ###
colData$condition = factor(colData$ShortLabel, levels = c("H2O_ETOH_8h","H2O_4OHT_8h","DOX_ETOH_8h","DOX_4OHT_8h"))
colData$Dox = factor(colData$Dox, levels = c("H2O","DOX"))
colData$Tam = factor(colData$Tam, levels = c("ETOH", "4OHT"))

### DESeq analysis ###
#full condition
dds <- DESeqDataSetFromMatrix(counts, colData, design = ~ condition)
dds <- DESeq(dds)

#Adding the extrated annotation that we did at the beginning and add it to the information that we just retrieved
mcols(dds) = cbind(annot,mcols(dds))

#### Exploratory analysis ####

## Pre-filtering the dataset
nrow(dds)
dds <- dds[rowSums(counts(dds)) > 10, ]
nrow(dds)

## PCA plot from rlog transformation
rld <- rlog(dds)

rv <- rowVars(assay(rld))
select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
genes <- mcols(rld)$GeneSymbol[select]

pca <- prcomp(t(assay(rld)[select, ]))
row.names(pca$rotation) <- genes

percentVar <- pca$sdev^2/sum(pca$sdev^2)
d <- data.frame(pca$x, pErk = as.character(colData$Tam), Oct4 = as.character(colData$Dox),stringsAsFactors = F)
d$pErk[d$pErk == "4OHT"] <- "+"
d$pErk[d$pErk == "ETOH"] <- "-"
d$Oct4[d$Oct4 == "H2O"] <- "+"
d$Oct4[d$Oct4 == "DOX"] <- "-"

PCA_plot <- ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "Oct4", shape = "pErk")) + 
  geom_point(size = 3) + xlab(paste0("PC1",": ", round(percentVar[which(colnames(pca$x)=="PC1")] * 100), "% variance")) + 
  ylab(paste0("PC2",": ", round(percentVar[which(colnames(pca$x)=="PC2")] * 100), "% variance")) + theme_bw()

PCA_plot <- PCA_plot + geom_hline(yintercept = 0, size=.2) + geom_vline(xintercept = 0, size=.2) + coord_fixed(ratio=1.5) 

#ggsave(plot = PCA_plot, filename = paste0(dir,"results/PCA_plot_8h.pdf"), width = 7, height = 7)

#PCA plot with loadings. Can make plots of other PC. nplot plots the n genes with most variance.
#pdf(paste0(dir,"/sample_myPCA.pdf"), 10,10)
#my_PCAbiplot(rld, X = "PC1", Y = "PC2", nplot = 20) + ggtitle("PCA on rlg transformation and loadings")
#dev.off()

## Sample distances from rlog transformation
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld$Rep2
colnames(sampleDistMatrix) <- rld$Rep2
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
annot_col <- colData[,c("Rep2", "Dox","Tam")]
rownames(annot_col) <- annot_col$Rep2
annot_col$pErk[annot_col$Tam == "4OHT"] <- "+"
annot_col$pErk[annot_col$Tam == "ETOH"] <- "-"
annot_col$Oct4[annot_col$Dox == "H2O"] <- "+"
annot_col$Oct4[annot_col$Dox == "DOX"] <- "-"
annot_col <- annot_col[,-c(1,2,3)]
ann_colors <- list(Oct4 = c("+"="grey45", "-"= "grey75"),pErk = c("+"="grey45", "-"= "grey75"))

#pdf(paste0(dir,"results/sample_distances_8h.pdf"),  width = 7, height = 7)
sample_distances <- pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, annotation_col = annot_col, annotation_colors = ann_colors)

#dev.off()

### DE analysis 4OHT vs Control

#DE_plot_func(dds, "condition_H2O_4OHT_8h_vs_H2O_ETOH_8h", 1, 0.05, "results/4OHT_vs_Control")


## Up and down regulated pErk regulated genes
res <- results(dds, name = "condition_H2O_4OHT_8h_vs_H2O_ETOH_8h")
res <- as.data.frame(res)
res$gene_name <- make.names(mcols(dds)$GeneSymbol, unique = T)
res[["padj"]][is.na(res[["padj"]])]<-0.99
res[["log2FoldChange"]][is.na(res[["log2FoldChange"]])]<-0

up<-filter(res,log2FoldChange>1, padj<0.05)
down<-filter(res,log2FoldChange<(-1),padj<0.05)

up_and_down_tam_vs_etoh <- filter(res, abs(log2FoldChange)> 1, padj < 0.05)$gene_name

tam_regulated <- data.frame (genes = rbind(up,down), DE = c(rep("up",nrow(up)), rep("down", nrow(down))))
#write.table(tam_regulated, file = paste0(dir,"results/4OHT_vs_Control_DE_genes.tsv"), sep = "\t",
#            quote = F, col.names = T, row.names = F)

## Hetmap of genes regulated by pErk
mat <- assay(normTransform(dds))
#mat <- assay(rlog(dds))
rownames(mat) <- res$gene_name
mat  <- mat[up_and_down_tam_vs_etoh,]
#mat  <- mat - rowMeans(mat)
mat <- mat[,c("DOX_8h_4OHT_a_S28_uniq", "DOX_8h_4OHT_b_S29_uniq", "DOX_8h_4OHT_d_S30_uniq", 
              "H2O_8h_4OHT_a_S25_uniq", "H2O_8h_4OHT_b_S26_uniq", "H2O_8h_4OHT_d_S27_uniq", 
              "DOX_8h_ETOH_a_S22_uniq", "DOX_8h_ETOH_b_S23_uniq", "DOX_8h_ETOH_d_S24_uniq",
              "H2O_8h_ETOH_a_S19_uniq", "H2O_8h_ETOH_b_S20_uniq", "H2O_8h_ETOH_d_S21_uniq")]

row_anno <- data.frame(row.names = up_and_down_tam_vs_etoh, DE = c(rep("Up",nrow(up)),
                                                                   rep("Down",nrow(down))))
anno <- as.data.frame(colData(dds)[, c("Dox","Tam")])

colnames(anno) <- c("Oct4", "pErk")
anno$Oct4 <- as.character(anno$Oct4)
anno$Oct4[which(anno$Oct4 == "DOX")] <- "-"
anno$Oct4[which(anno$Oct4 == "H2O")] <- "+"
anno$Oct4 <- as.factor(anno$Oct4)

anno$pErk <- as.character(anno$pErk)
anno$pErk[which(anno$pErk == "4OHT")] <- "+"
anno$pErk[which(anno$pErk == "ETOH")] <- "-"
anno$pErk <- as.factor(anno$pErk)

ann_colors <- list(Oct4 = c("+"="grey45", "-"= "grey75"),pErk = c("+"="grey45", "-"= "grey75"),
                   DE = c("Up" = "firebrick2", "Down" = "skyblue2"))
hm<- pheatmap(mat, annotation_col = anno, annotation_colors = ann_colors, 
              annotation_row = row_anno,
              show_rownames = F, show_colnames = F, 
              cluster_cols = F, cluster_rows = T, 
              main = "Dysregulation of Fgf-Erk dependent genes\nin the absence of Oct4", scale = "row",
              treeheight_row = 0, treeheight_col = 0)

#ggsave(plot = hm, filename = paste0(dir,"results/PCA_plot_8h.pdf"), width = 7, height = 7)




### DE analysis DOX vs Control
## MA and volcano plot with LFC shrinkage

#DE_plot_func(dds, "condition_DOX_ETOH_8h_vs_H2O_ETOH_8h", 1, 0.05, "results/DOX_vs_Control")

## Up and down regulated genes
res <- results(dds, name = "condition_DOX_ETOH_8h_vs_H2O_ETOH_8h")
res <- as.data.frame(res)
res$gene_name <- mcols(dds)$GeneSymbol
res[["padj"]][is.na(res[["padj"]])]<-0.99
res[["log2FoldChange"]][is.na(res[["log2FoldChange"]])]<-0

up<-filter(res,log2FoldChange>1, padj<0.05)
down<-filter(res,log2FoldChange<(-1),padj<0.05)

up_and_down_dox_vs_h2o <- filter(res, abs(log2FoldChange)> 1, padj < 0.05)$gene_name

dox_regulated <- data.frame (genes = rbind(up,down), DE = c(rep("up",nrow(up)), rep("down", nrow(down))))
#write.table(dox_regulated, file = paste0(dir,"results/DOX_vs_Control_DE_genes.tsv"), sep = "\t",
#            quote = F, col.names = T, row.names = F)

### DE analysis DOX-4OHT vs Control
## MA and volcano plot with LFC shrinkage

#DE_plot_func(dds, "condition_DOX_4OHT_8h_vs_H2O_ETOH_8h", 1, 0.05, "results/DOX_4OHT_vs_Control")

## Up and down regulated genes
res <- results(dds, name = "condition_DOX_4OHT_8h_vs_H2O_ETOH_8h")
res <- as.data.frame(res)
res$gene_name <- mcols(dds)$GeneSymbol
res[["padj"]][is.na(res[["padj"]])]<-0.99
res[["log2FoldChange"]][is.na(res[["log2FoldChange"]])]<-0

up<-filter(res,log2FoldChange>1, padj<0.05)
down<-filter(res,log2FoldChange<(-1),padj<0.05)

up_and_down_tam_dox_vs_h2o <- filter(res, abs(log2FoldChange)> 1, padj < 0.05)$gene_name

tam_dox_regulated <- data.frame (genes = rbind(up,down), DE = c(rep("up",nrow(up)), rep("down", nrow(down))))
#write.table(dox_regulated, file = paste0(dir,"results/DOX_4OHT_vs_Control_DE_genes.tsv"), sep = "\t",
#            quote = F, col.names = T, row.names = F)

## data frame for venn diagram
tam_dox_up <- tam_dox_regulated$genes.gene_name[tam_dox_regulated$DE == "up"]
dox_up <- dox_regulated$genes.gene_name[dox_regulated$DE == "up"]
tam_up <- tam_regulated$genes.gene_name[tam_regulated$DE == "up"]
maximum <- max(length(tam_dox_up), length(dox_up), length(tam_up))
up <- data.frame(DOX_TAM = c(tam_dox_up,rep("",maximum - length(tam_dox_up))),
                 DOX = c(dox_up,rep("",maximum - length(dox_up))),
                 TAM = c(tam_up,rep("",maximum - length(tam_up)))
                 )
#write.table(up, file = paste0(dir,"results/chow_ruskey_UP_data.tsv"), sep = "\t",
#            quote = F, col.names = T, row.names = F)


tam_dox_down <- tam_dox_regulated$genes.gene_name[tam_dox_regulated$DE == "down"]
dox_down <- dox_regulated$genes.gene_name[dox_regulated$DE == "down"]
tam_down <- tam_regulated$genes.gene_name[tam_regulated$DE == "down"]
maximum <- max(length(tam_dox_down), length(dox_down), length(tam_down))
down <- data.frame(DOX_TAM = c(tam_dox_down,rep("",maximum - length(tam_dox_down))),
                 DOX = c(dox_down,rep("",maximum - length(dox_down))),
                 TAM = c(tam_down,rep("",maximum - length(tam_down)))
)
#write.table(down, file = paste0(dir,"results/chow_ruskey_DOWN_data.tsv"), sep = "\t",
#            quote = F, col.names = T, row.names = F)

################# DOX_4OHT_vs_H2O_4OHT #######################

colData$condition = factor(colData$ShortLabel, levels = c("H2O_4OHT_8h","H2O_ETOH_8h","DOX_ETOH_8h","DOX_4OHT_8h"))
colData$Dox = factor(colData$Dox, levels = c("H2O","DOX"))
colData$Tam = factor(colData$Tam, levels = c("ETOH", "4OHT"))

### DESeq analysis ###
#full condition
dds <- DESeqDataSetFromMatrix(counts, colData, design = ~ condition)
dds <- DESeq(dds)

#Adding the extrated annotation that we did at the beginning and add it to the information that we just retrieved
mcols(dds) = cbind(annot,mcols(dds))

## Pre-filtering the dataset
nrow(dds)
dds <- dds[rowSums(counts(dds)) > 10, ]
nrow(dds)

#### results ####
## MA and volcano plot after LFC shrinkage
#DE_plot_func(dds, "condition_DOX_4OHT_8h_vs_H2O_4OHT_8h", 1, 0.05, "results/DOX_4OHT_vs_H2O_4OHT")

res <- results(dds, name = "condition_DOX_4OHT_8h_vs_H2O_4OHT_8h")
res <- as.data.frame(res)
res$gene_name <- mcols(dds)$GeneSymbol
res[["padj"]][is.na(res[["padj"]])]<-0.99
res[["log2FoldChange"]][is.na(res[["log2FoldChange"]])]<-0

## Up and down regulated genes
up<-filter(res,log2FoldChange>1, padj<0.05)$gene_name
down<-filter(res,log2FoldChange<(-1),padj<0.05)$gene_name

up_and_down_dox_tam_vs_h2o_tam <- filter(res, abs(log2FoldChange)> 1, padj < 0.05)
gene_order <- order(up_and_down_dox_tam_vs_h2o_tam$log2FoldChange, decreasing = T)
up_and_down_dox_tam_vs_h2o_tam <- up_and_down_dox_tam_vs_h2o_tam$gene_name[gene_order] 




## Gene clustering
#ztam_dox_and_tam <- intersect(up_and_down_dox_tam_vs_h2o_tam, up_and_down_tam_vs_etoh)
#tam_dox_and_tam_only <- setdiff(tam_dox_and_tam, up_and_down_dox_vs_h2o)
#tam_dox_and_tam_and_dox <- intersect(tam_dox_and_tam, up_and_down_dox_vs_h2o)

genes <-mcols(rld)$GeneSymbol %in% tam_dox_and_tam
mat  <- assay(rld)[genes,]
mat  <- mat - rowMeans(mat)
mat <- mat[,c("DOX_8h_4OHT_a_S28_uniq", "DOX_8h_4OHT_b_S29_uniq", "DOX_8h_4OHT_d_S30_uniq", 
              "H2O_8h_4OHT_a_S25_uniq", "H2O_8h_4OHT_b_S26_uniq", "H2O_8h_4OHT_d_S27_uniq", 
              "DOX_8h_ETOH_a_S22_uniq", "DOX_8h_ETOH_b_S23_uniq", "DOX_8h_ETOH_d_S24_uniq",
              "H2O_8h_ETOH_a_S19_uniq", "H2O_8h_ETOH_b_S20_uniq", "H2O_8h_ETOH_d_S21_uniq")]

anno <- as.data.frame(colData(rld)[, c("Dox","Tam")])

colnames(anno) <- c("Oct4", "pErk")
anno$Oct4 <- as.character(anno$Oct4)
anno$Oct4[which(anno$Oct4 == "DOX")] <- "-"
anno$Oct4[which(anno$Oct4 == "H2O")] <- "+"
anno$Oct4 <- as.factor(anno$Oct4)


anno$pErk <- as.character(anno$pErk)
anno$pErk[which(anno$pErk == "4OHT")] <- "+"
anno$pErk[which(anno$pErk == "ETOH")] <- "-"
anno$pErk <- as.factor(anno$pErk)


rownames(mat) <- mcols(rld)$GeneSymbol[genes]

ann_colors <- list(Oct4 = c("+"="grey45", "-"= "grey75"),pErk = c("+"="grey45", "-"= "grey75"))

hm<- pheatmap(mat, annotation_col = anno, annotation_colors = ann_colors, show_rownames = F, show_colnames = F, 
              cluster_cols = F, cutree_row = 5, main = "Dysregulation of Fgf-Erk dependent genes in the absence of Oct4", cellwidth = 40)

#ggsave(file = paste0(dir,"tam_and_dox-tam_heatmap.pdf"), plot = hm, units = "in", width = 10, height = 10)

## different clusters

clusters <- hclust(dist(mat, method = "euclidean"),method = "complete")

cluster_cut <- cutree(clusters, k = 5)

df <- data.frame (name = names(cluster_cut), cluster = cluster_cut)
write.table(df, file = paste0(dir,"genes_per_cluster.tsv"), sep = "\t", row.names = F, quote = F)

genes <- mcols(rld)$GeneSymbol %in% names(cluster_cut[which(cluster_cut == 1)])
mat  <- assay(rld)[genes,]
mat  <- mat - rowMeans(mat)
mat <- mat[,c("DOX_8h_4OHT_a_S28_uniq", "DOX_8h_4OHT_b_S29_uniq", "DOX_8h_4OHT_d_S30_uniq", 
              "H2O_8h_4OHT_a_S25_uniq", "H2O_8h_4OHT_b_S26_uniq", "H2O_8h_4OHT_d_S27_uniq", 
              "DOX_8h_ETOH_a_S22_uniq", "DOX_8h_ETOH_b_S23_uniq", "DOX_8h_ETOH_d_S24_uniq",
              "H2O_8h_ETOH_a_S19_uniq", "H2O_8h_ETOH_b_S20_uniq", "H2O_8h_ETOH_d_S21_uniq")]
anno <- as.data.frame(colData(rld)[, c("Dox","Tam")])
rownames(mat) <- mcols(rld)$GeneSymbol[genes]
pheatmap(mat, annotation_col = anno, cluster_cols = F, main = paste0("Gene cluster ",1), cellwidth = 30)

genes <- mcols(rld)$GeneSymbol %in% names(cluster_cut[which(cluster_cut == 2)])
mat  <- assay(rld)[genes,]
mat  <- mat - rowMeans(mat)
mat <- mat[,c("DOX_8h_4OHT_a_S28_uniq", "DOX_8h_4OHT_b_S29_uniq", "DOX_8h_4OHT_d_S30_uniq", 
              "H2O_8h_4OHT_a_S25_uniq", "H2O_8h_4OHT_b_S26_uniq", "H2O_8h_4OHT_d_S27_uniq", 
              "DOX_8h_ETOH_a_S22_uniq", "DOX_8h_ETOH_b_S23_uniq", "DOX_8h_ETOH_d_S24_uniq",
              "H2O_8h_ETOH_a_S19_uniq", "H2O_8h_ETOH_b_S20_uniq", "H2O_8h_ETOH_d_S21_uniq")]
anno <- as.data.frame(colData(rld)[, c("Dox","Tam")])
rownames(mat) <- mcols(rld)$GeneSymbol[genes]
pheatmap(mat, annotation_col = anno, cluster_cols = F, main = paste0("Gene cluster ",2), cellwidth = 30)

genes <- mcols(rld)$GeneSymbol %in% names(cluster_cut[which(cluster_cut == 3)])
mat  <- assay(rld)[genes,]
mat  <- mat - rowMeans(mat)
mat <- mat[,c("DOX_8h_4OHT_a_S28_uniq", "DOX_8h_4OHT_b_S29_uniq", "DOX_8h_4OHT_d_S30_uniq", 
              "H2O_8h_4OHT_a_S25_uniq", "H2O_8h_4OHT_b_S26_uniq", "H2O_8h_4OHT_d_S27_uniq", 
              "DOX_8h_ETOH_a_S22_uniq", "DOX_8h_ETOH_b_S23_uniq", "DOX_8h_ETOH_d_S24_uniq",
              "H2O_8h_ETOH_a_S19_uniq", "H2O_8h_ETOH_b_S20_uniq", "H2O_8h_ETOH_d_S21_uniq")]
anno <- as.data.frame(colData(rld)[, c("Dox","Tam")])
rownames(mat) <- mcols(rld)$GeneSymbol[genes]
pheatmap(mat, annotation_col = anno, cluster_cols = F, main = paste0("Gene cluster ",3), cellwidth = 30)


genes <- mcols(rld)$GeneSymbol %in% names(cluster_cut[which(cluster_cut == 4)])
mat  <- assay(rld)[genes,]
mat  <- mat - rowMeans(mat)
mat <- mat[,c("DOX_8h_4OHT_a_S28_uniq", "DOX_8h_4OHT_b_S29_uniq", "DOX_8h_4OHT_d_S30_uniq", 
              "H2O_8h_4OHT_a_S25_uniq", "H2O_8h_4OHT_b_S26_uniq", "H2O_8h_4OHT_d_S27_uniq", 
              "DOX_8h_ETOH_a_S22_uniq", "DOX_8h_ETOH_b_S23_uniq", "DOX_8h_ETOH_d_S24_uniq",
              "H2O_8h_ETOH_a_S19_uniq", "H2O_8h_ETOH_b_S20_uniq", "H2O_8h_ETOH_d_S21_uniq")]
anno <- as.data.frame(colData(rld)[, c("Dox","Tam")])
rownames(mat) <- mcols(rld)$GeneSymbol[genes]
pheatmap(mat, annotation_col = anno, cluster_cols = F, main = paste0("Gene cluster ",4), cellwidth = 30, cellheight = 20)


genes <- mcols(rld)$GeneSymbol %in% names(cluster_cut[which(cluster_cut == 5)])
mat  <- assay(rld)[genes,]
mat  <- mat - rowMeans(mat)
mat <- mat[,c("DOX_8h_4OHT_a_S28_uniq", "DOX_8h_4OHT_b_S29_uniq", "DOX_8h_4OHT_d_S30_uniq", 
              "H2O_8h_4OHT_a_S25_uniq", "H2O_8h_4OHT_b_S26_uniq", "H2O_8h_4OHT_d_S27_uniq", 
              "DOX_8h_ETOH_a_S22_uniq", "DOX_8h_ETOH_b_S23_uniq", "DOX_8h_ETOH_d_S24_uniq",
              "H2O_8h_ETOH_a_S19_uniq", "H2O_8h_ETOH_b_S20_uniq", "H2O_8h_ETOH_d_S21_uniq")]
anno <- as.data.frame(colData(rld)[, c("Dox","Tam")])
rownames(mat) <- mcols(rld)$GeneSymbol[genes]
pheatmap(mat, annotation_col = anno, cluster_cols = F, main = paste0("Gene cluster ",5), cellwidth = 30)


#### Heatmap of up and down regulated genes for H2O_4OHT vs H2O_ETOH condition

genes <-mcols(rld)$GeneSymbol %in% up_and_down_tam_vs_etoh
mat  <- assay(rld)[genes,]
mat  <- mat - rowMeans(mat)
mat <- mat[,c("DOX_8h_4OHT_a_S28_uniq", "DOX_8h_4OHT_b_S29_uniq", "DOX_8h_4OHT_d_S30_uniq", 
              "H2O_8h_4OHT_a_S25_uniq", "H2O_8h_4OHT_b_S26_uniq", "H2O_8h_4OHT_d_S27_uniq", 
              "DOX_8h_ETOH_a_S22_uniq", "DOX_8h_ETOH_b_S23_uniq", "DOX_8h_ETOH_d_S24_uniq",
              "H2O_8h_ETOH_a_S19_uniq", "H2O_8h_ETOH_b_S20_uniq", "H2O_8h_ETOH_d_S21_uniq")]
anno <- as.data.frame(colData(rld)[, c("Dox","Tam")])
rownames(mat) <- mcols(rld)$GeneSymbol[genes]

pheatmap(mat, annotation_col = anno, cluster_cols = F, main = "Differentially expressed genes by 4OHT (Tamoxifen)", cellwidth = 30, show_rownames = F)
