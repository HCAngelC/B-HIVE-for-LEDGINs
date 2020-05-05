# load the pachages
library(DOSE)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(bimaRt)
library(tidyr)
library(forcats)

# 1. Retrieve HIV with and without RNA expression
hiv_allsamples <- read.table("hiv_allsamples_jurkatchip.txt", header = T) #hiv_allsamples_jurkatchip.txt is available at xxx.

## SupT1 
SupT1_0ic_gene_S296_woRNA <- hiv_allsamples %>%
dplyr::filter(rep == "S296" & cat != "RE" & rna == 0) %>%
dplyr::select(gene_name, expscore) %>%
na.omit()

SupT1_2ic_gene_S297_woRNA <- hiv_allsamples %>%
dplyr::filter(rep == "S297" & cat != "RE" & rna == 0) %>%
dplyr::select(gene_name, expscore) %>%
na.omit()

SupT1_5ic_gene_S298_woRNA <- hiv_allsamples %>%
dplyr::filter(rep == "S298" & cat != "RE" & rna == 0) %>%
dplyr::select(gene_name, expscore) %>%
na.omit()

SupT1_10ic_gene_S299_woRNA <- hiv_allsamples %>%
dplyr::filter(rep == "S299" & cat != "RE" & rna == 0) %>%
dplyr::select(gene_name, expscore) %>%
na.omit()

## Jurkat  
Jurkat_0ic_gene_S300_woRNA <- hiv_allsamples %>%
dplyr::filter(rep == "S300" & cat != "RE" & rna == 0) %>%
dplyr::select(gene_name, expscore) %>%
na.omit()

Jurkat_2ic_gene_S301_woRNA <- hiv_allsamples %>%
dplyr::filter(rep == "S301" & cat != "RE" & rna == 0) %>%
dplyr::select(gene_name, expscore) %>%
na.omit()

Jurkat_5ic_gene_S302_woRNA <- hiv_allsamples %>%
dplyr::filter(rep == "S302" & cat != "RE" & rna == 0) %>%
dplyr::select(gene_name, expscore) %>%
na.omit()

Jurkat_10ic_gene_S303_woRNA <- hiv_allsamples %>%
dplyr::filter(rep == "S303" & cat != "RE" & rna == 0) %>%
dplyr::select(gene_name, expscore) %>%
na.omit()

# aggregate for each unique gene
## SupT1
SupT1_0ic_gen_S296_woRNA_agg <- aggregate(SupT1_0ic_gen_S296_woRNA$expscore, by = list(SupT1_0ic_gen_S296_woRNA$gene_name), FUN = mean)
names(SupT1_0ic_gen_S296_woRNA_agg) <- c("gene_name", "exp_mean")
head(SupT1_0ic_gen_S296_woRNA_agg)

SupT1_2ic_gen_S297_woRNA_agg <- aggregate(SupT1_2ic_gen_S297_woRNA$expscore, by = list(SupT1_2ic_gen_S297_woRNA$gene_name), FUN = mean)
names(SupT1_2ic_gen_S297_woRNA_agg) <- c("gene_name", "exp_mean")
head(SupT1_2ic_gen_S297_woRNA_agg)

SupT1_5ic_gen_S298_woRNA_agg <- aggregate(SupT1_5ic_gen_S298_woRNA$expscore, by = list(SupT1_5ic_gen_S298_woRNA$gene_name), FUN = mean)
names(SupT1_5ic_gen_S298_woRNA_agg) <- c("gene_name", "exp_mean")
head(SupT1_5ic_gen_S298_woRNA_agg)

SupT1_10ic_gen_S299_woRNA_agg <- aggregate(SupT1_10ic_gen_S299_woRNA$expscore, by = list(SupT1_10ic_gen_S299_woRNA$gene_name), FUN = mean)
names(SupT1_10ic_gen_S299_woRNA_agg) <- c("gene_name", "exp_mean")
head(SupT1_10ic_gen_S299_woRNA_agg)

## Jurkat
Jurkat_0ic_gen_S300_woRNA_agg <- aggregate(Jurkat_0ic_gen_S300_woRNA$expscore, by = list(Jurkat_0ic_gen_S300_woRNA$gene_name), FUN = mean)
names(Jurkat_0ic_gen_S300_woRNA_agg) <- c("gene_name", "exp_mean")
head(Jurkat_0ic_gen_S300_woRNA_agg)

Jurkat_2ic_gen_S301_woRNA_agg <- aggregate(Jurkat_2ic_gen_S301_woRNA$expscore, by = list(Jurkat_2ic_gen_S301_woRNA$gene_name), FUN = mean)
names(Jurkat_2ic_gen_S301_woRNA_agg) <- c("gene_name", "exp_mean")
head(Jurkat_2ic_gen_S301_woRNA_agg)

Jurkat_5ic_gen_S302_woRNA_agg <- aggregate(Jurkat_5ic_gen_S302_woRNA$expscore, by = list(Jurkat_5ic_gen_S302_woRNA$gene_name), FUN = mean)
names(Jurkat_5ic_gen_S302_woRNA_agg) <- c("gene_name", "exp_mean")
head(Jurkat_5ic_gen_S302_woRNA_agg)

Jurkat_10ic_gen_S303_woRNA_agg <- aggregate(Jurkat_10ic_gen_S303_woRNA$expscore, by = list(Jurkat_10ic_gen_S303_woRNA$gene_name), FUN = mean)
names(Jurkat_10ic_gen_S303_woRNA_agg) <- c("gene_name", "exp_mean")
head(Jurkat_10ic_gen_S303_woRNA_agg)

# entrezID
## SupT1
SupT1_0ic_gen_S296_woRNA_agg$gene_name <- as.character(SupT1_0ic_gen_S296_woRNA_agg$gene_name)
SupT1_0ic_gen_S296_woRNA_agg.entrezID <- bitr(SupT1_0ic_gen_S296_woRNA_agg$gene_name, fromType = "SYMBOL", toType =c("ENSEMBL","ENTREZID"), OrgDb = "org.Hs.eg.db")

SupT1_2ic_gen_S297_woRNA_agg$gene_name <- as.character(SupT1_2ic_gen_S297_woRNA_agg$gene_name)
SupT1_2ic_gen_S297_woRNA_agg.entrezID <- bitr(SupT1_2ic_gen_S297_woRNA_agg$gene_name, fromType = "SYMBOL", toType =c("ENSEMBL","ENTREZID"), OrgDb = "org.Hs.eg.db")

SupT1_5ic_gen_S298_woRNA_agg$gene_name <- as.character(SupT1_5ic_gen_S298_woRNA_agg$gene_name)
SupT1_5ic_gen_S298_woRNA_agg.entrezID <- bitr(SupT1_5ic_gen_S298_woRNA_agg$gene_name, fromType = "SYMBOL", toType =c("ENSEMBL","ENTREZID"), OrgDb = "org.Hs.eg.db")

SupT1_10ic_gen_S299_woRNA_agg$gene_name <- as.character(SupT1_10ic_gen_S299_woRNA_agg$gene_name)
SupT1_10ic_gen_S299_woRNA_agg.entrezID <- bitr(SupT1_10ic_gen_S299_woRNA_agg$gene_name, fromType = "SYMBOL", toType =c("ENSEMBL","ENTREZID"), OrgDb = "org.Hs.eg.db")

## Jurkat
Jurkat_0ic_gen_S300_woRNA_agg$gene_name <- as.character(Jurkat_0ic_gen_S300_woRNA_agg$gene_name)
Jurkat_0ic_gen_S300_woRNA_agg.entrezID <- bitr(Jurkat_0ic_gen_S300_woRNA_agg$gene_name, fromType = "SYMBOL", toType =c("ENSEMBL","ENTREZID"), OrgDb = "org.Hs.eg.db")

Jurkat_2ic_gen_S301_woRNA_agg$gene_name <- as.character(Jurkat_2ic_gen_S301_woRNA_agg$gene_name)
Jurkat_2ic_gen_S301_woRNA_agg.entrezID <- bitr(Jurkat_2ic_gen_S301_woRNA_agg$gene_name, fromType = "SYMBOL", toType =c("ENSEMBL","ENTREZID"), OrgDb = "org.Hs.eg.db")

Jurkat_5ic_gen_S302_woRNA_agg$gene_name <- as.character(Jurkat_5ic_gen_S302_woRNA_agg$gene_name)
Jurkat_5ic_gen_S302_woRNA_agg.entrezID <- bitr(Jurkat_5ic_gen_S302_woRNA_agg$gene_name, fromType = "SYMBOL", toType =c("ENSEMBL","ENTREZID"), OrgDb = "org.Hs.eg.db")

Jurkat_10ic_gen_S303_woRNA_agg$gene_name <- as.character(Jurkat_10ic_gen_S303_woRNA_agg$gene_name)
Jurkat_10ic_gen_S303_woRNA_agg.entrezID <- bitr(Jurkat_10ic_gen_S303_woRNA_agg$gene_name, fromType = "SYMBOL", toType =c("ENSEMBL","ENTREZID"), OrgDb = "org.Hs.eg.db")

# KEGG pathway enrichment analysis
## SupT1
kkenrich.SupT1_0ic_gen_S296_woRNA_agg.entrezID <- enrichKEGG(gene = SupT1_0ic_gen_S296_woRNA_agg.entrezID, organism = 'hsa', pvalueCutoff = 1)
kkenrich.SupT1_2ic_gen_S297_woRNA_agg.entrezID <- enrichKEGG(gene = SupT1_2ic_gen_S297_woRNA_agg.entrezID, organism = 'hsa', pvalueCutoff = 1)
kkenrich.SupT1_5ic_gen_S298_woRNA_agg.entrezID <- enrichKEGG(gene = SupT1_5ic_gen_S298_woRNA_agg.entrezID, organism = 'hsa', pvalueCutoff = 1)
kkenrich.SupT1_10ic_gen_S299_woRNA_agg.entrezID <- enrichKEGG(gene = SupT1_10ic_gen_S299_woRNA_agg.entrezID, organism = 'hsa', pvalueCutoff = 1)

## Jurkat
kkenrich.Jurkat_0ic_gen_S300_woRNA_agg.entrezID <- enrichKEGG(gene = Jurkat_0ic_gen_S300_woRNA_agg.entrezID, organism = 'hsa', pvalueCutoff = 1)
kkenrich.Jurkat_2ic_gen_S301_woRNA_agg.entrezID <- enrichKEGG(gene = Jurkat_2ic_gen_S301_woRNA_agg.entrezID, organism = 'hsa', pvalueCutoff = 1)
kkenrich.Jurkat_5ic_gen_S302_woRNA_agg.entrezID <- enrichKEGG(gene = Jurkat_5ic_gen_S302_woRNA_agg.entrezID, organism = 'hsa', pvalueCutoff = 1)
kkenrich.Jurkat_10ic_gen_S303_woRNA_agg.entrezID <- enrichKEGG(gene = Jurkat_10ic_gen_S303_woRNA_agg.entrezID, organism = 'hsa', pvalueCutoff = 1)

########################################################
# import rna-seq fc file
rnaseq_supt1_fc <- read.table("rnaseq_supt1_fc.txt", header = T)
rnaseq_jurkat_fc <- read.table("rnaseq_jurkat_fc.txt", header = T)

## associate rna-seq file with woRNA.entrezID file
SupT1_0ic_gen_S296_woRNA_agg.entrezID.rnaseq <- merge(SupT1_0ic_gen_S296_woRNA_agg.entrezID, rnaseq_supt1_fc, b = "SYMBOL")
SupT1_2ic_gen_S297_woRNA_agg.entrezID.rnaseq <- merge(SupT1_2ic_gen_S297_woRNA_agg.entrezID, rnaseq_supt1_fc, b = "SYMBOL")
SupT1_5ic_gen_S298_woRNA_agg.entrezID.rnaseq <- merge(SupT1_5ic_gen_S298_woRNA_agg.entrezID, rnaseq_supt1_fc, b = "SYMBOL")
SupT1_10ic_gen_S299_woRNA_agg.entrezID.rnaseq <- merge(SupT1_10ic_gen_S299_woRNA_agg.entrezID, rnaseq_supt1_fc, b = "SYMBOL")

Jurkat_0ic_gen_S300_woRNA_agg.entrezID.rnaseq <- merge(Jurkat_0ic_gen_S300_woRNA_agg.entrezID, rnaseq_jurkat_fc, b = "SYMBOL")
Jurkat_2ic_gen_S301_woRNA_agg.entrezID.rnaseq <- merge(Jurkat_2ic_gen_S301_woRNA_agg.entrezID, rnaseq_jurkat_fc, b = "SYMBOL")
Jurkat_5ic_gen_S302_woRNA_agg.entrezID.rnaseq <- merge(Jurkat_5ic_gen_S302_woRNA_agg.entrezID, rnaseq_jurkat_fc, b = "SYMBOL")
Jurkat_10ic_gen_S303_woRNA_agg.entrezID.rnaseq <- merge(Jurkat_10ic_gen_S303_woRNA_agg.entrezID, rnaseq_jurkat_fc, b = "SYMBOL")

### prepare geneList based on endogenous gene expression
#### SupT1
geneList_supt1_0ic_gen_S296_woRNA <- log10(SupT1_0ic_gen_S296_woRNA_agg.entrezID.rnaseq[,6] + 1)
names(geneList_supt1_0ic_gen_S296_woRNA) <- as.character(SupT1_0ic_gen_S296_woRNA_agg.entrezID.rnaseq[,3])
geneList_supt1_0ic_gen_S296_woRNA <- sort(geneList_supt1_0ic_gen_S296_woRNA, decreasing = T)

geneList_supt1_2ic_gen_S297_woRNA <- log10(SupT1_2ic_gen_S297_woRNA_agg.entrezID.rnaseq[,6] + 1)
names(geneList_supt1_2ic_gen_S297_woRNA) <- as.character(SupT1_2ic_gen_S297_woRNA_agg.entrezID.rnaseq[,3])
geneList_supt1_2ic_gen_S297_woRNA <- sort(geneList_supt1_2ic_gen_S297_woRNA, decreasing = T)

geneList_supt1_5ic_gen_S298_woRNA <- log10(SupT1_5ic_gen_S298_woRNA_agg.entrezID.rnaseq[,6] + 1)
names(geneList_supt1_5ic_gen_S298_woRNA) <- as.character(SupT1_5ic_gen_S298_woRNA_agg.entrezID.rnaseq[,3])
geneList_supt1_5ic_gen_S298_woRNA <- sort(geneList_supt1_5ic_gen_S298_woRNA, decreasing = T)

geneList_supt1_10ic_gen_S299_woRNA <- log10(SupT1_10ic_gen_S299_woRNA_agg.entrezID.rnaseq[,6] + 1)
names(geneList_supt1_10ic_gen_S299_woRNA) <- as.character(SupT1_10ic_gen_S299_woRNA_agg.entrezID.rnaseq[,3])
geneList_supt1_10ic_gen_S299_woRNA <- sort(geneList_supt1_10ic_gen_S299_woRNA, decreasing = T)

#### Jurkat
geneList_jurkat_0ic_gen_S300_woRNA <- log10(Jurkat_0ic_gen_S300_woRNA_agg.entrezID.rnaseq[,6] + 1)
names(geneList_jurkat_0ic_gen_S300_woRNA) <- as.character(Jurkat_0ic_gen_S300_woRNA_agg.entrezID.rnaseq[,3])
geneList_jurkat_0ic_gen_S300_woRNA <- sort(geneList_jurkat_0ic_gen_S300_woRNA, decreasing = T)

geneList_jurkat_2ic_gen_S301_woRNA <- log10(Jurkat_2ic_gen_S301_woRNA_agg.entrezID.rnaseq[,6] + 1)
names(geneList_jurkat_2ic_gen_S301_woRNA) <- as.character(Jurkat_2ic_gen_S301_woRNA_agg.entrezID.rnaseq[,3])
geneList_jurkat_2ic_gen_S301_woRNA <- sort(geneList_jurkat_2ic_gen_S301_woRNA, decreasing = T)

geneList_jurkat_5ic_gen_S302_woRNA <- log10(Jurkat_5ic_gen_S302_woRNA_agg.entrezID.rnaseq[,6] + 1)
names(geneList_jurkat_5ic_gen_S302_woRNA) <- as.character(Jurkat_5ic_gen_S302_woRNA_agg.entrezID.rnaseq[,3])
geneList_jurkat_5ic_gen_S302_woRNA <- sort(geneList_jurkat_5ic_gen_S302_woRNA, decreasing = T)

geneList_jurkat_10ic_gen_S303_woRNA <- log10(Jurkat_10ic_gen_S303_woRNA_agg.entrezID.rnaseq[,6] + 1)
names(geneList_jurkat_10ic_gen_S303_woRNA) <- as.character(Jurkat_10ic_gen_S303_woRNA_agg.entrezID.rnaseq[,3])
geneList_jurkat_10ic_gen_S303_woRNA <- sort(geneList_jurkat_10ic_gen_S303_woRNA, decreasing = T)

########################################################
# cnetplot
## SupT1
SupT1_0ic_gen_S296_woRNA_agg.entrezID.kk.reSet <- setReadable(kkenrich.SupT1_0ic_gen_S296_woRNA_agg.entrezID, "org.Hs.eg.db", "ENTREZID")
cnetplot(SupT1_0ic_gen_S296_woRNA_agg.entrezID.kk.reSet, foldChange = geneList_supt1_0ic_gen_S296_woRNA)

SupT1_2ic_gen_S297_woRNA_agg.entrezID.kk.reSet <- setReadable(kkenrich.SupT1_2ic_gen_S297_woRNA_agg.entrezID, "org.Hs.eg.db", "ENTREZID")
cnetplot(SupT1_2ic_gen_S297_woRNA_agg.entrezID.kk.reSet, foldChange = geneList_supt1_2ic_gen_S297_woRNA)

SupT1_5ic_gen_S298_woRNA_agg.entrezID.kk.reSet <- setReadable(kkenrich.SupT1_5ic_gen_S298_woRNA_agg.entrezID, "org.Hs.eg.db", "ENTREZID")
cnetplot(SupT1_5ic_gen_S298_woRNA_agg.entrezID.kk.reSet, foldChange = geneList_supt1_5ic_gen_S298_woRNA)

SupT1_10ic_gen_S299_woRNA_agg.entrezID.kk.reSet <- setReadable(kkenrich.SupT1_10ic_gen_S299_woRNA_agg.entrezID, "org.Hs.eg.db", "ENTREZID")
cnetplot(SupT1_10ic_gen_S299_woRNA_agg.entrezID.kk.reSet, foldChange = geneList_supt1_10ic_gen_S299_woRNA)

## Jurkat
Jurkat_0ic_gen_S300_woRNA_agg.entrezID.kk.reSet <- setReadable(kkenrich.Jurkat_0ic_gen_S300_woRNA_agg.entrezID, "org.Hs.eg.db", "ENTREZID")
cnetplot(Jurkat_0ic_gen_S300_woRNA_agg.entrezID.kk.reSet, foldChange = geneList_jurkat_0ic_gen_S300_woRNA)

Jurkat_2ic_gen_S301_woRNA_agg.entrezID.kk.reSet <- setReadable(kkenrich.Jurkat_2ic_gen_S301_woRNA_agg.entrezID, "org.Hs.eg.db", "ENTREZID")
cnetplot(Jurkat_2ic_gen_S301_woRNA_agg.entrezID.kk.reSet, foldChange = geneList_jurkat_2ic_gen_S301_woRNA)

Jurkat_5ic_gen_S302_woRNA_agg.entrezID.kk.reSet <- setReadable(kkenrich.Jurkat_5ic_gen_S302_woRNA_agg.entrezID, "org.Hs.eg.db", "ENTREZID")
cnetplot(Jurkat_5ic_gen_S302_woRNA_agg.entrezID.kk.reSet, foldChange = geneList_jurkat_5ic_gen_S302_woRNA)

Jurkat_10ic_gen_S303_woRNA_agg.entrezID.kk.reSet <- setReadable(kkenrich.Jurkat_10ic_gen_S303_woRNA_agg.entrezID, "org.Hs.eg.db", "ENTREZID")
cnetplot(Jurkat_10ic_gen_S303_woRNA_agg.entrezID.kk.reSet, foldChange = geneList_jurkat_10ic_gen_S303_woRNA)
