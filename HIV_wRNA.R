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
SupT1_0ic_gene_S296_wRNA <- hiv_allsamples %>%
dplyr::filter(rep == "S296" & cat != "RE" & rna != 0) %>%
dplyr::select(gene_name, expscore) %>%
na.omit()

SupT1_2ic_gene_S297_wRNA <- hiv_allsamples %>%
dplyr::filter(rep == "S297" & cat != "RE" & rna != 0) %>%
dplyr::select(gene_name, expscore) %>%
na.omit()

SupT1_5ic_gene_S298_wRNA <- hiv_allsamples %>%
dplyr::filter(rep == "S298" & cat != "RE" & rna != 0) %>%
dplyr::select(gene_name, expscore) %>%
na.omit()

SupT1_10ic_gene_S299_wRNA <- hiv_allsamples %>%
dplyr::filter(rep == "S299" & cat != "RE" & rna != 0) %>%
dplyr::select(gene_name, expscore) %>%
na.omit()

## Jurkat  
Jurkat_0ic_gene_S300_wRNA <- hiv_allsamples %>%
dplyr::filter(rep == "S300" & cat != "RE" & rna != 0) %>%
dplyr::select(gene_name, expscore) %>%
na.omit()

Jurkat_2ic_gene_S301_wRNA <- hiv_allsamples %>%
dplyr::filter(rep == "S301" & cat != "RE" & rna != 0) %>%
dplyr::select(gene_name, expscore) %>%
na.omit()

Jurkat_5ic_gene_S302_wRNA <- hiv_allsamples %>%
dplyr::filter(rep == "S302" & cat != "RE" & rna != 0) %>%
dplyr::select(gene_name, expscore) %>%
na.omit()

Jurkat_10ic_gene_S303_wRNA <- hiv_allsamples %>%
dplyr::filter(rep == "S303" & cat != "RE" & rna != 0) %>%
dplyr::select(gene_name, expscore) %>%
na.omit()

# aggregate for each unique gene
## SupT1
SupT1_0ic_gen_S296_wRNA_agg <- aggregate(SupT1_0ic_gen_S296_wRNA$expscore, by = list(SupT1_0ic_gen_S296_wRNA$gene_name), FUN = mean)
names(SupT1_0ic_gen_S296_wRNA_agg) <- c("gene_name", "exp_mean")
head(SupT1_0ic_gen_S296_wRNA_agg)

SupT1_2ic_gen_S297_wRNA_agg <- aggregate(SupT1_2ic_gen_S297_wRNA$expscore, by = list(SupT1_2ic_gen_S297_wRNA$gene_name), FUN = mean)
names(SupT1_2ic_gen_S297_wRNA_agg) <- c("gene_name", "exp_mean")
head(SupT1_2ic_gen_S297_wRNA_agg)

SupT1_5ic_gen_S298_wRNA_agg <- aggregate(SupT1_5ic_gen_S298_wRNA$expscore, by = list(SupT1_5ic_gen_S298_wRNA$gene_name), FUN = mean)
names(SupT1_5ic_gen_S298_wRNA_agg) <- c("gene_name", "exp_mean")
head(SupT1_5ic_gen_S298_wRNA_agg)

SupT1_10ic_gen_S299_wRNA_agg <- aggregate(SupT1_10ic_gen_S299_wRNA$expscore, by = list(SupT1_10ic_gen_S299_wRNA$gene_name), FUN = mean)
names(SupT1_10ic_gen_S299_wRNA_agg) <- c("gene_name", "exp_mean")
head(SupT1_10ic_gen_S299_wRNA_agg)

## Jurkat
Jurkat_0ic_gen_S300_wRNA_agg <- aggregate(Jurkat_0ic_gen_S300_wRNA$expscore, by = list(Jurkat_0ic_gen_S300_wRNA$gene_name), FUN = mean)
names(Jurkat_0ic_gen_S300_wRNA_agg) <- c("gene_name", "exp_mean")
head(Jurkat_0ic_gen_S300_wRNA_agg)

Jurkat_2ic_gen_S301_wRNA_agg <- aggregate(Jurkat_2ic_gen_S301_wRNA$expscore, by = list(Jurkat_2ic_gen_S301_wRNA$gene_name), FUN = mean)
names(Jurkat_2ic_gen_S301_wRNA_agg) <- c("gene_name", "exp_mean")
head(Jurkat_2ic_gen_S301_wRNA_agg)

Jurkat_5ic_gen_S302_wRNA_agg <- aggregate(Jurkat_5ic_gen_S302_wRNA$expscore, by = list(Jurkat_5ic_gen_S302_wRNA$gene_name), FUN = mean)
names(Jurkat_5ic_gen_S302_wRNA_agg) <- c("gene_name", "exp_mean")
head(Jurkat_5ic_gen_S302_wRNA_agg)

Jurkat_10ic_gen_S303_wRNA_agg <- aggregate(Jurkat_10ic_gen_S303_wRNA$expscore, by = list(Jurkat_10ic_gen_S303_wRNA$gene_name), FUN = mean)
names(Jurkat_10ic_gen_S303_wRNA_agg) <- c("gene_name", "exp_mean")
head(Jurkat_10ic_gen_S303_wRNA_agg)

