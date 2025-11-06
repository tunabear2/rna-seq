library(DESeq2)
library(dplyr)
library(tibble)
library(ggplot2)
library(tximport)
library(readr)
library(apeglm)
library(ggrepel)
library(ComplexHeatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(fgsea)
library(pathview)
library(AnnotationDbi)
library(biomaRt)
library(circlize)
library(msigdbr)
library(writexl)

samples <- c("prom1_TdT_TW_1_rsem","TdT_TW_1_rsem","prom1_TdT_TW_2_rsem","TdT_TW_2_rsem","prom1_TdT_TW_3_rsem","TdT_TW_3_rsem")

rsem_dir <- "C:/Users/rkawk/rsem/rsem_test"
files <- list.files(rsem_dir, pattern = "\\.genes\\.results$", full.names = TRUE)
stopifnot(length(files) > 1)

condition <- factor(rep(c("prom1_TdT_TW","TdT_TW"), times = length(samples)/2),
                    levels = c("prom1_TdT_TW","TdT_TW"))

coldata <- data.frame(sample = samples, condition = condition, row.names = samples)

condition
data.frame(samples, condition)

txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE,
                countsFromAbundance = "no")
txi$length[txi$length == 0] <- 1

dds <- DESeqDataSetFromTximport(txi = txi, colData = coldata, design = ~ condition)

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "TdT_TW", "prom1_TdT_TW"))
res <- lfcShrink(dds, coef = "condition_TdT_TW_vs_prom1_TdT_TW", type = "apeglm")
