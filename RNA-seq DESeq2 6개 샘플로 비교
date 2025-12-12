# DESeq2 사용을 위해 RSTUDIO 설치하기
# Rtools 다운로드 하기
https://cran.r-project.org/bin/windows/Rtools/
R 버전에 맞게 다운로드

# R에서 진행

# Bioconductor의 core packages 설치

# DESeq2 과정
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")  #DESeq2 패키지 설치

install.packages("ggplot2") # volcano 등 기본 그래프
install.packages("dplyr")  # 전처리 편의
install.packages("tibble") # 전처리 편의
install.packages("data.table") # 대용량 빠른 처리
install.packages(c("ggrepel","patchwork","cowplot")) # volcano 라벨 겹치 방지와 그림 합치기/레이아웃용
install.packages("msigdbr")
install.packages("magick")
install.packages("ragg")

BiocManager::install("pathview") # KEGG pathway 그림에 발현값 오버레이
BiocManager::install("apeglm") # log2FoldChange를 안정화 하는데 사용. volcano 품질이 올라감.
a
BiocManager::install("clusterProfiler") # GO, KEGG, GSEA 분석을 위한 도구
a
BiocManager::install("org.Hs.eg.db") # 인간 유전자 annotation DB
a
BiocManager::install("AnnotationDbi") # gene ID 변환용 도구
a
BiocManager::install("DOSE") # enrichment 분석용 내부 엔진 도구
a
BiocManager::install("fgsea") # GSEA 전용 패키지
a
BiocManager::install("enrichplot") # GO/KEGG/GSEA 결과 시각화용
a
BiocManager::install("pathview") # KEGG pathway 그림에 발현값 오버레이
a
BiocManager::install("ComplexHeatmap") # Heatmap 시각화
a
BiocManager::install("tximport") # rsem을 통해 생성된 데이터를 DESeq2 분석을 위한 table 데이터로 바꾸는 tool 한마디로 genes.results 불러오기
a
BiocManager::install("readr") # 나중에 csv,tsv 생성을 위한 tool
a
BiocManager::install("biomaRt") # Gene id symbol 설정.

library(DESeq2)
library(dplyr)
library(tibble)
library(ggplot2)
library(tximport)
library(readr)
library(writexl)
library(readxl)
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
library(magick)
library(ragg)
-----------------------------------------------------------------------------------------------------------------------
(2) S1 vs S2

rsem_dir <- "/Users/dwyun/rsem/my_data"
files <- list.files(rsem_dir, pattern = "\\.genes\\.results$", full.names = TRUE)
stopifnot(length(files) > 1)

samples <- sub("\\.genes\\.results$", "", basename(files))
names(files) <- samples

group <- ifelse(grepl("prom1_TdT_TW", samples), "prom1_TdT_TW", "TdT_TW")
condition <- factor(group, levels = c("TdT_TW","prom1_TdT_TW"))

coldata <- data.frame(sample = samples, condition = condition, row.names = samples)

condition
data.frame(samples, condition)

txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE,
                countsFromAbundance = "no")
txi$length[txi$length == 0] <- 1

dds <- DESeqDataSetFromTximport(txi = txi, colData = coldata, design = ~ condition)

dds <- DESeq(dds)

res <- results(dds,contrast = c("condition", "prom1_TdT_TW", "TdT_TW"))
rownames(res) <- sub("\\.\\d+$", "", rownames(res))

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

map_orgdb <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(res_df$gene_id),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID")
) %>% 
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

res_df <- res_df %>%
  dplyr::left_join(map_orgdb, by = c("gene_id" = "ENSEMBL"))

res_df <- res_df %>%
  mutate(
    SYMBOL = coalesce(na_if(SYMBOL, ""), gene_id),
    ENTREZID = coalesce(na_if(as.character(ENTREZID), ""), gene_id)
  )
res_df <- res_df %>%
  dplyr::select(gene_id, SYMBOL, ENTREZID, everything())

write_xlsx(as.data.frame(res_df), path = "/Users/dwyun/rsem/my_data/DESeq2 sub/S1 vs S2.xlsx")

# 필터링 값에 따라 필터링 된 gene list 뽑을거면.

lfc_cut <- 1
fdr_cut <- 0.05

sig_df <- res_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    abs(log2FoldChange) >= lfc_cut
  )
write_xlsx(sig_df,
           path = "/Users/dwyun/rsem/my_data/DESeq2 sub/S1 vs S2 sig_genes.xlsx")

# up gene이나 down gene만 뽑고 싶을 때
sig_up <- res_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    log2FoldChange >= lfc_cut
  )

sig_down <- res_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    log2FoldChange <= -lfc_cut
  )
------------------------------------------------------------------------------------------------------
(1) S1 vs S3
setwd("/Users/dwyun/rsem/my_data")
rsem_dir <- "/Users/dwyun/rsem/my_data"
files <- list.files(rsem_dir, pattern = "\\.genes\\.results$", full.names = TRUE)
stopifnot(length(files) > 1)

samples <- sub("\\.genes\\.results$", "", basename(files))
names(files) <- samples

group <- ifelse(grepl("prom1_TdT_TW", samples), "prom1_TdT_TW", "prom1_TdT")
condition <- factor(group, levels = c("prom1_TdT","prom1_TdT_TW"))

coldata <- data.frame(sample = samples, condition = condition, row.names = samples)

condition
data.frame(samples, condition)

txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE,
                countsFromAbundance = "no")
txi$length[txi$length == 0] <- 1

dds <- DESeqDataSetFromTximport(txi = txi, colData = coldata, design = ~ condition)

dds <- DESeq(dds)

res <- results(dds,contrast = c("condition", "prom1_TdT_TW", "prom1_TdT"))
rownames(res) <- sub("\\.\\d+$", "", rownames(res))

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

map_orgdb <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(res_df$gene_id),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID")
) %>% 
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

res_df <- res_df %>%
  dplyr::left_join(map_orgdb, by = c("gene_id" = "ENSEMBL"))

res_df <- res_df %>%
  mutate(
    SYMBOL = coalesce(na_if(SYMBOL, ""), gene_id),
    ENTREZID = coalesce(na_if(as.character(ENTREZID), ""), gene_id)
  )

res_df <- res_df %>%
  dplyr::select(gene_id, SYMBOL, ENTREZID, everything())

write_xlsx(as.data.frame(res_df), path = "/Users/dwyun/rsem/my_data/DESeq2 sub/S1 vs S3.xlsx")

# 필터링 값에 따라 필터링 된 gene list 뽑을거면.

lfc_cut <- 1
fdr_cut <- 0.05

sig_df <- res_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    abs(log2FoldChange) >= lfc_cut
  )
write_xlsx(sig_df,
           path = "/Users/dwyun/rsem/my_data/DESeq2 sub/S1 vs S3 sig_genes.xlsx")

# up gene이나 down gene만 뽑고 싶을 때
sig_up <- res_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    log2FoldChange >= lfc_cut
  )

sig_down <- res_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    log2FoldChange <= -lfc_cut
  )

--------------------------------------------------------------------------------------------------------------
 (3) S1 vs S4

setwd("/Users/dwyun/rsem/my_data")
rsem_dir <- "/Users/dwyun/rsem/my_data"
files <- list.files(rsem_dir, pattern = "\\.genes\\.results$", full.names = TRUE)
stopifnot(length(files) > 1)

samples <- sub("\\.genes\\.results$", "", basename(files))
names(files) <- samples

group <- ifelse(grepl("prom1_TdT_TW", samples), "prom1_TdT_TW", "TdT")
condition <- factor(group, levels = c("TdT","prom1_TdT_TW"))

coldata <- data.frame(sample = samples, condition = condition, row.names = samples)

condition
data.frame(samples, condition)

txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE,
                countsFromAbundance = "no")
txi$length[txi$length == 0] <- 1

dds <- DESeqDataSetFromTximport(txi = txi, colData = coldata, design = ~ condition)

dds <- DESeq(dds)

res <- results(dds,contrast = c("condition", "prom1_TdT_TW", "TdT"))
rownames(res) <- sub("\\.\\d+$", "", rownames(res))

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

map_orgdb <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(res_df$gene_id),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID")
) %>% 
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

res_df <- res_df %>%
  dplyr::left_join(map_orgdb, by = c("gene_id" = "ENSEMBL"))

res_df <- res_df %>%
  mutate(
    SYMBOL = coalesce(na_if(SYMBOL, ""), gene_id),
    ENTREZID = coalesce(na_if(as.character(ENTREZID), ""), gene_id)
  )

res_df <- res_df %>%
  dplyr::select(gene_id, SYMBOL, ENTREZID, everything())

write_xlsx(as.data.frame(res_df), path = "/Users/dwyun/rsem/my_data/DESeq2 sub/S1 vs S4.xlsx")

# 필터링 값에 따라 필터링 된 gene list 뽑을거면.

lfc_cut <- 1
fdr_cut <- 0.05

sig_df <- res_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    abs(log2FoldChange) >= lfc_cut
  )
write_xlsx(sig_df,
           path = "/Users/dwyun/rsem/my_data/DESeq2 sub/S1 vs S4 sig_genes.xlsx")

# up gene이나 down gene만 뽑고 싶을 때
sig_up <- res_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    log2FoldChange >= lfc_cut
  )

sig_down <- res_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    log2FoldChange <= -lfc_cut
  )


----------------------------------------------------------------------------------------------------------
 (4) S2 vs S3

setwd("/Users/dwyun/rsem/my_data")
rsem_dir <- "/Users/dwyun/rsem/my_data"
files <- list.files(rsem_dir, pattern = "\\.genes\\.results$", full.names = TRUE)
stopifnot(length(files) > 1)

samples <- sub("\\.genes\\.results$", "", basename(files))
names(files) <- samples

group <- ifelse(grepl("TdT_TW", samples), "TdT_TW", "prom1_TdT")
condition <- factor(group, levels = c("prom1_TdT","TdT_TW"))

coldata <- data.frame(sample = samples, condition = condition, row.names = samples)

condition
data.frame(samples, condition)

txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE,
                countsFromAbundance = "no")
txi$length[txi$length == 0] <- 1

dds <- DESeqDataSetFromTximport(txi = txi, colData = coldata, design = ~ condition)

dds <- DESeq(dds)

res <- results(dds,contrast = c("condition", "TdT_TW", "prom1_TdT"))
rownames(res) <- sub("\\.\\d+$", "", rownames(res))

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

map_orgdb <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(res_df$gene_id),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID")
) %>% 
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

res_df <- res_df %>%
  dplyr::left_join(map_orgdb, by = c("gene_id" = "ENSEMBL"))

res_df <- res_df %>%
  mutate(
    SYMBOL = coalesce(na_if(SYMBOL, ""), gene_id),
    ENTREZID = coalesce(na_if(as.character(ENTREZID), ""), gene_id)
  )

res_df <- res_df %>%
  dplyr::select(gene_id, SYMBOL, ENTREZID, everything())

write_xlsx(as.data.frame(res_df), path = "/Users/dwyun/rsem/my_data/DESeq2 sub/S2 vs S3.xlsx")

# 필터링 값에 따라 필터링 된 gene list 뽑을거면.

lfc_cut <- 1
fdr_cut <- 0.05

sig_df <- res_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    abs(log2FoldChange) >= lfc_cut
  )
write_xlsx(sig_df,
           path = "/Users/dwyun/rsem/my_data/DESeq2 sub/S2 vs S3 sig_genes.xlsx")

# up gene이나 down gene만 뽑고 싶을 때
sig_up <- res_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    log2FoldChange >= lfc_cut
  )

sig_down <- res_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    log2FoldChange <= -lfc_cut
  )

------------------------------------------------------------------------------------------------
(5) S2 vs S4

setwd("/Users/dwyun/rsem/my_data")
rsem_dir <- "/Users/dwyun/rsem/my_data"
files <- list.files(rsem_dir, pattern = "\\.genes\\.results$", full.names = TRUE)
stopifnot(length(files) > 1)

samples <- sub("\\.genes\\.results$", "", basename(files))
names(files) <- samples

group <- ifelse(grepl("TdT_TW", samples), "TdT_TW", "TdT")
condition <- factor(group, levels = c("TdT","TdT_TW"))

coldata <- data.frame(sample = samples, condition = condition, row.names = samples)

condition
data.frame(samples, condition)

txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE,
                countsFromAbundance = "no")
txi$length[txi$length == 0] <- 1

dds <- DESeqDataSetFromTximport(txi = txi, colData = coldata, design = ~ condition)

dds <- DESeq(dds)

res <- results(dds,contrast = c("condition", "TdT_TW", "TdT"))
rownames(res) <- sub("\\.\\d+$", "", rownames(res))

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

map_orgdb <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(res_df$gene_id),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID")
) %>% 
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

res_df <- res_df %>%
  dplyr::left_join(map_orgdb, by = c("gene_id" = "ENSEMBL"))

res_df <- res_df %>%
  mutate(
    SYMBOL = coalesce(na_if(SYMBOL, ""), gene_id),
    ENTREZID = coalesce(na_if(as.character(ENTREZID), ""), gene_id)
  )

res_df <- res_df %>%
  dplyr::select(gene_id, SYMBOL, ENTREZID, everything())

write_xlsx(as.data.frame(res_df), path = "/Users/dwyun/rsem/my_data/DESeq2 sub/S2 vs S4.xlsx")

# 필터링 값에 따라 필터링 된 gene list 뽑을거면.

lfc_cut <- 1
fdr_cut <- 0.05

sig_df <- res_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    abs(log2FoldChange) >= lfc_cut
  )
write_xlsx(sig_df,
           path = "/Users/dwyun/rsem/my_data/DESeq2 sub/S2 vs S4 sig_genes.xlsx")

# up gene이나 down gene만 뽑고 싶을 때
sig_up <- res_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    log2FoldChange >= lfc_cut
  )

sig_down <- res_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    log2FoldChange <= -lfc_cut
  )

-------------------------------------------------------------------------------------------
(6) S3 vs S4

setwd("/Users/dwyun/rsem/my_data")
rsem_dir <- "/Users/dwyun/rsem/my_data"
files <- list.files(rsem_dir, pattern = "\\.genes\\.results$", full.names = TRUE)
stopifnot(length(files) > 1)

samples <- sub("\\.genes\\.results$", "", basename(files))
names(files) <- samples

group <- ifelse(grepl("prom1_TdT", samples), "prom1_TdT", "TdT")
condition <- factor(group, levels = c("TdT","prom1_TdT"))

coldata <- data.frame(sample = samples, condition = condition, row.names = samples)

condition
data.frame(samples, condition)

txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE,
                countsFromAbundance = "no")
txi$length[txi$length == 0] <- 1

dds <- DESeqDataSetFromTximport(txi = txi, colData = coldata, design = ~ condition)

dds <- DESeq(dds)

res <- results(dds,contrast = c("condition", "prom1_TdT", "TdT"))
rownames(res) <- sub("\\.\\d+$", "", rownames(res))

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

map_orgdb <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(res_df$gene_id),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID")
) %>% 
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

res_df <- res_df %>%
  dplyr::left_join(map_orgdb, by = c("gene_id" = "ENSEMBL"))

res_df <- res_df %>%
  mutate(
    SYMBOL = coalesce(na_if(SYMBOL, ""), gene_id),
    ENTREZID = coalesce(na_if(as.character(ENTREZID), ""), gene_id)
  )

res_df <- res_df %>%
  dplyr::select(gene_id, SYMBOL, ENTREZID, everything())

write_xlsx(as.data.frame(res_df), path = "/Users/dwyun/rsem/my_data/DESeq2 sub/S3 vs S4.xlsx")

# 필터링 값에 따라 필터링 된 gene list 뽑을거면.

lfc_cut <- 1
fdr_cut <- 0.05

sig_df <- res_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    abs(log2FoldChange) >= lfc_cut
  )
write_xlsx(sig_df,
           path = "/Users/dwyun/rsem/my_data/DESeq2 sub/S3 vs S4 sig_genes.xlsx")

# up gene이나 down gene만 뽑고 싶을 때
sig_up <- res_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    log2FoldChange >= lfc_cut
  )

sig_down <- res_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    log2FoldChange <= -lfc_cut
  )

# MA Plot
plotMA(res, alpha = 0.01, ylim=c(-6,6))
