-------------------

setwd("/Users/dwyun/rsem/my_data")  # 예: setwd("/Users/tuna/rsem")

# 2) 패키지 로드
library(tximport)
library(DESeq2)
library(readr)
# 3) genes.results 파일 목록 가져오기
files <- list.files(pattern = "\\.genes\\.results$", full.names = TRUE)
files
length(files)   # 12개 나오는지 확인
# 4) 샘플 이름 만들기 (확장자만 제거)
samples <- sub("\\.genes\\.results$", "", basename(files))
samples
names(files) <- samples
# 5) condition 벡터 생성
condition <- dplyr::case_when(
  grepl("prom1_TdT_TW", samples) ~ "S1",
  grepl("prom1_TdT",    samples) ~ "S3",
  grepl("TdT_TW",       samples) ~ "S2",
  grepl("TdT_",         samples) ~ "S4",
  TRUE ~ NA_character_
)

condition

# 6) factor로 정리 (원하는 순서대로 level 지정)
condition <- factor(condition,
                    levels = c("S1", "S3", "S2", "S4"))

# 7) coldata 만들기
coldata <- data.frame(
  sample = samples,
  condition = condition,
  row.names = samples
)
coldata

# 8) tximport로 RSEM 결과 불러오기
txi <- tximport(files,
                type = "rsem",
                txIn  = FALSE,
                txOut = FALSE,
                countsFromAbundance = "no")

# length에 0이 있으면 DESeq2가 에러 내니까 1로 바꿔주기
txi$length[txi$length == 0] <- 1

# 9) DESeqDataSet 생성
dds_all <- DESeqDataSetFromTximport(txi = txi,
                                    colData = coldata,
                                    design  = ~ condition)

dds_all

# 10) DESeq 실행
dds_all <- DESeq(dds_all)

# TdT_TW vs prom1_TdT
S2_vs_S3 <- results(dds_all,
                    contrast = c("condition", "S2", "S3"))
# TdT_TW vs TdT
S2_vs_S4 <- results(dds_all,
                    contrast = c("condition", "S2", "S4"))

# prom1_TdT vs TdT
S3_vs_S4 <- results(dds_all,
                    contrast = c("condition", "S3", "S4"))
# prom1_TdT_TW vs TdT 
S1_vs_S4 <- results(dds_all,
                    contrast = c("condition", "S1", "S4"))
# prom1_TdT_TW vs prom1_TdT
S1_vs_S3 <- results(dds_all,
                    contrast = c("condition", "S1", "S3"))
# prom1_TdT_TW vs TdT_TW
S1_vs_S2 <- results(dds_all,
                    contrast = c("condition", "S1", "S2"))

# S1 vs S2
S1_vs_S2_df <- as.data.frame(S1_vs_S2)
rownames(S1_vs_S2_df) <- sub("\\.\\d+$", "", rownames(S1_vs_S2_df))

S1_vs_S2_df$gene_id <- rownames(S1_vs_S2_df)

map_orgdb <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(S1_vs_S2_df$gene_id),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID")
) %>% 
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

S1_vs_S2_df <- S1_vs_S2_df %>%
  dplyr::left_join(map_orgdb, by = c("gene_id" = "ENSEMBL"))

S1_vs_S2_df <- S1_vs_S2_df %>%
  mutate(
    SYMBOL = coalesce(na_if(SYMBOL, ""), gene_id),
    ENTREZID = coalesce(na_if(as.character(ENTREZID), ""), gene_id)
  )

S1_vs_S2_df <- S1_vs_S2_df %>%
  dplyr::select(gene_id, SYMBOL, ENTREZID, everything())

write_xlsx(as.data.frame(S1_vs_S2_df), path = "/Users/dwyun/rsem/my_data/DESeq2 data/S1 vs S2.xlsx")

lfc_cut <- 1
fdr_cut <- 0.05

S1_vs_S2_sig_df <- S1_vs_S2_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    abs(log2FoldChange) >= lfc_cut
  )
write_xlsx(S1_vs_S2_sig_df,
           path = "/Users/dwyun/rsem/my_data/DEseq2 data/S1 vs S2 sig_genes.xlsx")

# S1 vs S3

S1_vs_S3_df <- as.data.frame(S1_vs_S3)
rownames(S1_vs_S3_df) <- sub("\\.\\d+$", "", rownames(S1_vs_S3_df))

S1_vs_S3_df$gene_id <- rownames(S1_vs_S3_df)

map_orgdb <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(S1_vs_S3_df$gene_id),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID")
) %>% 
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

S1_vs_S3_df <- S1_vs_S3_df %>%
  dplyr::left_join(map_orgdb, by = c("gene_id" = "ENSEMBL"))

S1_vs_S3_df <- S1_vs_S3_df %>%
  mutate(
    SYMBOL = coalesce(na_if(SYMBOL, ""), gene_id),
    ENTREZID = coalesce(na_if(as.character(ENTREZID), ""), gene_id)
  )

S1_vs_S3_df <- S1_vs_S3_df %>%
  dplyr::select(gene_id, SYMBOL, ENTREZID, everything())

write_xlsx(as.data.frame(S1_vs_S3_df), path = "/Users/dwyun/rsem/my_data/DESeq2 data/S1 vs S3.xlsx")

lfc_cut <- 1
fdr_cut <- 0.05

S1_vs_S3_sig_df <- S1_vs_S3_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    abs(log2FoldChange) >= lfc_cut
  )
write_xlsx(S1_vs_S3_sig_df,
           path = "/Users/dwyun/rsem/my_data/DEseq2 data/S1 vs S3 sig_genes.xlsx")
# S1 vs S4

S1_vs_S4_df <- as.data.frame(S1_vs_S4)
rownames(S1_vs_S4_df) <- sub("\\.\\d+$", "", rownames(S1_vs_S4_df))

S1_vs_S4_df$gene_id <- rownames(S1_vs_S4_df)

map_orgdb <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(S1_vs_S4_df$gene_id),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID")
) %>% 
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

S1_vs_S4_df <- S1_vs_S4_df %>%
  dplyr::left_join(map_orgdb, by = c("gene_id" = "ENSEMBL"))

S1_vs_S4_df <- S1_vs_S4_df %>%
  mutate(
    SYMBOL = coalesce(na_if(SYMBOL, ""), gene_id),
    ENTREZID = coalesce(na_if(as.character(ENTREZID), ""), gene_id)
  )

S1_vs_S4_df <- S1_vs_S4_df %>%
  dplyr::select(gene_id, SYMBOL, ENTREZID, everything())

write_xlsx(as.data.frame(S1_vs_S4_df), path = "/Users/dwyun/rsem/my_data/DESeq2 data/S1 vs S4.xlsx")

lfc_cut <- 1
fdr_cut <- 0.05

S1_vs_S4_sig_df <- S1_vs_S4_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    abs(log2FoldChange) >= lfc_cut
  )
write_xlsx(S1_vs_S4_sig_df,
           path = "/Users/dwyun/rsem/my_data/DEseq2 data/S1 vs S4 sig_genes.xlsx")
# S2 vs S3

S2_vs_S3_df <- as.data.frame(S2_vs_S3)
rownames(S2_vs_S3_df) <- sub("\\.\\d+$", "", rownames(S2_vs_S3_df))

S2_vs_S3_df$gene_id <- rownames(S2_vs_S3_df)

map_orgdb <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(S2_vs_S3_df$gene_id),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID")
) %>% 
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

S2_vs_S3_df <- S2_vs_S3_df %>%
  dplyr::left_join(map_orgdb, by = c("gene_id" = "ENSEMBL"))

S2_vs_S3_df <- S2_vs_S3_df %>%
  mutate(
    SYMBOL = coalesce(na_if(SYMBOL, ""), gene_id),
    ENTREZID = coalesce(na_if(as.character(ENTREZID), ""), gene_id)
  )

S2_vs_S3_df <- S2_vs_S3_df %>%
  dplyr::select(gene_id, SYMBOL, ENTREZID, everything())

write_xlsx(as.data.frame(S2_vs_S3_df), path = "/Users/dwyun/rsem/my_data/DESeq2 data/S2 vs S3.xlsx")

lfc_cut <- 1
fdr_cut <- 0.05

S2_vs_S3_sig_df <- S2_vs_S3_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    abs(log2FoldChange) >= lfc_cut
  )
write_xlsx(S2_vs_S3_sig_df,
           path = "/Users/dwyun/rsem/my_data/DEseq2 data/S2 vs S3 sig_genes.xlsx")
# S2 vs S4

S2_vs_S4_df <- as.data.frame(S2_vs_S4)
rownames(S2_vs_S4_df) <- sub("\\.\\d+$", "", rownames(S2_vs_S4_df))

S2_vs_S4_df$gene_id <- rownames(S2_vs_S4_df)

map_orgdb <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(S2_vs_S4_df$gene_id),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID")
) %>% 
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

S2_vs_S4_df <- S2_vs_S4_df %>%
  dplyr::left_join(map_orgdb, by = c("gene_id" = "ENSEMBL"))

S2_vs_S4_df <- S2_vs_S4_df %>%
  mutate(
    SYMBOL = coalesce(na_if(SYMBOL, ""), gene_id),
    ENTREZID = coalesce(na_if(as.character(ENTREZID), ""), gene_id)
  )

S2_vs_S4_df <- S2_vs_S4_df %>%
  dplyr::select(gene_id, SYMBOL, ENTREZID, everything())

write_xlsx(as.data.frame(S2_vs_S4_df), path = "/Users/dwyun/rsem/my_data/DESeq2 data/S2 vs S4.xlsx")

lfc_cut <- 1
fdr_cut <- 0.05

S2_vs_S4_sig_df <- S2_vs_S4_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    abs(log2FoldChange) >= lfc_cut
  )
write_xlsx(S2_vs_S4_sig_df,
           path = "/Users/dwyun/rsem/my_data/DEseq2 data/S2 vs S4 sig_genes.xlsx")

# S3 vs S4

S3_vs_S4_df <- as.data.frame(S3_vs_S4)
rownames(S3_vs_S4_df) <- sub("\\.\\d+$", "", rownames(S3_vs_S4_df))

S3_vs_S4_df$gene_id <- rownames(S3_vs_S4_df)

map_orgdb <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(S3_vs_S4_df$gene_id),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID")
) %>% 
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

S3_vs_S4_df <- S3_vs_S4_df %>%
  dplyr::left_join(map_orgdb, by = c("gene_id" = "ENSEMBL"))

S3_vs_S4_df <- S3_vs_S4_df %>%
  mutate(
    SYMBOL = coalesce(na_if(SYMBOL, ""), gene_id),
    ENTREZID = coalesce(na_if(as.character(ENTREZID), ""), gene_id)
  )

S3_vs_S4_df <- S3_vs_S4_df %>%
  dplyr::select(gene_id, SYMBOL, ENTREZID, everything())

write_xlsx(as.data.frame(S3_vs_S4_df), path = "/Users/dwyun/rsem/my_data/DESeq2 data/S3 vs S4.xlsx")

lfc_cut <- 1
fdr_cut <- 0.05

S3_vs_S4_sig_df <- S3_vs_S4_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    abs(log2FoldChange) >= lfc_cut
  )
write_xlsx(S3_vs_S4_sig_df,
           path = "/Users/dwyun/rsem/my_data/DEseq2 data/S3 vs S4 sig_genes.xlsx")
