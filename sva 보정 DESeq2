rsem_dir <- "/Users/dwyun/rsem/my_data/"   # 네 경로로
files <- list.files(rsem_dir, pattern = "\\.genes\\.results$", full.names = TRUE)

samples <- sub("\\.genes\\.results$", "", basename(files))
names(files) <- samples

# condition 정의 (네 파일명 패턴에 맞게)
condition <- case_when(
  grepl("^TdT_TW", samples) ~ "TdT_TW",
  grepl("^TdT_", samples) & !grepl("TW", samples) ~ "TdT",
  grepl("^prom1_TdT_TW", samples) ~ "prom1_TdT_TW",
  grepl("^prom1_TdT_", samples) & !grepl("TW", samples) ~ "prom1_TdT",
  TRUE ~ NA_character_
)

condition <- factor(condition,
                    levels = c("TdT", "TdT_TW", "prom1_TdT", "prom1_TdT_TW"))

coldata <- data.frame(sample = samples,
                      condition = condition,
                      row.names = samples)

# tximport
txi <- tximport(files, type = "rsem",
                txIn = FALSE, txOut = FALSE,
                countsFromAbundance = "no")

# sva용 count 행렬
counts <- round(txi$counts)

# low count 필터 (선택)
keep <- rowSums(counts) >= 1
counts <- counts[keep, ]

# coldata 순서 맞추기 (중요)
coldata <- coldata[colnames(counts), , drop = FALSE]

# full model: 관심 있는 조건 포함
mod  <- model.matrix(~ condition, data = coldata)
# null model: 조건 없이
mod0 <- model.matrix(~ 1, data = coldata)

svobj <- svaseq(as.matrix(counts), mod, mod0)

svobj$n.sv      # 몇 개 나왔는지 확인
head(svobj$sv)

for (i in seq_len(svobj$n.sv)) {
  coldata[[paste0("SV", i)]] <- svobj$sv[, i]
}

coldata

# 디자인식: ~ SV1 + SV2 + ... + condition
sv_terms <- paste0("SV", seq_len(svobj$n.sv), collapse = " + ")
design_formula <- as.formula(paste("~", sv_terms, "+ condition"))
design_formula
# 예: ~ SV1 + SV2 + condition

dds_sva <- DESeqDataSetFromMatrix(countData = counts,
                                  colData  = coldata,
                                  design   = design_formula)

dds_sva <- DESeq(dds_sva)

S2_vs_S3 <- results(dds_sva,
                    contrast = c("condition", "TdT_TW", "prom1_TdT"))
# TdT_TW vs TdT
S2_vs_S4 <- results(dds_sva,
                    contrast = c("condition", "TdT_TW", "TdT"))

# prom1_TdT vs TdT
S3_vs_S4 <- results(dds_sva,
                    contrast = c("condition", "prom1_TdT", "TdT"))
# prom1_TdT_TW vs TdT 
S1_vs_S4 <- results(dds_sva,
                    contrast = c("condition", "prom1_TdT_TW", "TdT"))
# prom1_TdT_TW vs prom1_TdT
S1_vs_S3 <- results(dds_sva,
                    contrast = c("condition", "prom1_TdT_TW", "prom1_TdT"))
# prom1_TdT_TW vs TdT_TW
S1_vs_S2 <- results(dds_sva,
                    contrast = c("condition", "prom1_TdT_TW", "TdT_TW"))

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

write_xlsx(as.data.frame(S1_vs_S2_df), path = "/Users/dwyun/rsem/my_data/DESeq2 sva/S1 vs S2.xlsx")

lfc_cut <- 1
fdr_cut <- 0.05

S1_vs_S2_sig_df <- S1_vs_S2_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    abs(log2FoldChange) >= lfc_cut
  )
write_xlsx(S1_vs_S2_sig_df,
           path = "/Users/dwyun/rsem/my_data/DEseq2 sva/S1 vs S2 sig_genes.xlsx")

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

write_xlsx(as.data.frame(S1_vs_S3_df), path = "/Users/dwyun/rsem/my_data/DESeq2 sva/S1 vs S3.xlsx")

lfc_cut <- 1
fdr_cut <- 0.05

S1_vs_S3_sig_df <- S1_vs_S3_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    abs(log2FoldChange) >= lfc_cut
  )
write_xlsx(S1_vs_S3_sig_df,
           path = "/Users/dwyun/rsem/my_data/DEseq2 sva/S1 vs S3 sig_genes.xlsx")
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

write_xlsx(as.data.frame(S1_vs_S4_df), path = "/Users/dwyun/rsem/my_data/DESeq2 sva/S1 vs S4.xlsx")

lfc_cut <- 1
fdr_cut <- 0.05

S1_vs_S4_sig_df <- S1_vs_S4_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    abs(log2FoldChange) >= lfc_cut
  )
write_xlsx(S1_vs_S4_sig_df,
           path = "/Users/dwyun/rsem/my_data/DEseq2 sva/S1 vs S4 sig_genes.xlsx")
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

write_xlsx(as.data.frame(S2_vs_S3_df), path = "/Users/dwyun/rsem/my_data/DESeq2 sva/S2 vs S3.xlsx")

lfc_cut <- 1
fdr_cut <- 0.05

S2_vs_S3_sig_df <- S2_vs_S3_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    abs(log2FoldChange) >= lfc_cut
  )
write_xlsx(S2_vs_S3_sig_df,
           path = "/Users/dwyun/rsem/my_data/DEseq2 sva/S2 vs S3 sig_genes.xlsx")
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

write_xlsx(as.data.frame(S2_vs_S4_df), path = "/Users/dwyun/rsem/my_data/DESeq2 sva/S2 vs S4.xlsx")

lfc_cut <- 1
fdr_cut <- 0.05

S2_vs_S4_sig_df <- S2_vs_S4_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    abs(log2FoldChange) >= lfc_cut
  )
write_xlsx(S2_vs_S4_sig_df,
           path = "/Users/dwyun/rsem/my_data/DEseq2 sva/S2 vs S4 sig_genes.xlsx")

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

write_xlsx(as.data.frame(S3_vs_S4_df), path = "/Users/dwyun/rsem/my_data/DESeq2 sva/S3 vs S4.xlsx")

lfc_cut <- 1
fdr_cut <- 0.05

S3_vs_S4_sig_df <- S3_vs_S4_df %>%
  filter(
    !is.na(padj),
    padj < fdr_cut,
    abs(log2FoldChange) >= lfc_cut
  )
write_xlsx(S3_vs_S4_sig_df,
           path = "/Users/dwyun/rsem/my_data/DEseq2 sva/S3 vs S4 sig_genes.xlsx")

-------------------------------
  공통 진 밴다이어그램 및 만들기
nstall.packages("ggvenn")  # 한 번만 설치
library(ggvenn)

genes_raw <- unique(prom1_TdT_TW_vs_TdT_TW_sig_genes$SYMBOL)     # 또는 gene_id
genes_sva <- unique(sva_prom1_TdT_TW_vs_TdT_TW_sig_genes$SYMBOL)

venn_list <- list(
  no_sva = genes_raw,
  sva    = genes_sva
)

ggvenn(venn_list,
       show_percentage = TRUE,
       show_elements = FALSE)  # TRUE로 하면 겹치는 gene 이름 다 나옴

# 공통 DEG
common_genes <- intersect(genes_raw, genes_sva)

# sva 안 쓴 분석에서만 DEG인 gene
only_raw <- setdiff(genes_raw, genes_sva)

# sva 쓴 분석에서만 DEG인 gene
only_sva <- setdiff(genes_sva, genes_raw)

length(common_genes); length(only_raw); length(only_sva)

common_df <- sva_prom1_TdT_TW_vs_TdT_TW_sig_genes %>% filter(SYMBOL %in% common_genes)  # 예: sva 쪽 정보 기준
# write.csv(common_df, "common_DEGs_sva_vs_nosva.csv", row.names = FALSE)
write.csv(common_df, "common_sva vs non_sva S1 vs S2.csv", row.names = FALSE)
