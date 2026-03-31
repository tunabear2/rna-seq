-=-------------------------------------------------------------------------------
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
## 1) DESeq 실행
dds_all <- DESeq(dds_all)

## 2) variance 계산 (vst 기반)
vsd <- vst(dds_all, blind = FALSE)
vmat <- assay(vsd)
gene_var <- apply(vmat, 1, var)

get_top_var_genes <- function(N) {
  names(sort(gene_var, decreasing = TRUE))[seq_len(N)]
}

var_sizes <- c(162, 1647, 2928, 5712)
var_gene_lists <- lapply(var_sizes, get_top_var_genes)
names(var_gene_lists) <- paste0("topVar_", var_sizes)

## 3) pair 정의
pairs <- list(
  c("S1","S2"),
  c("S1","S3"),
  c("S1","S4"),
  c("S2","S3"),
  c("S2","S4"),
  c("S3","S4")
)

## 4) DEG 기준
alpha_cut <- 0.05
lfc_cut   <- 1

## 5) 모든 plot에서 공통 y축 범위 계산 (outlier 조금 무시)
get_global_ylim <- function(dds, pairs, q = 0.995) {
  vals <- c()
  for (p in pairs) {
    res <- results(dds, contrast = c("condition", p[1], p[2]))
    vals <- c(vals, res$log2FoldChange)
  }
  m <- as.numeric(quantile(abs(vals), q, na.rm = TRUE))
  m <- ceiling(m)
  c(-m, m)
}
ylim_global <- get_global_ylim(dds_all, pairs, q = 0.995)

## =========================
## [B] MA용 데이터 만들기
## =========================
make_ma_df <- function(dds, A, B, var_genes, alpha_cut, lfc_cut) {
  
  res <- results(dds, contrast = c("condition", A, B), alpha = alpha_cut)
  
  df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    mutate(
      # padj 유의 여부(=p-value 기준 표시용)
      is_padj_sig = !is.na(padj) & padj < alpha_cut,
      
      # DEG 정의: padj + |log2FC|
      is_deg = is_padj_sig & abs(log2FoldChange) >= lfc_cut,
      
      # variance 상위 여부
      is_var = gene_id %in% var_genes,
      
      # 4분류
      group = case_when(
        is_deg &  is_var ~ "overlap",
        is_deg & !is_var ~ "deg_only",
        !is_deg & is_var ~ "wgcna_only",
        TRUE ~ "other"
      ),
      
      x = log10(baseMean + 1)
    ) %>%
    filter(!is.na(log2FoldChange), !is.na(baseMean))
  
  # facet strip에 (개수) 붙이기
  cnt <- df %>% count(group, name = "n")
  df <- df %>%
    left_join(cnt, by = "group") %>%
    mutate(
      group_label = factor(
        paste0(group, "(", n, ")"),
        levels = paste0(c("deg_only","other","overlap","wgcna_only"),
                        "(", cnt$n[match(c("deg_only","other","overlap","wgcna_only"), cnt$group)], ")")
      )
    )
  
  df
}

## =========================
## [C] Plot 함수
##  - other: 검정(연하게)
##  - deg_only/overlap/wgcna_only: 색
##  - wgcna_only는 padj 유의/비유의로 진하기 차이
##  - y축 범위 통일 + |log2FC| 컷오프선(±1)
## =========================
plot_ma_facet <- function(df, title_text = "", ylim = c(-10,10), lfc_cut = 1) {
  
  ggplot(df, aes(x = x, y = log2FoldChange)) +
    
    # other: 검정 배경
    geom_point(
      data = df %>% filter(group == "other"),
      color = "black", alpha = 0.20, size = 0.5
    ) +
    
    # deg_only: 빨강
    geom_point(
      data = df %>% filter(group == "deg_only"),
      color = "#E64B35FF", alpha = 0.85, size = 0.8
    ) +
    
    # overlap: 청록
    geom_point(
      data = df %>% filter(group == "overlap"),
      color = "#4DBBD5FF", alpha = 0.85, size = 0.8
    ) +
    
    # wgcna_only 중 padj 유의: 보라(진하게)
    geom_point(
      data = df %>% filter(group == "wgcna_only", is_padj_sig),
      color = "#7E57C2", alpha = 0.90, size = 0.8
    ) +
    
    # wgcna_only 중 padj 비유의: 보라(연하게)
    geom_point(
      data = df %>% filter(group == "wgcna_only", !is_padj_sig),
      color = "#7E57C2", alpha = 0.20, size = 0.6
    ) +
    
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = c(-lfc_cut, lfc_cut), linetype = "dotted") +
    
    facet_wrap(~ group_label, ncol = 2) +
    coord_cartesian(ylim = ylim) +
    labs(
      title = title_text,
      x = "log10(baseMean + 1)",
      y = "log2 fold change"
    ) +
    theme_bw()
}

## =========================
## [D] 전체 실행: PDF로 저장
## =========================
pdf("MA_facets_DEG_vs_topVar_colored_fixedY.pdf", width = 12, height = 7)

for (vn in names(var_gene_lists)) {
  
  var_genes <- var_gene_lists[[vn]]
  
  for (p in pairs) {
    A <- p[1]; B <- p[2]
    
    df <- make_ma_df(dds_all, A, B, var_genes, alpha_cut, lfc_cut)
    
    g <- plot_ma_facet(
      df,
      title_text = paste0(A, " vs ", B, " | ", vn,
                          " | DEG: padj<", alpha_cut, ", |log2FC|>=", lfc_cut),
      ylim = ylim_global,
      lfc_cut = lfc_cut
    )
    
    print(g)
  }
}

dev.off()

install.packages("ggvenn")   # 한번만
library(ggvenn)

# 예: wgcna_genes <- topVar_5712 같은 리스트
#     deg_genes  <- ... (DEG 결과 gene_id 벡터)

install.packages("ggvenn")  # 처음 1번만
library(ggvenn)

make_venn_plot <- function(df, title_text = "") {
  
  # ✅ 각 MA plot(df)에서 실제로 사용된 gene set 추출
  wgcna_set <- unique(df$gene_id[df$is_var])  # topVar에 해당하는 gene들
  deg_set   <- unique(df$gene_id[df$is_deg])  # padj+|log2FC| 만족하는 DEG들
  
  venn_list <- list(
    WGCNA = wgcna_set,
    DEG   = deg_set
  )
  
  ggvenn(
    venn_list,
    show_elements = FALSE,
    show_percentage = TRUE,
    digits = 1,
    set_name_size = 6,
    text_size = 5,
    stroke_size = 1
  ) + ggtitle(title_text) + theme_void()
}
pdf("Venn_WGCNA_vs_DEG_each_MA.pdf", width = 6, height = 5)

for (vn in names(var_gene_lists)) {
  
  var_genes <- var_gene_lists[[vn]]
  
  for (p in pairs) {
    A <- p[1]; B <- p[2]
    
    # MA plot 만들 때 쓰던 df 그대로 재사용
    df <- make_ma_df(dds_all, A, B, var_genes, alpha_cut, lfc_cut)
    
    vp <- make_venn_plot(
      df,
      title_text = paste0(A, " vs ", B, " | ", vn, " | DEG: padj<", alpha_cut, ", |log2FC|>=", lfc_cut)
    )
    
    print(vp)
  }
}

dev.off()

    
    print(g)
  }
}

dev.off()
