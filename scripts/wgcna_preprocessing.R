install.packages("WGCNA")

library(DESeq2)
library(tximport)
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()


rsem_dir <- "/Users/dwyun/rsem/data"   # ← rsem .genes.results 파일들이 있는 폴더

# .genes.results 파일 목록
files <- list.files(rsem_dir,
                    pattern = "\\.genes\\.results$",
                    full.names = TRUE)
stopifnot(length(files) > 1)

# sample 이름 = 파일명에서 .genes.results 제거한 것
samples <- sub("\\.genes\\.results$", "", basename(files))
names(files) <- samples

samples

condition <- dplyr::case_when(
  grepl("prom1_TdT_TW", samples) ~ "prom1_TdT_TW",
  grepl("prom1_TdT",     samples) ~ "prom1_TdT",
  grepl("TdT_TW",        samples) ~ "TdT_TW",
  grepl("TdT",           samples) ~ "TdT",
  TRUE ~ "unknown"
)

condition <- factor(
  condition,
  levels = c("prom1_TdT_TW", "TdT_TW", "prom1_TdT", "TdT")  # 순서는 원하는대로 조정
)

coldata <- data.frame(
  sample    = samples,
  condition = condition,
  row.names = samples
)

coldata

txi <- tximport(
  files,
  type = "rsem",
  txIn = FALSE,
  txOut = FALSE,
  countsFromAbundance = "no"
)

# 길이가 0인 유전자는 에러 방지용으로 1로 교체
txi$length[txi$length == 0] <- 1

dds <- DESeqDataSetFromTximport(
  txi     = txi,
  colData = coldata,
  design  = ~ condition
)

dds <- DESeq(dds)

# vst 변환 (WGCNA용 expression 값)
vsd <- vst(dds, blind = FALSE)
expr_mat <- assay(vsd)   # genes x samples
