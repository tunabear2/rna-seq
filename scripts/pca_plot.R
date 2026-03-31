library(DESeq2)
library(tximport)
library(readr)
library(dplyr)
library(tibble)

# 1) genes.results 파일 목록
rsem_dir <- "/Users/dwyun/rsem/data/"   # 너가 쓰는 폴더로 맞춰줘
files <- list.files(rsem_dir,
                    pattern = "\\.genes\\.results$",
                    full.names = TRUE)
stopifnot(length(files) > 1)

# 2) 샘플 이름 (확장자 제거)
samples <- sub("\\.genes\\.results$", "", basename(files))
names(files) <- samples
samples

condition <- case_when(
  grepl("prom1_TdT_TW", samples) ~ "S1",
  grepl("prom1_TdT",     samples) & !grepl("TW", samples) ~ "S3",
  grepl("TdT_TW",        samples) ~ "S2",
  grepl("TdT",           samples) & !grepl("TW", samples) ~ "S4",
  TRUE ~ NA_character_
)

stopifnot(!any(is.na(condition)))  # 매핑 안 된 샘플 없는지 체크

condition <- factor(condition, levels = c("S1","S2","S3","S4"))

coldata_all <- data.frame(
  sample = samples,
  condition = condition,
  row.names = samples
)

coldata_all
table(coldata_all$condition)  # S1~S4에 샘플이 잘 들어갔는지 확인

# tximport로 RSEM genes.results 읽기
txi_all <- tximport(files,
                    type = "rsem",
                    txIn = FALSE,
                    txOut = FALSE,
                    countsFromAbundance = "no")

# 길이가 0인 gene은 1로 바꿔서 에러 막기
txi_all$length[txi_all$length == 0] <- 1

# 4개 condition(S1~S4)을 포함하는 DESeqDataSet
dds_all <- DESeqDataSetFromTximport(
  txi = txi_all,
  colData = coldata_all,
  design = ~ condition
)

dds_all
levels(dds_all$condition)   # "S1" "S2" "S3" "S4" 나오는지 확인

vsd <- vst(dds_all, blind = FALSE)
plotPCA(vsd, intgroup = "condition")

library(ggplot2)
library(ggrepel)

pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE, ntop = 1000)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2,
                    color = condition,        # S1,S2,S3,S4 색으로 구분
                    label = name)) +          # 샘플 이름
  geom_point(size = 3) +
  geom_text_repel(size = 3) +                 # 점 옆에 이름 표시
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()
