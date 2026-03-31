prom1_TdT_TW_1_genes_results$gene_id <- sub("\\.\\d+$", "", prom1_TdT_TW_1_genes_results$gene_id)
write_xlsx(prom1_TdT_TW_1_genes_results, "/Users/dwyun/prom1_TdT_TW_1_genes_results.xlsx")

metric <- "TPM"
name1 <- "my_sample"
name2 <- "rokit_sample"

collapse_fun <- mean
df1 <- my_prom1_TdT_TW_1 %>% group_by(gene_id) %>% summarise(tpm1 = collapse_fun(.data[[metric]]), .groups = "drop")
df2 <- rokitg_prom1_TdT_TW_1 %>% group_by(gene_id) %>% summarise(tpm2 = collapse_fun(.data[[metric]]), .groups = "drop")

cmd <- inner_join(df1, df2, by = "gene_id") %>%
  filter(is.finite(tpm1), is.finite(tpm2))

# 값의 스케일이 크면 로그가 보기 좋습니다.
cmd <- cmd %>% mutate(x = log10(tpm1 + 1), y = log10(tpm2 + 1))
cor_pearson  <- cor(cmd$x, cmd$y, method = "pearson")
cor_spearman <- cor(cmd$x, cmd$y, method = "spearman")
cor_pearson; cor_spearman


p <- ggplot(cmd, aes(x = x, y = y)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = paste0(name1, " ", metric, " (log10+1)"),
       y = paste0(name2, " ", metric, " (log10+1)"),
       title = paste0(metric, " scatter: ", name1, " vs ", name2),
       subtitle = sprintf("Pearson r=%.3f, Spearman ρ=%.3f", cor_pearson, cor_spearman)) +
  theme_bw()

print(p)
# 파일 저장
ggsave("scatter_TPM_Sample1_vs_Sample2.png", p, width = 6, height = 5, dpi = 300)

library(dplyr)
library(readr)

# ── 파라미터(필요시 조정) ─────────────────────────────
present_thr <- 0.1   # "있다"로 볼 TPM 기준 (>=)
absent_thr  <- 0.01   # "거의 없다"로 볼 TPM 기준 (<=)
mad_k       <- 3     # 잔차 MAD 배수 (2.5~4 권장 범위)

# ── A) 단독 발현(outlier-unique) ─────────────────────
only_name1 <- cmd %>%
  filter(tpm1 >= present_thr, tpm2 <= absent_thr) %>%
  transmute(gene_id,
            !!paste0(name1, "_TPM") := tpm1,
            !!paste0(name2, "_TPM") := tpm2,
            type = paste0("only_", name1))

only_name2 <- cmd %>%
  filter(tpm2 >= present_thr, tpm1 <= absent_thr) %>%
  transmute(gene_id,
            !!paste0(name1, "_TPM") := tpm1,
            !!paste0(name2, "_TPM") := tpm2,
            type = paste0("only_", name2))

only_unique <- bind_rows(only_name1, only_name2)

# ── B) 대각선에서 벗어난 통계적 아웃라이어 ───────────
cmd2 <- cmd %>% mutate(resid = y - x)  # y=x 기준 잔차
med_r <- median(cmd2$resid, na.rm = TRUE)
mad_r <- mad(cmd2$resid, constant = 1, na.rm = TRUE)
thr_hi <- med_r + mad_k * mad_r
thr_lo <- med_r - mad_k * mad_r

out_above <- cmd2 %>%
  filter(resid >= thr_hi) %>%
  transmute(gene_id,
            !!paste0(name1, "_TPM") := tpm1,
            !!paste0(name2, "_TPM") := tpm2,
            resid,
            type = paste0("above_diag (", name2, " >> ", name1, ")"))

out_below <- cmd2 %>%
  filter(resid <= thr_lo) %>%
  transmute(gene_id,
            !!paste0(name1, "_TPM") := tpm1,
            !!paste0(name2, "_TPM") := tpm2,
            resid,
            type = paste0("below_diag (", name1, " >> ", name2, ")"))

outliers_offdiag <- bind_rows(out_above, out_below)

# ── 정렬해서 확인(선택)

only_name1 %>% arrange(desc(!!sym(paste0(name1, "_TPM")))) %>% print(n = 164)
only_name2 %>% arrange(desc(!!sym(paste0(name2, "_TPM")))) %>% print(n = 64)
data.frame(only_name1)
write_xlsx(only_name1, path = "/Users/dwyun/only_my_sample.xlsx")
write_xlsx(only_name2, path = "/Users/dwyun/only_rokit_sample.xlsx")
out_above   %>% arrange(desc(resid)) %>% head()
out_below   %>% arrange(resid)       %>% head()
getwd()
