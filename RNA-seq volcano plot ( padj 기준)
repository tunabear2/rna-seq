1. "prom1_TdT_TW", "prom1_TdT"

lfc_cut <- 1
fdr_cut <- 0.05

df_plot <- res_df %>%
  mutate(
    log2FC =as.numeric(log2FoldChange),
    padj = as.numeric(padj),
    pvalue = as.numeric(pvalue),
    sig = case_when(
      !is.na(padj) & padj < fdr_cut & log2FC >= lfc_cut ~ "Up",
      !is.na(padj) & padj < fdr_cut & log2FC <= -lfc_cut ~ "Down",
      TRUE ~ "NS"
    ),
    label = ifelse(is.na(SYMBOL) | SYMBOL=="", gene_id, SYMBOL)
  )

topN <- 15
label_top <- df_plot %>%
  filter(sig != "NS", !is.na(padj)) %>%
  arrange(padj, desc(abs(log2FC))) %>%
  slice_head(n = topN) %>%
  pull(SYMBOL)

df_plot <- df_plot %>%
  mutate(delabel = ifelse(SYMBOL %in% label_top, SYMBOL, NA_character_)
  )

volcanoplot <- ggplot(df_plot, aes(x = log2FC, y = -log10(padj),
                                   color = sig, label = delabel)) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), color = "gray", linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_cut), color = "gray", linetype = "dashed") +
  geom_point(size = 2, alpha = 0.9, stroke = 0) +
  scale_color_manual(values = c("Down" = "#79A9E1", "NS" = "grey70", "Up" = "#8B2B2B"),
                     labels = c("Down" = "prom1_TdT",
                                "NS" = "Not significant",
                                "Up" = "prom1_TdT_TW")) +
  labs(color = "Group",
       x = expression("log"[2]*"FC"),
       y = expression("-log"[10]*"FDR (padj)"),
       title = "Volcano plot: prom1_TdT vs prom1_TdT_TW") +
  coord_cartesian(ylim = c(0, NA), xlim = c(-10, 10)) +
  scale_x_continuous(breaks = seq(-10, 10, by = 2)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  ) +
  geom_text_repel(max.overlaps = Inf, na.rm = TRUE, box.padding = 0.4, point.padding = 0.2,
                  size = 3, segment.color = "grey60", min.segment.length = 0)

print(volcanoplot)

ttt <- data.frame(rank = seq_along(label_top), label_top = label_top)
write_xlsx(ttt, "C:/Users/rkawk/rsem/label_top(prom1_TdT vs prom1_TdT_TW).xlsx")

------------------------------------------------------------------------------------------------------------------------
2. "prom1_TdT_TW", "TdT_TW"


lfc_cut <- 1
fdr_cut <- 0.05

df_plot <- res_df %>%
  mutate(
    log2FC =as.numeric(log2FoldChange),
    padj = as.numeric(padj),
    pvalue = as.numeric(pvalue),
    sig = case_when(
      !is.na(padj) & padj < fdr_cut & log2FC >= lfc_cut ~ "Up",
      !is.na(padj) & padj < fdr_cut & log2FC <= -lfc_cut ~ "Down",
      TRUE ~ "NS"
    ),
    label = ifelse(is.na(SYMBOL) | SYMBOL=="", gene_id, SYMBOL)
  )

topN <- 15
label_top <- df_plot %>%
  filter(sig != "NS", !is.na(padj)) %>%
  arrange(padj, desc(abs(log2FC))) %>%
  slice_head(n = topN) %>%
  pull(SYMBOL)

df_plot <- df_plot %>%
  mutate(delabel = ifelse(SYMBOL %in% label_top, SYMBOL, NA_character_)
  )

volcanoplot <- ggplot(df_plot, aes(x = log2FC, y = -log10(padj),
                                   color = sig, label = delabel)) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), color = "gray", linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_cut), color = "gray", linetype = "dashed") +
  geom_point(size = 2, alpha = 0.9, stroke = 0) +
  scale_color_manual(values = c("Down" = "#79A9E1", "NS" = "grey70", "Up" = "#8B2B2B"),
                     labels = c("Down" = "TdT_TW",
                                "NS" = "Not significant",
                                "Up" = "prom1_TdT_TW")) +
  labs(color = "Group",
       x = expression("log"[2]*"FC"),
       y = expression("-log"[10]*"FDR (padj)"),
       title = "Volcano plot: TdT_TW vs prom1_TdT_TW") +
  coord_cartesian(ylim = c(0, NA), xlim = c(-10, 10)) +
  scale_x_continuous(breaks = seq(-10, 10, by = 2)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  ) +
  geom_text_repel(max.overlaps = Inf, na.rm = TRUE, box.padding = 0.4, point.padding = 0.2,
                  size = 3, segment.color = "grey60", min.segment.length = 0)

print(volcanoplot)

ttt <- data.frame(rank = seq_along(label_top), label_top = label_top)
write_xlsx(ttt, "C:/Users/rkawk/rsem/label_top(TdT_TW vs prom1_TdT_TW).xlsx")

------------------------------------------------------------------------------------------------------------------
3. TdT vs TdT_TW


lfc_cut <- 1
fdr_cut <- 0.05

df_plot <- res_df %>%
  mutate(
    log2FC =as.numeric(log2FoldChange),
    padj = as.numeric(padj),
    pvalue = as.numeric(pvalue),
    sig = case_when(
      !is.na(padj) & padj < fdr_cut & log2FC >= lfc_cut ~ "Up",
      !is.na(padj) & padj < fdr_cut & log2FC <= -lfc_cut ~ "Down",
      TRUE ~ "NS"
    ),
    label = ifelse(is.na(SYMBOL) | SYMBOL=="", gene_id, SYMBOL)
  )

topN <- 15
label_top <- df_plot %>%
  filter(sig != "NS", !is.na(padj)) %>%
  arrange(padj, desc(abs(log2FC))) %>%
  slice_head(n = topN) %>%
  pull(SYMBOL)

df_plot <- df_plot %>%
  mutate(delabel = ifelse(SYMBOL %in% label_top, SYMBOL, NA_character_)
  )

volcanoplot <- ggplot(df_plot, aes(x = log2FC, y = -log10(padj),
                                   color = sig, label = delabel)) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), color = "gray", linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_cut), color = "gray", linetype = "dashed") +
  geom_point(size = 2, alpha = 0.9, stroke = 0) +
  scale_color_manual(values = c("Down" = "#79A9E1", "NS" = "grey70", "Up" = "#8B2B2B"),
                     labels = c("Down" = "TdT_TW",
                                "NS" = "Not significant",
                                "Up" = "TdT")) +
  labs(color = "Group",
       x = expression("log"[2]*"FC"),
       y = expression("-log"[10]*"FDR (padj)"),
       title = "Volcano plot: TdT vs TdT_TW") +
  coord_cartesian(ylim = c(0, NA), xlim = c(-10, 10)) +
  scale_x_continuous(breaks = seq(-10, 10, by = 2)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  ) +
  geom_text_repel(max.overlaps = Inf, na.rm = TRUE, box.padding = 0.4, point.padding = 0.2,
                  size = 3, segment.color = "grey60", min.segment.length = 0)

print(volcanoplot)

ttt <- data.frame(rank = seq_along(label_top), label_top = label_top)
write_xlsx(ttt, "C:/Users/rkawk/rsem/label_top(TdT vs TdT_TW).xlsx")
------------------------------------------------------------------------------------------------------------
4. "prom1_TdT", "TdT"

lfc_cut <- 1
fdr_cut <- 0.05

df_plot <- res_df %>%
  mutate(
    log2FC =as.numeric(log2FoldChange),
    padj = as.numeric(padj),
    pvalue = as.numeric(pvalue),
    sig = case_when(
      !is.na(padj) & padj < fdr_cut & log2FC >= lfc_cut ~ "Up",
      !is.na(padj) & padj < fdr_cut & log2FC <= -lfc_cut ~ "Down",
      TRUE ~ "NS"
    ),
    label = ifelse(is.na(SYMBOL) | SYMBOL=="", gene_id, SYMBOL)
  )

topN <- 15
label_top <- df_plot %>%
  filter(sig != "NS", !is.na(padj)) %>%
  arrange(padj, desc(abs(log2FC))) %>%
  slice_head(n = topN) %>%
  pull(SYMBOL)

df_plot <- df_plot %>%
  mutate(delabel = ifelse(SYMBOL %in% label_top, SYMBOL, NA_character_)
  )

volcanoplot <- ggplot(df_plot, aes(x = log2FC, y = -log10(padj),
                                   color = sig, label = delabel)) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), color = "gray", linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_cut), color = "gray", linetype = "dashed") +
  geom_point(size = 2, alpha = 0.9, stroke = 0) +
  scale_color_manual(values = c("Down" = "#79A9E1", "NS" = "grey70", "Up" = "#8B2B2B"),
                     labels = c("Down" = "prom1_TdT",
                                "NS" = "Not significant",
                                "Up" = "TdT")) +
  labs(color = "Group",
       x = expression("log"[2]*"FC"),
       y = expression("-log"[10]*"FDR (padj)"),
       title = "Volcano plot: TdT vs prom1_TdT") +
  coord_cartesian(ylim = c(0, NA), xlim = c(-10, 10)) +
  scale_x_continuous(breaks = seq(-10, 10, by = 2)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  ) +
  geom_text_repel(max.overlaps = Inf, na.rm = TRUE, box.padding = 0.4, point.padding = 0.2,
                  size = 3, segment.color = "grey60", min.segment.length = 0)

print(volcanoplot)

ttt <- data.frame(rank = seq_along(label_top), label_top = label_top)
write_xlsx(ttt, "C:/Users/rkawk/rsem/label_top(TdT vs prom1_TdT).xlsx")
---------------------------------------------------------------------------------------------------
5. TdT vs prom1_TdT_TW

lfc_cut <- 1
fdr_cut <- 0.05

df_plot <- res_df %>%
  mutate(
    log2FC =as.numeric(log2FoldChange),
    padj = as.numeric(padj),
    pvalue = as.numeric(pvalue),
    sig = case_when(
      !is.na(padj) & padj < fdr_cut & log2FC >= lfc_cut ~ "Up",
      !is.na(padj) & padj < fdr_cut & log2FC <= -lfc_cut ~ "Down",
      TRUE ~ "NS"
    ),
    label = ifelse(is.na(SYMBOL) | SYMBOL=="", gene_id, SYMBOL)
  )

topN <- 20
label_top <- df_plot %>%
  filter(sig != "NS", !is.na(padj)) %>%
  arrange(padj, desc(abs(log2FC))) %>%
  slice_head(n = topN) %>%
  pull(SYMBOL)

df_plot <- df_plot %>%
  mutate(delabel = ifelse(SYMBOL %in% label_top, SYMBOL, NA_character_)
  )

volcanoplot <- ggplot(df_plot, aes(x = log2FC, y = -log10(padj),
                                   color = sig, label = delabel)) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), color = "gray", linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_cut), color = "gray", linetype = "dashed") +
  geom_point(size = 2, alpha = 0.9, stroke = 0) +
  scale_color_manual(values = c("Down" = "#79A9E1", "NS" = "grey70", "Up" = "#8B2B2B"),
                     labels = c("Down" = "TdT",
                                "NS" = "Not significant",
                                "Up" = "prom1_TdT_TW")) +
  labs(color = "Group",
       x = expression("log"[2]*"FC"),
       y = expression("-log"[10]*"FDR (padj)"),
       title = "Volcano plot: TdT vs prom1_TdT_TW") +
  coord_cartesian(ylim = c(0, NA), xlim = c(-10, 10)) +
  scale_x_continuous(breaks = seq(-10, 10, by = 2)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  ) +
  geom_text_repel(max.overlaps = Inf, na.rm = TRUE, box.padding = 0.4, point.padding = 0.2,
                  size = 3, segment.color = "grey60", min.segment.length = 0)

print(volcanoplot)

ttt <- data.frame(rank = seq_along(label_top), label_top = label_top)
write_xlsx(ttt, "C:/Users/rkawk/rsem/label_top(TdT vs prom1_TdT_TW).xlsx")
------------------------------------------------------------------------------------------------------------
(6) TdT_TW vs prom1_TdT

lfc_cut <- 1
fdr_cut <- 0.05

df_plot <- res_df %>%
  mutate(
    log2FC =as.numeric(log2FoldChange),
    padj = as.numeric(padj),
    pvalue = as.numeric(pvalue),
    sig = case_when(
      !is.na(padj) & padj < fdr_cut & log2FC >= lfc_cut ~ "Up",
      !is.na(padj) & padj < fdr_cut & log2FC <= -lfc_cut ~ "Down",
      TRUE ~ "NS"
    ),
    label = ifelse(is.na(SYMBOL) | SYMBOL=="", gene_id, SYMBOL)
  )

topN <- 20
label_top <- df_plot %>%
  filter(sig != "NS", !is.na(padj)) %>%
  arrange(padj, desc(abs(log2FC))) %>%
  slice_head(n = topN) %>%
  pull(SYMBOL)

df_plot <- df_plot %>%
  mutate(delabel = ifelse(SYMBOL %in% label_top, SYMBOL, NA_character_)
  )

volcanoplot <- ggplot(df_plot, aes(x = log2FC, y = -log10(padj),
                                   color = sig, label = delabel)) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), color = "gray", linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_cut), color = "gray", linetype = "dashed") +
  geom_point(size = 2, alpha = 0.9, stroke = 0) +
  scale_color_manual(values = c("Down" = "#79A9E1", "NS" = "grey70", "Up" = "#8B2B2B"),
                     labels = c("Down" = "TdT_TW",
                                "NS" = "Not significant",
                                "Up" = "prom1_TdT")) +
  labs(color = "Group",
       x = expression("log"[2]*"FC"),
       y = expression("-log"[10]*"FDR (padj)"),
       title = "Volcano plot: TdT_TW vs prom1_TdT") +
  coord_cartesian(ylim = c(0, NA), xlim = c(-10, 10)) +
  scale_x_continuous(breaks = seq(-10, 10, by = 2)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  ) +
  geom_text_repel(max.overlaps = Inf, na.rm = TRUE, box.padding = 0.4, point.padding = 0.2,
                  size = 3, segment.color = "grey60", min.segment.length = 0)

print(volcanoplot)

ttt <- data.frame(rank = seq_along(label_top), label_top = label_top)
write_xlsx(ttt, "C:/Users/rkawk/rsem/label_top(TdT_TW vs prom1_TdT).xlsx")
----------------------------------------------------------------------------------------------------------
5. log2FC 버전은 이 부분만 수정

label_top <- df_plot %>%
  filter(sig != "NS") %>%
  slice_max(abs(log2FC), n = topN, with_ties = FALSE) %>%
  pull(SYMBOL)

topN <- 15
label_top <- df_plot %>%
  filter(sig != "NS", !is.na(padj)) %>%
  arrange(padj, desc(abs(log2FC))) %>%
  slice_head(n = topN) %>%
  pull(SYMBOL)
