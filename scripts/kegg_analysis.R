1. prom1_TdT_TW vs prom1_TdT

outdir <- "rsem/results/kegg"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

bg_entrez <- res_df$ENTREZID %>% unique() %>% na.omit()

ent <- as.character(res_df$ENTREZID)
keep <- !is.na(ent) & grepl("^[0-9]+$", ent)
bg_entrez <- unique(ent[keep])

sig_entrez <- res_df %>%
  filter(!is.na(ENTREZID), !is.na(padj)) %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 0.585) %>%
  pull(ENTREZID) %>% unique()

length(bg_entrez); length(sig_entrez)

ekegg <- enrichKEGG(
  gene      = sig_entrez,
  universe  = bg_entrez,
  organism  = "hsa",
  keyType   = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  minGSSize     = 10
)

ekegg_r <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

if (!is.null(ekegg_r) && nrow(as.data.frame(ekegg_r)) > 0) {
  p_dot <- dotplot(ekegg_r, showCategory = 20) + ggtitle("KEGG ORA — dotplot")
  ggsave(file.path(outdir, "prom1_TdT_TW vs prom1_TdT kegg_ora_dotplot.png"), p_dot, width = 10, height = 8, dpi = 300)
  
  p_bar <- barplot(ekegg_r, showCategory = 20) + ggtitle("KEGG ORA — barplot")
  ggsave(file.path(outdir, "prom1_TdT_TW vs prom1_TdT kegg_ora_barplot.png"), p_bar, width = 10, height = 8, dpi = 300)
  
  # 네트워크류(선택)
  if (nrow(as.data.frame(ekegg_r)) >= 3) {
    p_cnet <- cnetplot(ekegg_r, showCategory = 10)
    ggsave(file.path(outdir, "prom1_TdT_TW vs prom1_TdT kegg_ora_cnetplot.png"), p_cnet, width = 10, height = 8, dpi = 300)
  }
  write.csv(as.data.frame(ekegg_r),
            file = file.path(outdir, "prom1_TdT_TW vs prom1_TdT kegg_ora_results.csv"),
            row.names = FALSE)
}
---------------------------------------------------------------------------------------------------------------------------------
2. prom1_TdT_TW vs TdT_TW

outdir <- "rsem/results/kegg"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

bg_entrez <- res_df$ENTREZID %>% unique() %>% na.omit()

ent <- as.character(res_df$ENTREZID)
keep <- !is.na(ent) & grepl("^[0-9]+$", ent)
bg_entrez <- unique(ent[keep])

sig_entrez <- res_df %>%
  filter(!is.na(ENTREZID), !is.na(padj)) %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 0.585) %>%
  pull(ENTREZID) %>% unique()

length(bg_entrez); length(sig_entrez)

ekegg <- enrichKEGG(
  gene      = sig_entrez,
  universe  = bg_entrez,
  organism  = "hsa",
  keyType   = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  minGSSize     = 10
)

ekegg_r <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

if (!is.null(ekegg_r) && nrow(as.data.frame(ekegg_r)) > 0) {
  p_dot <- dotplot(ekegg_r, showCategory = 20) + ggtitle("KEGG ORA — dotplot")
  ggsave(file.path(outdir, "prom1_TdT_TW vs TdT_TW kegg_ora_dotplot.png"), p_dot, width = 10, height = 8, dpi = 300)
  
  p_bar <- barplot(ekegg_r, showCategory = 20) + ggtitle("KEGG ORA — barplot")
  ggsave(file.path(outdir, "prom1_TdT_TW vs TdT_TW kegg_ora_barplot.png"), p_bar, width = 10, height = 8, dpi = 300)
  
  # 네트워크류(선택)
  if (nrow(as.data.frame(ekegg_r)) >= 3) {
    p_cnet <- cnetplot(ekegg_r, showCategory = 10)
    ggsave(file.path(outdir, "prom1_TdT_TW vs TdT_TW kegg_ora_cnetplot.png"), p_cnet, width = 10, height = 8, dpi = 300)
  }
  write.csv(as.data.frame(ekegg_r),
            file = file.path(outdir, "prom1_TdT_TW vs TdT_TW kegg_ora_results.csv"),
            row.names = FALSE)
}
-------------------------------------------------------------------------------------------------------------------------
3. TdT_TW vs TdT

outdir <- "rsem/results/kegg"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

bg_entrez <- res_df$ENTREZID %>% unique() %>% na.omit()

ent <- as.character(res_df$ENTREZID)
keep <- !is.na(ent) & grepl("^[0-9]+$", ent)
bg_entrez <- unique(ent[keep])

sig_entrez <- res_df %>%
  filter(!is.na(ENTREZID), !is.na(padj)) %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 0.585) %>%
  pull(ENTREZID) %>% unique()

length(bg_entrez); length(sig_entrez)

ekegg <- enrichKEGG(
  gene      = sig_entrez,
  universe  = bg_entrez,
  organism  = "hsa",
  keyType   = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  minGSSize     = 10
)

ekegg_r <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

if (!is.null(ekegg_r) && nrow(as.data.frame(ekegg_r)) > 0) {
  p_dot <- dotplot(ekegg_r, showCategory = 20) + ggtitle("KEGG ORA — dotplot")
  ggsave(file.path(outdir, "TdT_TW vs TdT kegg_ora_dotplot.png"), p_dot, width = 10, height = 8, dpi = 300)
  
  p_bar <- barplot(ekegg_r, showCategory = 20) + ggtitle("KEGG ORA — barplot")
  ggsave(file.path(outdir, "TdT_TW vs TdT kegg_ora_barplot.png"), p_bar, width = 10, height = 8, dpi = 300)
  
  # 네트워크류(선택)
  if (nrow(as.data.frame(ekegg_r)) >= 3) {
    p_cnet <- cnetplot(ekegg_r, showCategory = 10)
    ggsave(file.path(outdir, "TdT_TW vs TdT kegg_ora_cnetplot.png"), p_cnet, width = 10, height = 8, dpi = 300)
  }
  write.csv(as.data.frame(ekegg_r),
            file = file.path(outdir, "TdT_TW vs TdT kegg_ora_results.csv"),
            row.names = FALSE)
}
