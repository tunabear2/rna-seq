DESeq 이후 데이터로 진행

# clusterProfiler 패키지에 enrichGO 사용
library(clusterProfiler)
library(org.Hs.eg.db)
# gene에는 test할 gene의 entrezid 값을 넣어야 함.
GO_test <- enrichGO(
  gene = DEG_union,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  keyType = "ENTREZID",
  universe = unique(res_df$ENTREZID),
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  minGSSize = 5,
  maxGSSize = 500,
  readable = TRUE
)
GO_test2 <- as.data.frame(GO_test)

write_xlsx(GO_test2, "GO_test.xlsx")

library(enrichplot)
# showCategory에는 숫자를 입력해도 되고 보고싶은 terms만 입력해서 볼수도 있음
# 보고싶은 terms만 볼거면 list를 만들어서 하기, terms <- c("~", "~~", "~~~") 이렇게
dotplot(GO_test, showCategory = 10)
barplot(GO_test, showCategory = 10)
cnetplot(GO_test, showCategory = "extracellular matrix structural constituent")
heatplot(GO_test, showCategory = 10)
