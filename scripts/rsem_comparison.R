infile  <- "C:/path/to/your/sample.genes.results"   # <- 파일 경로 바꿔줘
outfile <- "C:/path/to/your/sample.genes.results.xlsx"

# RSEM genes.results 읽기
df <- read_tsv(infile, col_types = cols())  # 자동 타입 감지

# 엑셀로 내보내기
write_xlsx(df, outfile)

cat("저장 완료:", outfile, "\n")
