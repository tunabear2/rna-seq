# 1개 히스토그램 그리기
# expr_mat까지 되어 있는 상태 기준
t <- rowVars(expr_mat)
t <- t[is.finite(t) & tt > 0] # variance 0 제거

ggplot(data.frame(var = t), aes(x = var)) +
  geom_histogram(bins = 125) +
  scale_x_log10() +
  theme_bw() +
  labs(
    title = "Variance distribution (expr_mat_sub)",
    x = "gene-wise variance (log10 scale)",
    y = "Gene count"
  )
# sub set 만들어서 histogram 만들기
expr_mat_sub <- expr_mat[rownames(expr_mat) %in% gene_union3, , drop = FALSE]
v <- rowVars(expr_mat_sub)
v <- v[is.finite(v) & v > 0] # variance 0 제거
# 두개 합쳐서 하나의 histogram으로 표현하기

df <- dplyr::bind_rows(
  data.frame(var = t, group = "All genes"),
  data.frame(var = v, group = "DEG_union")
)

ggplot(df, aes(x = var, fill = group)) +
  geom_histogram(bins = 125, position = "identity", alpha = 0.45) +
  scale_x_log10() +
  theme_bw() +
  labs(
    title = "Variance distribution: All genes vs DEG_union",
    x = "gene-wise variance (log10 scale)",
    y = "Gene count",
    fill = NULL
  )

# 누적 히스토그램 그리고 cutoff도 그리기

var_vec <- t  # 또는 v

# 안전 필터
var_vec <- var_vec[is.finite(var_vec) & var_vec > 0]

# reverse cumulative 계산: 각 variance 값에서 ">= 그 값"인 gene 수
x <- sort(var_vec)                 # 오름차순
y <- rev(seq_along(x))             # N, N-1, ..., 1  (x[i] 이상 개수)

df_rc <- data.frame(var = x, n_ge = y)

ggplot(df_rc, aes(x = var, y = n_ge)) +
  geom_line() +
  scale_x_log10() +
  theme_bw() +
  labs(
    title = "Reverse cumulative curve of gene-wise variance",
    x = "Variance (log10 scale)",
    y = "Number of genes (>= variance)"
  )

cutoffs <- c(0.1, 0.15, 0.2, 0.5)

data.frame(
  cutoff = cutoffs,
  n_genes_ge = sapply(cutoffs, function(c) sum(var_vec >= c))
)
