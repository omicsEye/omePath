# this script inlcudes scoring functions supported by deepath

test2groups <- function(data,
         metadata,
         meta = 'Group',
         case_lable = 'case',
         control_label = 'control',
         test_type = 'wilcox.test',
         paired = F) {
  case <- rownames(metadata[metadata$Group == case_lable, ])
  control <- rownames(metadata[metadata$Group == control_label, ])
  Case <- data[case,]
  Control <- data[control,]
  stats_table <-
    setNames(
      data.frame(matrix(
        ncol = 5, nrow = dim(Case)[2]
      )),
      c(
        "logFC",
        "statistic",
        "P.Value",
        "adj.P.Val",
        "fdr"
      )
    )
  rownames(stats_table) <- colnames(case)
  for (i in 1:dim(stats_table)[1]) {
    if (all(is.na(case[, i])) || all(is.na(control[, i]))) {
      #print(i)
      next
    }
    stats_table[i, 'logFC'] <-
      log2(mean(case[, i], na.rm = TRUE)) - log2(mean(control[, i], na.rm = TRUE))
    tryCatch({
      if (test_type == 't.test') {
        stats_table[i, 'P.Value'] <-
          t.test(case[, i], control[, i], paired = paired)$p.value
        stats_table[i, 'statistic'] <-
          t.test(case[, i], control[, i], paired = paired)$statistic
      } else{
        stats_table[i, 'P.Value'] <-
          wilcox.test(case[, i], control[, i], paired = paired)$p.value
        stats_table[i, 'statistic'] <-
          wilcox.test(case[, i], control[, i], paired = paired)$statistic
      }
    }, error = function(e) {
      stats_table[i, 'P.Value'] <- NA
      stats_table[i, 'statistic'] <- NA
    })
  }
  
  stats_table$fdr <-
    p.adjust(
      as.vector(stats_table$P.Value),
      method = 'BH',
      n = length(stats_table$P.Value)
    )
  stats_table$feature <- rownames(stats_table)
  return (stats_table)
}
