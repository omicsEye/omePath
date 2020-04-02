# this script inlcudes scoring functions supported by deepath

test2groups <- function(data,
         metadata,
         meta = 'Group',
         case_label = 'case',
         control_label = 'control',
         test_type = 'wilcox.test',
         paired = F) {
  case_samples <- rownames(metadata[metadata[meta] == case_label, ])
  control_samples <- rownames(metadata[metadata[meta] == control_label, ])
  #data <- as.data.frame(t(data))
  case <- data[case_samples,]
  control <- data[control_samples,]
  #tranpose data to row as samples and columns as features
  
  stats_table <-
    setNames(
      data.frame(matrix(
        ncol = 5, nrow = dim(case)[2]
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
    if (all(is.na(case[,i])) || all(is.na(control[,i]))) {
      #print(i)
      next
    }
    #i <- 1
    stats_table[i, 'logFC'] <-
      log2(mean(case[,i], na.rm = TRUE)) - log2(mean(control[,i], na.rm = TRUE))
    tryCatch({
      if (test_type == 't.test') {
        temp_result <- t.test(case[,i], control[,i], paired = F)
        stats_table[i, 'P.Value'] <- temp_result$p.value
        stats_table[i, 'statistic'] <- temp_result$statistic
      } else{
        temp_result <- wilcox.test(case[,i], control[,i], paired = F)
        stats_table[i, 'P.Value'] <- temp_result$p.value
        stats_table[i, 'statistic'] <- temp_result$statistic
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
  return
  (stats_table)
}
