# load in the required libraries, report an error if they are not installed
for (lib in c('gsEasy', 'future', 'ggplot2')) {
  if (!suppressPackageStartupMessages(require(lib, character.only = TRUE)))
    stop(paste("Please install the R package: ", lib))
}

plan(multisession)

OSEA <- function(stats_table,
                 output = "~/",
                 score_col = 'logFC',
                 pval_threshold = 0.05,
                 fdr_threshold = NA,
                 Pathway.Subject = NA,
                 method = 'gset',
                 min_member = 2,
                 mapper_file = NA,
                 do_plot = TRUE,
                 pathway_col = "Pathway",
                 feature_col = "Feature") {
  # load mapping files
  if (is.character(mapper_file)) {
    mapper_pathway2feature <-
      data.frame(data.table::fread(
        mapper_file,
        header = TRUE,
        check.names = FALSE,
        sep = "\t"
      ))
    if (nrow(mapper_pathway2feature) == 1) {
      # read again to get row name
      mapper_pathway2feature <-
        utils::read.table(
          mapper_file,
          header = TRUE,
          check.names = FALSE,
          row.names = 1
        )
    }
  } else {
    mapper_pathway2feature <- mapper_file
  }
  if (!is.na(Pathway.Subject)) {
    if ('Subject' %in% colnames(mapper_pathway2feature))
      mapper_pathway2feature <-
        mapper_pathway2feature[mapper_pathway2feature$Subject == Pathway.Subject, ]
    else
      print(paste('Pathway.Subject is not in mapping file!!!'))
  } else{
    Pathway.Subject <- ''
  }
  if (is.character(stats_table)) {
    stats <-
      data.frame(data.table::fread(stats_table, header = TRUE, sep = "\t"),
                 row.names = 1)
    if (nrow(stats_table) == 1) {
      # read again to get row name
      stats <-
        utils::read.table(
          stats_table,
          header = TRUE,
          row.names = 1,
          sep = "\t",
          fill = FALSE,
          comment.char = "" ,
          check.names = FALSE
        )
    }
  } else {
    stats <- stats_table
  }
  stats <- stats[!is.na(stats[score_col]), , drop = F]
  stats$score_rank <- NA
  if (score_col == 'P.Value') {
    stats <- stats |>
      dplyr::arrange(dplyr::desc(score_col))
  } else{
    stats <- stats |>
      dplyr::arrange(score_col)
  }
  stats$score_rank <- seq.int(nrow(stats))
  
  enrichment_stats <-
    data.frame(
      pathway = c(),
      pathway_members = c(),
      pval = c(),
      fdr = c(),
      n = c(),
      N = c(),
      set_enrichment_score = c()
    )
  
  # calculate p value for all sets (pathway terms)
  logging::loginfo(
    "Loading files is done, now calculating enrichments p-value! This may take 10 or more minutes!!!"
  )
  pathways <-
    unique(unlist(lapply(mapper_pathway2feature[pathway_col], as.character)))
  
  for (i in 1:length(pathways)) {
    current_member <- pathways[i] # "Purine Metabolism" #
    pathway_members <-
      mapper_pathway2feature[as.character(mapper_pathway2feature[, pathway_col]) == current_member, feature_col]
    pathway_members_in_study <-
      intersect(pathway_members, rownames(stats))
    N <- length(pathway_members)
    n <- length(pathway_members_in_study)
    if (n < min_member)
      next
    if (method == 'ks') {
      indx <- match(pathway_members_in_study, rownames(stats))
      if (length(indx) > 0) {
        stats_val <- stats::ks.test(indx, stats$score_rank)
        pval <- stats_val$p.value
      } else
        pval <- 1.0
      set_enrichment_score <-
        stats_val$statistic #gsEasy::gset(S = pathway_members_in_study, r = rownames(stats), raw_score = T)
    } else if (method == 'gset') {
      pval <-
        gsEasy::gset(S = pathway_members_in_study, r = rownames(stats))
      set_enrichment_score <-
        gsEasy::gset(S = pathway_members_in_study,
                     r = rownames(stats),
                     raw_score = T)
    } else if (method == 'wilcox') {
      indx <- match(pathway_members_in_study, rownames(stats))
      if (length(indx) > 0) {
        stats_val <- stats::wilcox.test(indx, stats$score_rank)
        pval <- stats_val$p.value
      } else
        pval <- 1.0
      set_enrichment_score <-
        stats_val$statistic#gsEasy::gset(S = pathway_members_in_study, r = rownames(stats), raw_score = T)
    }
    enrichment_stats <-
      rbind(
        enrichment_stats,
        data.frame(
          pathway = current_member,
          pathway_members = paste(as.character(pathway_members_in_study), collapse = ";"),
          pval = pval,
          fdr = NA,
          n = n,
          N = N,
          set_enrichment_score = set_enrichment_score
        )
      )
  }
  # Add q value to the enrichment stats and write it to the output
  enrichment_stats$fdr <-
    stats::p.adjust(enrichment_stats$pval, method = 'BH')
  if (dim(enrichment_stats)[1] == 0)
    stop(
      "No pathway found with minimum number of features!!!\nPlease make sure your mapping files includes your features!!!"
    )
  enrichment_stats <- enrichment_stats |>
    dplyr::arrange(pval)
  return(enrichment_stats)
} #end of function
