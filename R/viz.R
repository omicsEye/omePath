enrichment_plot <- function(
                 stats_table,
                 enrichment_stats,
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
        read.table(
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
        read.table(
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
    stats <- stats[order(-stats[score_col]), ]
  } else{
    stats <- stats[order(stats[score_col]), ]
  }
  stats$score_rank <- seq.int(nrow(stats))
  
  pathways <-
    unique(unlist(lapply(mapper_pathway2feature[pathway_col], as.character)))
  
  
  rank_plots <- NA
  score_plots <- NA
  site_colours <- set_coloures()
  site_names <- set_names()
  #pathways_names <- names(mapper)#unique(mapper_pathway2feature[1])
  rank_plots <-
    vector(mode = "list", length = dim(enrichment_stats)[1])
  score_plots <-
    vector(mode = "list", length = dim(enrichment_stats)[1])
  if (dim(enrichment_stats)[1] == 0)
    return(list(enrichment_stats, rank_plots, score_plots))
  logging::loginfo("Plotting data for %s, %s", "feature", "enrichment")
  pdf(
    paste(output, '/enrichment_plots.pdf', sep = ''),
    width = 2.4,
    height = 2.25,
    onefile = TRUE
  )
  names(rank_plots) <- row.names(enrichment_stats)
  names(score_plots) <- row.names(enrichment_stats)
  for (i in 1:dim(enrichment_stats)[1]) {
    if (is.na(fdr_threshold)) {
      if (as.numeric(enrichment_stats[i, 'pval']) > pval_threshold) {
        next
      }
    } else if (as.numeric(enrichment_stats[i, 'fdr']) > fdr_threshold) {
      next
    }
    #i <- 3
    #print(c("p-value: ", enrichment_stats[i, 'pval']))
    pathway_members_in_study <-
      intersect(mapper_pathway2feature[as.character(mapper_pathway2feature[, pathway_col]) == as.character(enrichment_stats[i, 'pathway']), feature_col], rownames(stats))
    indx <- match(pathway_members_in_study, rownames(stats))
    
    tryCatch({
      stats$Set <- 'all_features'
      stats$Set[indx] <- 'pathway'
      
      #set the color of rugs only for member of pathway o interest
      stats$Rug <- NA
      stats$Rug[indx] <- 'pathway'
      
      ### plot the enrichment based on score_rank ####################
      density_plot <-
        ggplot2::ggplot(stats, ggplot2::aes(
          x = stats$score_rank,
          fill = Set,
          color = Set
        ))
      density_plot <- density_plot +
        ggplot2::geom_density(alpha = 0.4, size = .15) +
        ggplot2::geom_rug(ggplot2::aes(x = stats$score_rank, color = Rug, y = 0),
                          alpha = 0.4,
                          size = .1) +
        ggplot2::scale_color_manual(
          values = unlist(site_colours),
          labels = unlist(site_names),
          name = "Set"
        ) +
        ggplot2::scale_fill_manual(
          values = unlist(site_colours),
          labels = unlist(site_names),
          name = "Set"
        ) + #theme(legend.position="none") +
        ggplot2::labs(title = enrichment_stats[i, 'pathway']) + theme_omicsEye() +
        #ggplot2::geom_hline(yintercept = set_enrichment_score, colour = 'black',  linetype = "dashed", size =.5) +
        ggplot2::guides(
          fill = ggplot2::guide_legend(
            title = "",
            keywidth = 0.25 ,
            keyheight = 0.25,
            default.unit = "cm"
          ),
          colour = ggplot2::guide_legend(
            title = "",
            keywidth = 0.25 ,
            keyheight = 0.25,
            default.unit = "cm"
          )
        ) +
        ggplot2::theme(legend.justification = c(0, 0),
                       legend.position = c(.15, .7)) +
        ggplot2::annotate(
          geom = "text",
          x = Inf,
          y = Inf,
          hjust = 1,
          vjust = 1,
          #x=0,  y=Inf, vjust=2,
          label = sprintf(
            "p-value: %.4f\nn: %s out of %s\n%s",
            enrichment_stats[i, 'pval'],
            enrichment_stats[i, 'n'],
            enrichment_stats[i, 'N'],
            Pathway.Subject
          ) ,
          color = "black",
          size = 2.25,
          fontface = "italic"
        ) +
        ggplot2::annotate(
          geom = "text",
          x = median(stats$score_rank),
          y = Inf,
          hjust = 1,
          vjust = 1,
          #x=0,  y=Inf, vjust=2,
          label = sprintf("Control") ,
          color = "black",
          size = 2.25,
          fontface = "italic"
        ) +
        ggplot2::geom_vline(ggplot2::aes(xintercept = median(stats$score_rank)),
                            color = "black",
                            size = 0.1) +
        
        ggplot2::xlab(sprintf("Rank of %s", score_col)) +
        ggplot2::ylab('Density')
      rank_plots[[enrichment_stats[i, 'pathway']]] <- density_plot
      #ggsave(filename=paste(outpath,'/score_',enrichment_stats[i, 'pathway'],'.pdf', sep =''),
      #       plot=density_plot, width = 60, height = 50, units = "mm", dpi = 350)
      stdout <-
        capture.output(print(density_plot), type = "message")
      #logging::logdebug(stdout)
      
      ### plot the enrichment based on score ####################
      density_plot <-
        ggplot2::ggplot(stats, ggplot2::aes(
          x = get(score_col),
          fill = Set,
          color = Set
        ))
      density_plot <- density_plot +
        ggplot2::geom_density(alpha = 0.4, size = .15) +
        ggplot2::geom_rug(ggplot2::aes(
          x = get(score_col),
          color = Rug,
          y = 0
        ),
        alpha = 0.4,
        size = .1) +
        ggplot2::scale_color_manual(
          values = unlist(site_colours),
          labels = unlist(site_names),
          name = "Set"
        ) +
        ggplot2::scale_fill_manual(
          values = unlist(site_colours),
          labels = unlist(site_names),
          name = "Set"
        ) + #theme(legend.position="none") +
        ggplot2::labs(title = enrichment_stats[i, 'pathway']) + theme_omicsEye() +
        #ggplot2::geom_hline(yintercept = set_enrichment_score, colour = 'black',  linetype = "dashed", size =.5) +
        ggplot2::guides(
          fill = ggplot2::guide_legend(
            title = "",
            keywidth = 0.25 ,
            keyheight = 0.25,
            default.unit = "cm"
          ),
          colour = ggplot2::guide_legend(
            title = "",
            keywidth = 0.25 ,
            keyheight = 0.25,
            default.unit = "cm"
          )
        ) +
        ggplot2::theme(legend.justification = c(0, 0),
                       legend.position = c(.15, .7)) +
        ggplot2::annotate(
          geom = "text",
          x = Inf,
          y = Inf,
          hjust = 1,
          vjust = 1,
          #x=0,  y=Inf, vjust=2,
          label = sprintf(
            "p-value: %.4f\nn: %s out of %s\n%s",
            enrichment_stats[i, 'pval'],
            enrichment_stats[i, 'n'],
            enrichment_stats[i, 'N'],
            Pathway.Subject
          ) ,
          color = "black",
          size = 2.25,
          fontface = "italic"
        ) +
        ggplot2::annotate(
          geom = "text",
          x = 0.0,
          y = Inf,
          hjust = 1,
          vjust = 1,
          #x=0,  y=Inf, vjust=2,
          label = sprintf("Control") ,
          color = "black",
          size = 2.25,
          fontface = "italic"
        ) +
        ggplot2::geom_vline(xintercept = 0,
                            color = "black",
                            size = 0.1) +
        ggplot2::xlab(score_col) +
        ggplot2::ylab('Density')
      score_plots[[enrichment_stats[i, 'pathway']]] <- density_plot
      #ggsave(filename=paste(outpath,'/score_',enrichment_stats[i, 'pathway'],'.pdf', sep =''),
      #       plot=density_plot, width = 60, height = 50, units = "mm", dpi = 350)
      stdout <-
        capture.output(print(density_plot), type = "message")
      #logging::logdebug(stdout)
      
      #}
      # pdf(NULL)
      # dev.control(displaylist="enable")
      # barcodeplot(stats$score_col, index=indx,
      #             main=paste(enrichment_stats[i, 'pathway'], "\nP-value: ",
      #                        enrichment_stats[i, 'pval'],', FDR:', enrichment_stats[i, 'fdr'], ', N: ',
      #                        enrichment_stats[i, 'n'],  sep=''),
      #             labels = c("Low","High"), xlab = barcodeplot_xlable, #col.bars = c('red', 'blue'),#'darkolivegreen4',#
      #             alpha = 1, span.worm = 0.45, quantiles = c(.0,.0))
      # t  <- recordPlot()
      # invisible(dev.off())
    }, error = function(e) {
      message(e)
    })
    
  } # end of the loop for plotting
  invisible(dev.off())
  saveRDS(rank_plots, file = paste(output,"/figures/",  "gg_enrichment_rank.RDS", sep = ""))
  saveRDS(score_plots, file = paste(output,"/figures/", "gg_enrichment_score.RDS", sep = ""))
  
  plot_results <- list()
  plot_results$rank_plots <- rank_plots
  plot_results$score_plots <- score_plots
  return(plot_results)
} #end of function