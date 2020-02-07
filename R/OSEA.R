
# load in the required libraries, report an error if they are not installed
for( lib in c('gsEasy','future', 'limma', 'ggplot2')) {
  if(! suppressPackageStartupMessages(require(lib, character.only=TRUE)) ) stop(paste("Please install the R package: ",lib))
}

plan(multisession)

OSEA <- function(stats_table,
                 output = "~/",
                 score_col = 'logFC',
                 pval_threshold = 0.05,
                 fdr_threshold = NA,
                 Pathway.Subject = 'Metabolic',
                 method = 'gset',
                 min_member=2,
                 mapper_file=NA,
                 do_plot = TRUE,
                 pathway_col = "Pathway",
                 feature_col = "Feature"){
  
  # load mapping files
  ## load all elements (HMDBID) in all sets (Ko's pathways)
  if (is.na(mapper_file)) {
    mapper_pathway2feature <- read.csv("data/smpdb_metabolites.csv")
  }else if (endsWith(mapper_file,".csv")){
    mapper_pathway2feature <- read.csv(mapper_file)
  }else {
    mapper_pathway2feature <- read.table(mapper_file, header = T, sep = '\t')# col.names = NA, row.names = T
  }
  if (!is.na(Pathway.Subject)) {
    if ('Subject' %in% colnames(mapper_pathway2feature))
      mapper_pathway2feature <- mapper_pathway2feature[mapper_pathway2feature$Subject == Pathway.Subject,]
    else print(paste('Pathway.Subject is not mapping file!!!'))
  }
  ## load mapping old<-->new HMDB IDs files
  #if (!is.na(mapper_oldhmdb_newhmdb_file)) {
  #  mapper_oldhmdb_newhmdb <- read.csv(mapper_oldhmdb_newhmdb_file)
  #}
  #else{
  #  mapper_oldhmdb_newhmdb <- read.csv("data/ParsedHMDB_v4.0.csv")
  #}
  
  # load the element file for the first stats
  # if a character string then this is a file name, else it is a data frame
  if (is.character(stats_table)) {
    stats <- data.frame(data.table::fread(stats_table, header = TRUE, sep = "\t"), row.names = 1)
    if (nrow(stats_table) == 1) {
      # read again to get row name
      stats <- read.table(stats_table, header = TRUE, row.names = 1,sep = "\t", fill = FALSE,
                          comment.char = "" , check.names = FALSE)
    }
  } else {
    stats <- stats_table
  }
  stats <- stats[!is.na(stats[score_col]),,drop = F]
  # setup stats file for features in the study
  #stats <- data.frame(score_col = c(), Rank = c())
  #stats <- rbind(stats, data.frame(score_col = stats[score_col]))
  #colnames(stats) <- c(score_col)
  stats$Rank <- NA
  #if (nchar(stats$feature[1]) == 9 )
  #  stats$feature <- mapper_oldhmdb_newhmdb[match(stats$feature, mapper_oldhmdb_newhmdb$OldHMDBIDs), 'HMDBIDs']
  #stats$feature <- rownames(stats)
  #stats$Rank <- NA
  if (score_col == 'P.Value') {
    stats <- stats[order(-stats[score_col]),]
  }else{
    stats <- stats[order(stats[score_col]),]
  }
  stats$Rank <- seq.int(nrow(stats))
  
  
  #stats <- stats#[sample(1:dim(stats)[1],200000,replace=F),]
  enrichment_stats <- data.frame(pathway = c(), pathway_members = c(),
                                 pval = c(), fdr = c(), n = c(), N = c(), set_enrichment_score = c())
  
  # calculate p value for all sets (pathway terms)
  logging::loginfo("Loading files is done, now calculating enrichments p-value! This may take 10 or more minutes!!!")
  pathways <- unique(unlist(lapply(mapper_pathway2feature$Pathway, as.character)))
  #pathways <- intersect(mapper_pathway2feature$Pathway, mapper_pathway2feature$Pathway)
  #i <- 1
  for (i in 1:length(pathways)) {
    current_member <- pathways[i] # "Purine Metabolism" #
    pathway_members <- mapper_pathway2feature[as.character(mapper_pathway2feature[,pathway_col]) == current_member, feature_col]
    pathway_members_in_study <- intersect(pathway_members, rownames(stats))
    N <- length(pathway_members)
    n <- length(pathway_members_in_study)
    if (n < min_member) next
    #print(c(i, current_member, length(pathway_members_in_study), min_member))
    #print(sprintf("Number of elements(HMDBID) in the set (pathway term): %s",length(pathway_members_in_study)))
    if (method == 'ks') {
      indx <- match(pathway_members_in_study, rownames(stats))
      if (length(indx) > 0) {
        stats_val <- ks.test(indx, stats$Rank)
        pval <- stats_val$p.value
      }else pval <- 1.0
      set_enrichment_score <- gsEasy::gset(S = pathway_members_in_study, r = rownames(stats), raw_score = T)
    }else if (method == 'gset') {
      pval <- gsEasy::gset(S = pathway_members_in_study, r = rownames(stats))
      set_enrichment_score <- gsEasy::gset(S = pathway_members_in_study, r = rownames(stats),
                                           raw_score = T)
    }else if (method == 'wilcox') {
      indx <- match(pathway_members_in_study, rownames(stats))
      if (length(indx) > 0) {
        stats_val <- wilcox.test(indx, stats$Rank)
        pval <- stats_val$p.value
      }else pval <- 1.0
      set_enrichment_score <- stats_val$statistic#gsEasy::gset(S = pathway_members_in_study, r = rownames(stats), raw_score = T)
    }
    enrichment_stats <- rbind(enrichment_stats, data.frame(pathway = current_member,
                                                           pathway_members = paste(as.character(pathway_members_in_study), collapse = ";"),
                                                           pval = pval, fdr = NA,
                                                           n = n, N = N,
                                                           set_enrichment_score = set_enrichment_score))
  }
  # Add q value to the enrichment stats and write it to th eoutput
  enrichment_stats$fdr <- p.adjust(enrichment_stats$pval, method = 'BH')
  if(dim(enrichment_stats)[1] == 0) stop("No pathwasy found with minumu number of features!!!\nPlease make sure your your mapping files includes your features!!!")
  enrichment_stats <- enrichment_stats[order(enrichment_stats["pval"]),]
  # write stats in the output
  write.table(enrichment_stats,  paste(output, '/enrichment_stats.txt', sep=''),
              sep = "\t", eol = "\n", col.names = NA, row.names = T)
  #make barcodeplot
  
  #print("Now time for  barcodeplots!")
  #print(enrichment_stats)
  ## variables for plots
  rank_plots <- NA
  score_plots <- NA
  site_colours <- deepath::set_coloures()
  site_names <- deepath::set_names()
  #pathways_names <- names(mapper)#unique(mapper_pathway2feature[1])
  rank_plots <- vector(mode = "list", length = dim(enrichment_stats)[1])
  score_plots <- vector(mode = "list", length = dim(enrichment_stats)[1])
  if (dim(enrichment_stats)[1] == 0) return(list(enrichment_stats, rank_plots, score_plots))
  logging::loginfo("Plotting data for %s, %s", "feature", "enrichment")
  pdf(paste(output,'/enrichment_plots.pdf', sep = ''), width = 2.5, height = 2.25, onefile = TRUE)
  names(rank_plots) <- row.names(enrichment_stats)
  names(score_plots) <- row.names(enrichment_stats)
  for (i in 1:dim(enrichment_stats)[1]) {
    if (is.na(fdr_threshold)) {
      if (as.numeric(enrichment_stats[i, 'pval']) > pval_threshold) {
        next
      }
    }else if (as.numeric(enrichment_stats[i, 'fdr']) > fdr_threshold) {
      next
    }
    #i <- 3
    #print(c("p-value: ", enrichment_stats[i, 'pval']))
    pathway_members_in_study <- intersect(mapper_pathway2feature[as.character(mapper_pathway2feature[,pathway_col]) == as.character(enrichment_stats[i, 'pathway']), feature_col], rownames(stats))
    indx <- match(pathway_members_in_study, rownames(stats))
    #print(c("N: ", enrichment_stats[i, 'n']) )
    
    tryCatch({
      #stats_val <- ks.test(indx, stats$Rank)
      #print(stats_val$p.value)
      #if (stats_val$p.value<0.000001){
      stats$Set <- 'all_features'
      stats$Set[indx] <- 'pathway'
      
      ### plot the enrichment based on rank ####################
      density_plot <- ggplot2::ggplot(stats)
      density_plot <- density_plot +
        ggplot2::geom_density(ggplot2::aes(x=Rank, fill = Set), alpha = 0.4) +
        ggplot2::geom_rug(ggplot2::aes(x = Rank, color = Set, y = 0),  alpha = 0.4, size =.1 )+
        ggplot2::scale_color_manual(values = unlist(site_colours), labels = unlist(site_names), name="Set")+
        ggplot2::scale_fill_manual(values = unlist(site_colours), labels = unlist(site_names), name="Set") + #theme(legend.position="none") +
        ggplot2::labs(title = enrichment_stats[i, 'pathway'] )+ theme_nature()+
        #ggplot2::geom_hline(yintercept = set_enrichment_score, colour = 'black',  linetype = "dashed", size =.5) +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "", keywidth=0.25 ,keyheight=0.25, default.unit="cm"),
                        colour = ggplot2::guide_legend(title = "",keywidth=0.25 ,keyheight=0.25, default.unit="cm"))+
        ggplot2::theme(legend.justification=c(0,0), legend.position=c(.15,.7))+
        ggplot2::annotate(geom="text", x= Inf, y = Inf, hjust=1,vjust=1,#x=0,  y=Inf, vjust=2,
                          label=sprintf("p-value: %.4f\nn: %s out of %s",enrichment_stats[i, 'pval'], enrichment_stats[i,'n'],
                                        enrichment_stats[i,'N']) ,
                          color="black", size= 2.25, fontface="italic")+
        ggplot2::annotate(geom="text", x= median(stats$Rank), y = Inf, hjust=1,vjust=1,#x=0,  y=Inf, vjust=2,
                          label=sprintf("Control") ,
                          color="black", size= 2.25, fontface="italic")+
        ggplot2::geom_vline( ggplot2::aes(xintercept = median(stats$Rank)), color="black", size = 0.1)+

        ggplot2::xlab(sprintf("Rank of %s", score_col))+
        ggplot2::ylab('Density')
      rank_plots[[enrichment_stats[i, 'pathway']]] <- density_plot
      #ggsave(filename=paste(outpath,'/score_',enrichment_stats[i, 'pathway'],'.pdf', sep =''),
      #       plot=density_plot, width = 60, height = 50, units = "mm", dpi = 350)
      stdout <- capture.output(print(density_plot),type="message")
      #logging::logdebug(stdout)
      
      ### plot the enrichment based on score ####################
      density_plot <- ggplot2::ggplot(stats)
      density_plot <- density_plot +
        ggplot2::geom_density(ggplot2::aes(x=get(score_col), fill = Set), alpha = 0.4) +
        ggplot2::geom_rug(ggplot2::aes(x = get(score_col), color = Set, y = 0),  alpha = 0.4, size =.3 )+
        ggplot2::scale_color_manual(values = unlist(site_colours), labels = unlist(site_names), name="Set")+
        ggplot2::scale_fill_manual(values = unlist(site_colours), labels = unlist(site_names), name="Set") + #theme(legend.position="none") +
        ggplot2::labs(title = enrichment_stats[i, 'pathway'] )+ theme_nature()+
        #ggplot2::geom_hline(yintercept = set_enrichment_score, colour = 'black',  linetype = "dashed", size =.5) +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "", keywidth=0.25 ,keyheight=0.25, default.unit="cm"),
                        colour = ggplot2::guide_legend(title = "",keywidth=0.25 ,keyheight=0.25, default.unit="cm"))+
        ggplot2::theme(legend.justification=c(0,0), legend.position=c(.15,.7))+
        ggplot2::annotate(geom="text", x= Inf, y = Inf, #hjust=1,vjust=1,#x=0,  y=Inf, vjust=2,
                          label=sprintf("p-value: %.4f\nn: %s out of %s",enrichment_stats[i, 'pval'], enrichment_stats[i,'n'],
                                        enrichment_stats[i,'N']) ,
                          color="black", size= 2.25, fontface="italic")+
        ggplot2::annotate(geom="text", x= 0, y = Inf, #hjust=1,vjust=1,#x=0,  y=Inf, vjust=2,
                          label=sprintf("Control") ,
                          color="black", size= 2.25, fontface="italic")+
        ggplot2::geom_vline( xintercept = 0, color="black", size = 0.1)+
        ggplot2::xlab(score_col)+
        ggplot2::ylab('Density')
      score_plots[[enrichment_stats[i, 'pathway']]] <- density_plot
      #ggsave(filename=paste(outpath,'/score_',enrichment_stats[i, 'pathway'],'.pdf', sep =''),
      #       plot=density_plot, width = 60, height = 50, units = "mm", dpi = 350)
      stdout <- capture.output(print(density_plot),type="message")
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
    }, error=function(e){
      message(e)
    })
    
  }# end of the loop for barcodeplot
  invisible(dev.off())
  result <- list()
  result$enrichment_stats <- enrichment_stats
  result$rank_plots <- rank_plots
  result$score_plots <- score_plots
  return (result)
}#end of function


