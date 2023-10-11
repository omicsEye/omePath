# Libraries --------------------------------------------------------------

library(igraph)

omepath_network <-  function(edges){
  # Data Preparation -------------------------------------------------------
  setwd("~/Downloads/data")
  
  edges <- read.delim(
    'data/RYGB_metabolite_weight.txt',
    sep = '\t',
    header = TRUE,
    fill = F,
    comment.char = "" ,
    check.names = F,
    row.names = 1
  )
  #edges <- edges[edges$weight > 1.7,]
  
  
  #Create graph for the algorithms
  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  
  # Community Detection ----------------------------------------------------
  
  # Louvain
  lc <- igraph::cluster_louvain(g)
  lc_com <- igraph::membership(lc)
  #igraph::communities(lc)
  #plot(lc, g)
  
  # Infomap
  imc <- igraph::cluster_infomap(g)
  igraph::membership(imc)
  igraph::communities(imc)
  #plot(lc, g)
  #igraph::plot.igraph(g)
  
  #-------------vizNetwrok
  library(visNetwork)
  
  label = c(unlist(unique(union(edges$from, edges$to))))
  id = seq(length(unlist(unique(union(edges$from, edges$to)))))
  nodes = data.frame(id, label)
  colnames(nodes) <- c("id", "label")
  
  #id has to be the same like from and to columns in edges
  nodes$id <- nodes$label
  nodes$title <- nodes$label
  #Edges
  edges <- edges[,1:3]# as.data.frame(lesmis[1])
  colnames(edges) <- c("from", "to", "width")
  
  #Create graph for Louvain
  graph <- igraph::graph_from_data_frame(edges, directed = FALSE)
  
  
  #Louvain Comunity Detection
  cluster <- igraph::cluster_louvain(graph)
  
  cluster_df <- as.data.frame(as.list(igraph::membership(cluster)))
  cluster_df <- as.data.frame(t(cluster_df))
  cluster_df$label <- rownames(cluster_df)
  
  #Create group column
  nodes <- left_join(nodes, cluster_df, by = "label")
  colnames(nodes)[4] <- "group"
  
  loaded_data <- load_data(input = 'data/Metabolites_netome_format.xlsx',
                           type = 'all', sheet = 1, ID = 'Metabolite')
  features_info <- loaded_data$feature_metadata
  nodes <- nodes[nodes$label %in% features_info$Metabolite,]
  edges <- edges[edges$from %in% features_info$Metabolite & edges$to %in% features_info$Metabolite ,]
  nodes$pathway <- features_info[match(nodes$label, features_info$Metabolite),"SUB_PATHWAY"] #SUPER_PATHWAY
  
  #comment out the next two statement for superpathway 
  pathways_of_interest = c("Aminosugar Metabolism",
                           "Diacylglycerol",
                           "Sphingomyelins"
  )
  nodes$pathway[!nodes$pathway %in% pathways_of_interest] <- NA
  
  nodes$group <- nodes$pathway
  
  
  vizgraph <- visNetwork(nodes, edges, width = "100%") %>%
    visIgraphLayout() %>%
    visNodes(
      shape = "dot",
      color = list(
        background = "#0085AF",
        border = "#013848",
        highlight = "#000000"
      ),
      shadow = list(enabled = F, size = 10)
    ) %>%
    visEdges(
      shadow = FALSE,
      color = list(color = "#0085AF", highlight = "#C62F4B")
    ) %>%
    visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T),
               selectedBy = "group", nodesIdSelection = F) %>% 
    #visLegend(width = 0.1, position = "right", main = "community") %>%
    addFontAwesome() %>%
    addExport(pdf = TRUE)  %>%
    visLayout(randomSeed = 90) %>%
    visExport(type = "pdf", name = "export-network-community",
              float = "left", label = "Save network", background = "white", style= "") 
  vizgraph
  write.table(
    nodes,
    'data/community-network.txt',
    sep = "\t",
    eol = "\n",
    quote = F,
    col.names = NA,
    row.names = T
  )
}

