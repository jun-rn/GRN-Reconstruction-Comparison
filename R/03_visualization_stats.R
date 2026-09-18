#Build igraph objects for each
g_relevance <- graph_from_data_frame(edges_relevance_top, directed = FALSE)
g_mi <- graph_from_data_frame(edges_mi_top, directed = FALSE)
g_aracne <- graph_from_data_frame(edges_aracne_top, directed = FALSE)

#save relevance network
jpeg("relevance_network.jpg", width = 800, height = 800, quality = 100)
plot(g_relevance,
     main = "Relevance Network",
     vertex.size = 3,
     vertex.label = NA,
     edge.color = "steelblue",
     vertex.color = "tomato")
dev.off()

#save mutual information network
jpeg("mutual_information.jpg", width = 800, height = 800, quality = 100)
plot(g_mi,
     main = "Mutual Information",
     vertex.size = 3,
     vertex.label = NA,
     edge.color = "darkgreen",
     vertex.color = "tomato")
dev.off()

#save ARACNE Network
jpeg("aracne_network.jpg", width = 800, height = 800, quality = 100)
plot(g_aracne,
     main = "ARACNE Network",
     vertex.size = 3,
     vertex.label = NA,
     edge.color = "purple",
     vertex.color = "tomato")
dev.off()

compare_stats <- data.frame(
  Method = c("Relevance Network", "Mutual Information", "ARACNE"),
  Nodes = c(vcount(g_relevance), vcount(g_mi), vcount(g_aracne)),
  Edges = c(ecount(g_relevance), ecount(g_mi), ecount(g_aracne)),
  Density = c(
    edge_density(g_relevance),
    edge_density(g_mi),
    edge_density(g_aracne)
    
  ),
  Avg_degree = c(
    mean(degree(g_relevance)),
    mean(degree(g_mi)),
    mean(degree(g_aracne))
  ),
  Clustering_coefficient = c(
    transitivity(g_relevance),
    transitivity(g_mi),
    transitivity(g_aracne)
  ),
  Avg_path_length = c(
    mean_distance(g_relevance),
    mean_distance(g_mi),
    mean_distance(g_aracne)
  )
)

print(compare_stats)