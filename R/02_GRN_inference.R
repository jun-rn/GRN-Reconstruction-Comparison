library(minet)
library(igraph)
library(ggraph)

expr_t <- t(as.matrix(expr_small))

#mutual information matrix
mi_matrix <- build.mim(expr_t, estimator = "spearman")

grn_relevance <- clr(mi_matrix)
#clr = content likelihood of relatedness
#relevance network approach

#convert to edge list
get_edges <- function(mat, method_name) {
  which(mat>0, arr.ind = TRUE) %>%
    as.data.frame() %>%
    mutate(
      gene1 = rownames(mat)[row],
      gene2 = colnames(mat)[col],
      weight = mat[cbind(row,col)],
      method = method_name
    ) %>%
    dplyr::select(gene1,gene2,weight,method) %>%
    filter(gene1 < gene2)
}

edges_relevance <- get_edges(grn_relevance, "Relevance Network")
cat("Relevance Network edges:", nrow(edges_relevance), "\n")

grn_mi <- mrnet(mi_matrix)
#mrnet = maximum relevance/minimum redundancy network
#represents pure mutual information network

edges_mi <- get_edges(grn_mi, "Mutual Information")
cat("Mutual Information edges:", nrow(edges_mi), "\n")

grn_aracne <- aracne(mi_matrix, eps = 0.05)

edges_aracne <- get_edges(grn_aracne, "ARACNE")
cat("ARACNE edges:", nrow(edges_aracne), "\n")

#Filter top edges
edges_relevance_top <- edges_relevance %>% arrange(desc(weight)) %>% head(5000)
edges_mi_top        <- edges_mi %>%        arrange(desc(weight)) %>% head(5000)
edges_aracne_top    <- edges_aracne %>%    arrange(desc(weight)) %>% head(5000)