#' Plotting the visual representation of the order of complex events.
#' This function generates a visual representation of the order of complex events using the bnlearn package.
#' It creates a plot where nodes represent variables (complex events) and their relationships are depicted through directed edges.
#' Optional highlighting of specific nodes is supported.
#' @param res A bnlearn object representing the learned structure of complex events.
#' @export
ooe.plot <- function(res) {
  library(bnlearn)
  nodes <- c()
  for (i in 1:length(res$nodes)) {
    nodes <- c(nodes, res$nodes[[i]]$mb)
    nodes <- unique(nodes)
  }
  
  inter_nodes <- c()
  for (i in 1:length(nodes)) {
    if (strsplit(nodes[i], "\\.")[[1]][1] != strsplit(nodes[i], "\\.")[[1]][3])
    {
      inter_nodes <- c(inter_nodes, nodes[i])
    }
  }
  if (length(inter_nodes) > 0) {
    highlight_list <- list(nodes = inter_nodes, fill = "orange")
    graphviz.plot(res, shape = "rectangle", highlight = highlight_list)
  } else {graphviz.plot(res, shape="rectangle")}
}
