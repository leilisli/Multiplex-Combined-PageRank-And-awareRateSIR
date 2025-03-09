
#This repository provides an R function called MCPR (Modified Combined PageRank) 
#that combines closeness centrality from one graph (graph1) with the personalized PageRank 
#computation on another graph (graph2). This function can be useful in scenarios where 
#you want to inject information about one network's structure (via closeness centrality) 
#into the PageRank computation of a different network.


#Overview

#MCPR does the following:
  #  Computes closeness centrality on graph1 using the igraph package (with normalization turned on).
  # Normalizes the resulting closeness values so they sum up to 1. This normalized vector serves as a personalized vector for PageRank.
  #  Runs PageRank on graph2 with the specified damping parameter, using the closeness-based vector from step 2 as its personalization.
  #  Returns the PageRank scores of the vertices in graph2.

#Installation
  #  1.Install R from CRAN if you don't already have it.
  #  2.Install the igraph package in R:   install.packages("igraph")
  #  3.Clone or download this repository and open the R file containing the MCPR function.

#Parameters:
 #   graph1: An igraph object used to compute closeness centrality.
 #   graph2: An igraph object on which the personalized PageRank is calculated.
  #  damping: Damping factor for PageRank (default is 0.85).

#Returns:
  #  A named numeric vector of PageRank scores for each vertex in graph2.

#-------------------------------------------MCPR-----------------------------------------------------#
MCPR <- function(graph1, graph2, damping = 0.85) {
  # Compute closeness centrality for graph1 (normalized)


  closeness_vals <- closeness(graph1, normalized = TRUE)
  
  # Normalize the closeness values so that they sum to 1
  if(sum(closeness_vals) == 0) {
    personalized_vector <- rep(1 / length(closeness_vals), length(closeness_vals))
  } else {
    personalized_vector <- closeness_vals / sum(closeness_vals)
  }
  
  # Compute PageRank on graph2 using the personalized vector from graph1's closeness
  pr <- page.rank(graph2, damping = damping, personalized = personalized_vector)$vector
  
  return(pr)
}

# Example: Create two graphs with identical vertices
g1 <- make_ring(10)
g2 <- make_star(10, mode = "undirected")

# Compute MCPR
MC_pr <- MCPR(g1, g2)
print(MC_pr)
