#This README provides an overview of the functions included in the code snippet and instructions on how to use them.
# All of these functions are written in R and rely on igraph for graph/network functionalities. 
#The code snippet also references a custom function from .../MCPR.R, so make sure you have that file available as well.

#Table of Contents

    #Prerequisites
    #Installation
    #Usage
    #Functions Overview
        #1. InitialImmunization(net, measure, res)
        #2. centrality(tmpG, measure)
        #3. classicImmunization(graphContact, measure1, numImun)
        #4. awareImmunization(graphContact, graphCom, numImun)
        #5. MultiplexPageRank(numImun)
    #Additional Notes

#Prerequisites

   # R (version 3.0 or higher recommended)
   # igraph R package
    #The MCPR.R script located in the Codes directory (the code calls source('.../MCPR.R')).

#Installation

    #Install R from The Comprehensive R Archive Network (CRAN).

    #Open R (or RStudio) and install the igraph package:
#install.packages("igraph")
#library(igraph)

#Make sure your project directory structure matches the code references. By default, it expects:
#Codes/ MCPR.R

#Place the code snippet you have in an R file (e.g., immunization.R), 
#then source('immunization.R') to load all its functions in your environment.

#Usage

    #Load the MCPR script:
#source('.../MCPR.R')
#Load the immunization code (assuming you saved it as immunization.R):
#source('immunization.R')

#Create or load your networks as igraph objects in R.
#Call the desired function by passing the appropriate arguments (see below for function details).

#Functions Overview
#1. InitialImmunization(net, measure, res)

    #Description:
    #Takes an igraph network net, calculates a chosen centrality measure for each node, ranks the nodes by that measure in descending order, and removes the top res nodes from the network.
    #Parameters:
        #net: An igraph network.
        #measure: An integer specifying which centrality measure to use. (See centrality for details.)
        #res: Number of nodes to immunize (remove).
    #Returns:
    #A modified igraph network (currentNet) after removing the top res nodes by the given centrality measure.

#-----------centrality(tmpG, measure)-----------#

    #Description:
    #A helper function to compute various centralities.
    #Parameters:
       # tmpG: An igraph graph.
        #measure: An integer code for the centrality measure:
           # Degree centrality
            #Betweenness centrality
            #Closeness centrality
           # Eigenvector centrality
            #PageRank centrality
    #Returns:
    #A numeric vector containing the centrality scores for all vertices in tmpG.

#----------classicImmunization(graphContact, measure1, numImun)--------#

    #Description:
    #Computes the chosen centrality measure (via centrality) for graphContact, sorts the nodes by their centrality in descending order, and returns the top numImun nodes (i.e., the nodes to immunize).
   # Parameters:
        #graphContact: An igraph network.
       # measure1: Centrality measure (same encoding as in centrality).
       # numImun: Number of nodes to immunize.
    #Returns:
    #A numeric vector (node IDs) representing the top numImun nodes.


#--------------awareImmunization(graphContact, graphCom, numImun)-------------#

    #Description:
    #Uses the custom MCPR function (sourced from .../MCPR.R) to calculate a multiplex-based ranking. It then removes the top numImun nodes from the ranking.
    #Parameters:
       # graphContact: An igraph network (physical/first layer).
        #graphCom: Another igraph network (virtual/second layer).
       # numImun: Number of nodes to immunize.
    #Returns:
    #A numeric vector (node IDs) representing the top numImun nodes to immunize based on the multiplex PageRank (MCPR).

#--------------MultiplexPageRank(numImun)------------#

    #Description:
    #Loads pre-computed centrality values from a CSV file (centrality.csv by default at C:/MyCode/DataSet/centrality.csv), then ranks them in descending order and returns the top numImun nodes.
    #Parameters:
       # numImun: Number of nodes to immunize (i.e., the top numImun entries from the CSV).
    #Returns:
    #A numeric vector (node IDs).
    #Note:
    #You may need to update the file path inside the function read.csv() to match where your file is located.
#------------------------------------------------------------------------#
#Additional Notes

    #File Paths:
    #In MultiplexPageRank(), the path to centrality.csv is hardcoded as "C:/MyCode/DataSet/centrality.csv". Update this to match the location of your CSV file.

    #Dependencies:
    #Make sure you have installed and loaded the igraph package. If you are using awareImmunization(), ensure that MCPR.R is properly sourced and that it does not have additional unmet dependencies.

    #Data Requirements:
    #Each function assumes the graphs are valid igraph objects with consistent node labeling. If you work with large networks, watch out for performance implications when calculating certain centrality measures (especially betweenness).

#The centrality.csv file contains the MultiplexPR values of the nodes, which you can obtain using the multiplexPageRank function in the multiplexPageRank.m file.
#-------------------------------------------------------------------------------------------------------#
source('.../MCPR.R')

InitialImmunization<- function(net,measure,res){
  #res=numImun
  CVec<-centrality(net,measure);
  CIndex<-order(CVec,decreasing=TRUE);
  ImmuneCandidate<-CIndex[1:(res)];
  currentNet<-net-ImmuneCandidate;
  return(currentNet);		
}

centrality <-function(tmpG,measure){
  if(measure==1)#degree centrality
    CVec<-degree(tmpG, v=V(tmpG), mode = c("all"),loop=FALSE);
  if(measure==2)#betweenness centrality
    CVec<-betweenness(tmpG, v=V(tmpG), directed = TRUE, weights = NULL);
  if(measure==3)#closeness centrality
    CVec<-closeness(tmpG, vids=V(tmpG), mode = c("all"));
  if(measure==4)#eigenvector centrality
    CVec<-evcent(tmpG)$vector;
  if(measure==5)#PageRank centraliy
    CVec<-page.rank(tmpG,directed = FALSE)$vector;	
     
  return (CVec)
}


classicImmunization <- function(graphContact, measure1, numImun){
  vcnt = vcount(graphContact)
  CVec <- centrality( graphContact, measure1 )
  CIndex <- order( CVec, decreasing = TRUE )
  deleted <- vector(mode="numeric", length=0)
  deleted = CIndex[1:numImun]
  #print(numImun)
  return(deleted)
}


awareImmunization <- function(graphContact, graphCom, numImun) {
  vcnt <- min(vcount(graphContact), vcount(graphCom))
  rankingVec <- rep(0, vcnt)
  deleted <- vector(mode = "numeric", length = 0)
  rankingVec <- MCPR(graphContact, graphCom, damping = 0.85)
  rankingIndex <- order(rankingVec, decreasing = TRUE)
  deleted <- rankingIndex[1:numImun]

  return(deleted)
}


MultiplexPageRank <- function(numImun)
{
  #return(ImmuneCandidate)
  
  MultiPlexPR <- read.csv("C:/MyCode/DataSet/centrality.csv",header = FALSE)
  NodeRankingMPR <- order( MultiPlexPR, decreasing = TRUE )
  deleted <- vector(mode="numeric", length=0)
  deleted = NodeRankingMPR[1:numImun]
  return(deleted)
}
