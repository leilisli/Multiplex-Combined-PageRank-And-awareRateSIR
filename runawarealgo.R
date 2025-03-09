#This README provides an explanation of the runawarealgo function, its workflow, and the input/output parameters. 
#The code is written in R and utilizes helper functions such as awareRateSIR and InitialImmunization (in separate files) 
#to model the process of immunizing (removing/protecting certain nodes) and subsequently simulating disease spread 
#in multiplex networks (physical and virtual layers).

#Table of Contents
    #Introduction
    #Prerequisites
   # Function Explanation: runawarealgo
    #Input Parameters
    #Function Output
    #Usage Workflow
    #Additional Notes

#Introduction
#In epidemic models run on networks, one of the methods to reduce the extent of a disease outbreak is to immunize a fraction of the nodes based on criteria such as centrality or awareness. This code first uses helper functions like awareImmunization to identify the key nodes for immunization, and then uses runAwareRateSIR to estimate the final epidemic size after removing those nodes.
#Prerequisites
    #R installed (version 3.0 or higher is recommended).
    #Necessary packages (such as igraph, foreach, etc.) must be installed and loaded to run helper functions like runAwareRateSIR and awareImmunization.
    #The following files must be present in the correct path so that they can be sourced via source(...):
        #.../awareRateSIR.R
        #.../InitialImmunization.R

#----------------Function Explanation: runawarealgo--------------------------#
#Immunization:
#In each iteration, the function removes a portion of the nodes (based on a specified percentage of all vertices)
#using awareImmunization. These nodes are then excluded from the simulation by setting their state to 3.

#Epidemic Simulation:
#After removing the target nodes, the function calls runAwareRateSIR to compute the epidemic size. 
#This simulation respects the maximum step limit (maxstep) and ultimately reports the total number of infected nodes.

#Main Loop:
#The procedure is repeated for a range of immunization percentages (start to finish, in increments of interval), 
#and the resulting epidemic sizes are stored in the vector res.

#Input Parameters:
#gcontact: Type: igraph ; The contact (physical) network on which the disease spreads.
#gcom: Type: igraph; The communication (virtual) network where awareness spreads.
#VertexCnt: Type:  numeric ; Number of vertices in the network (assumes both networks have the same size).
#interval: Type: numeric; The incremental step for immunization percentage. For example, 0.1 means increasing by 10% at each step.
#start: Type: numeric; The starting percentage of immunization, e.g., 0.0 or 0.1.
#finish: Type: numeric; The ending percentage of immunization, e.g., 0.5 or 0.7.
#iter: Type: numeric; Number of simulation runs passed to runAwareRateSIR.
#maxstep: Type: numeric; Maximum number of time steps in the runAwareRateSIR simulation.
#Beta: Type: numeric; The disease transmission probability in the contact network.
#mu: Type: numeric; The recovery probability for an infected node per time step.
#fi: Type: numeric; A factor adjusting the infection rate based on awareness levels.
#gama: Type: numeric; Probability of forgetting or reducing awareness in each time step.
#alpha: Type: numeric; Probability of awareness spreading in the communication network.

#Note: In the code, numseed is defined as floor(0.02 * VertexCnt), i.e., the initial infected nodes in the epidemic simulation make up 2% of all vertices.

#Function Output
#The function returns a list with two elements:
    #bool: A boolean value set to FALSE at the end.
    #result: A numeric vector of length ((finish - start) / interval) + 1, containing the final epidemic size after immunizing at different percentages.

#Additional Notes
    #Choosing the Immunization Method: The code has a commented line referencing MultiplexPageRank.
    #You can switch to that method instead of awareImmunization depending on your chosen immunization strategy.
    #Performance: For large networks, repeatedly removing nodes and recalculating neighborhoods can be time-consuming.
    #Dependency on Helper Functions: Ensure that awareImmunization and runAwareRateSIR are properly loaded.
#----------------------------------------------------------------#
source('.../awareRateSIR.R')
source('.../InitialImmunization.R')

runawarealgo <- function (gcontact, gcom, VertexCnt, interval, start, finish,
                           iter, maxstep, Beta,
                          mu, fi, gama, alpha, color, bool){
  
  lennum <- ((finish - start ) / interval);
  lennum <- lennum + 1;
  res <- rep(0, lennum );
  
  numseed <- floor(0.02 * VertexCnt)
  
  numseq <- seq(start, finish, interval)
  #numImun <- floor( interval * VertexCnt )
  numImun <- floor(numseq * VertexCnt);
  for ( i in (1:length(res)) ){
    
    ###########################set id attribute for all vertex###############
    Upgcontact = gcontact
    for(n in 1:vcount(Upgcontact)){
      Upgcontact <- set.vertex.attribute(Upgcontact, 'nameid', n, n)
    }
    Upgcom = gcom
    for(m in 1:vcount(Upgcom)){
      Upgcom <- set.vertex.attribute(Upgcom, 'nameid', m, m)
    }
    ###########################################################################
    state = rep( 0, VertexCnt );
    awareness = rep(-1, VertexCnt);
    
    deleted <- awareImmunization( Upgcontact, Upgcom, numImun[i]);
    
     #deleted <-  MultiplexPageRank(numImun[i]);
    
    delname <- rep(0, length(deleted))
    for(j in 1:length(deleted)){
      delname[j] <- get.vertex.attribute(Upgcontact, 'nameid', deleted[j])
    }
    Upgcontact<- Upgcontact-deleted;
    Upgcom <- Upgcom-deleted
    
    state[delname] = 3
    awareness[delname] = 0
    res[i] <- runAwareRateSIR(Upgcontact, Upgcom,state, awareness, maxstep, Beta, mu,
                              fi, gama, alpha, iter, numseed)
    cat("round:")
    cat(i)
    cat(" ", res[i])
    cat("\n")
    
  }
  print("epidemic size after aware immunization of  population is:")

  bool = FALSE
  return( list("bool"=bool, "result"=res) )
}


