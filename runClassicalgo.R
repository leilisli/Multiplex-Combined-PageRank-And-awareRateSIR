#This README provides an overview of the awareRateSIR model implementation and the associated runAwareRateSIR function in R. 
#These functions simulate an SIR-like epidemic spreading process on a multiplex network 
#(composed of a contact layer and a communication layer) that also incorporates awareness dynamics.

#Table of Contents

    #1.Introduction
    #2.Functions Overview
        #1. awareRateSIR(...)
        #2. runAwareRateSIR(...)
    #3.Parameters and Their Meaning
    #4.Usage Example
    #5. Additional Notes

#Introduction

#In many real-world scenarios, individuals (nodes) can become aware of a disease (or other spreading phenomenon) 
#through both direct contact (physical/in-person) and communication or social media (virtual/online). 
#This code extends a standard SIR model to account for:

    #Disease spreading on the contact network (graphContact).
    #Awareness spreading on the communication network (graphCom).
    #State transitions that reflect infection, awareness, and recovery processes.
   # Stochastic processes (using probabilities) for both infection and awareness transmission.

#By integrating these mechanisms, the awareRateSIR function and its runner runAwareRateSIR 
#allow us to observe how awareness influences the spread of infection and vice versa.

#Introduction

#In many real-world scenarios, individuals (nodes) can become aware
 #of a disease (or other spreading phenomenon) through both direct contact (physical/in-person) 
#and communication or social media (virtual/online). This code extends a standard SIR model to account for:

    #Disease spreading on the contact network (graphContact).
    #Awareness spreading on the communication network (graphCom).
    #State transitions that reflect infection, awareness, and recovery processes.
    #Stochastic processes (using probabilities) for both infection and awareness transmission.

#By integrating these mechanisms, the awareRateSIR function and its runner runAwareRateSIR 
#allow us to observe how awareness influences the spread of infection and vice versa.

#Functions Overview
#--------------------1. awareRateSIR(graphContact, graphCom, state, awareness, maxstep, Beta, mu, fi, gama, alpha, numseed)----------------#

    #Purpose:
    #Simulates one instance of the SIR-like model with awareness dynamics for up to maxstep discrete time steps.

   # Procedure (in brief):
        #Seed Infection: Randomly chooses numseed nodes to infect initially.
        #Iterative Updates:
           # Infection Spread on graphContact with probability Beta, affected by an individual’s awareness level (reducing effective infection rate).
           # Recovery occurs with probability mu.
           # Awareness Spread on graphCom with probability alpha.
            #Awareness Fading: aware nodes lose awareness (or become less aware) with probability gama.
        #Track Epidemic Size (epidemicsize) as the number of infections (cumulative).
        #Stop Condition: Runs until either there are no more infected nodes or maxstep steps have been reached.

    #Returns:
        #A single numeric value (epidemicsize), indicating how many nodes got infected (cumulatively) by the time the simulation ends.

#------------------------2. runAwareRateSIR(graphContact, graphCom, state, awareness, maxstep, Beta, mu, fi, gama, alpha, iter, numseed)---------------------#

    #Purpose:
    #Runs multiple iterations of the awareRateSIR simulation to gather average results.

    #Procedure:
        #Parallel / Loop Execution: Repeats the awareRateSIR process iter times (potentially in parallel if you have a parallel backend configured via foreach and %dopar%).
        #Collect Results: Aggregates the epidemic sizes from all runs into tp.
        #Compute Final Output:
            #Takes the average of the epidemic sizes that exceed a small threshold (e.g., 0.001 * vcount(...)) to avoid runs with negligible spread.
            #If no runs exceed that threshold, returns 0.

    #Returns:
       # A single numeric value: the mean epidemic size (above the threshold) across iter runs.
#------------------------------------------------------------------------------#
#Parameters and Their Meaning

#Below is a brief explanation of the main parameters you’ll see in both functions:

    #graphContact:
       # An igraph object representing the contact (physical) network. Edges here facilitate infection.

    #graphCom:
        #An igraph object representing the communication (virtual) network. Edges here facilitate the spread of awareness.

    #state:
       # A numeric vector of the same length as the number of nodes (in both networks). Possible values typically mean:
            #0 = Susceptible (unaware and uninfected)
            #1 = Susceptible (aware but uninfected)
            #2 = Infected
           # 3 = Recovered (or removed)

    #awareness:
        #A numeric vector (same size as state), tracking how aware each node is. Negative values often indicate a freshly infected state or an awareness cooldown, while non-negative values indicate levels or timers of awareness.

    #maxstep:
       # Maximum number of simulation steps.

    #Beta:
       # Probability (base rate) of infection transmission on contact edges.

    #mu:
        #Probability that an infected node recovers in a given time step.

   # fi (or phi):
        #Factor to modify the infection rate among aware individuals.

   # gama (gamma):
        #Probability that an aware node’s awareness fades or decrements at each time step.

    #alpha:
    #Probability of awareness being transmitted from one aware node to its neighbors in graphCom.

    #numseed:
    #Number of initially infected nodes chosen randomly at the start of each simulation.

    #iter:
   # Number of simulation runs (for runAwareRateSIR).

#Additional Notes

    #Parallelization:
        #The runAwareRateSIR function uses foreach and %dopar%. Make sure you have set up your parallel 
#environment correctly (e.g., doParallel, doSNOW, etc.). If you prefer a single-core approach, replace %dopar% with %do% or a standard for loop.

    #Thresholding:
        #The final result in runAwareRateSIR considers only those runs where the epidemic size exceeded 0.001 * vcount(...) (or 0.01 * vcount(...) in the earlier commented code). This is to avoid counting trivial outbreaks as part of the average.

    #Modifying State Definitions:
       # Feel free to customize the numeric encoding of state or to add additional states. However, ensure you update all relevant logic in awareRateSIR consistently.

   # Large Networks:
        #Be aware that for large networks, the repeated calls to neighborhood(...) can be expensive. You may need efficient data structures or adjacency lists to speed up computations.

#--------------------------------------------------------------------------------------#
source('.../awareRateSIR.R')
source('.../InitialImmunization.R')

runClassicalgo<- function (gcontact, gcom, VertexCnt, interval, start, finish,
                           measure1, iter, maxstep, Beta,
                           mu, fi, gama, alpha){
  
  # state = rep( 0, VertexCnt );
  # awareness = rep(-1, VertexCnt);
  
  lennum <- ((finish - start ) / interval);
  lennum <- lennum + 1;
  res2 <- rep(0, lennum );
  
  numseed <- floor(0.02 * VertexCnt)
  
  numseq <- seq(start, finish, interval)
  #numImun <- floor( interval * VertexCnt )
  numImun <- floor(numseq* VertexCnt);
  
  for ( i in (1:length(numseq)) ){
    
    
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
    
    deleted2 <- classicImmunization(Upgcontact, measure1, numImun[i] )

    delname2 <- rep(0, length(deleted2))
    for(j in 1:length(deleted2)){
      delname2[j] <- get.vertex.attribute(Upgcontact, 'nameid', deleted2[j])

    }
    Upgcontact<- Upgcontact-deleted2;
    Upgcom <- Upgcom-deleted2
    
    state[delname2] = 3
    res2[i] <- runAwareRateSIR(Upgcontact,Upgcom, state, awareness, maxstep, Beta, mu,
                               fi, gama, alpha, iter, numseed)
    cat("round:")
    cat(i)
    cat(" ", res2[i])
    cat("\n")
  }
  
  print("epidemic size after classic immunization of population is:")
  
 # write.xlsx(res2,"D:/afterclassic.xlsx",sheetName = "Newdata")
  return(res2)
  
}
