
#This code implements a simulation of the SIR disease spread model (Susceptible, Infected, Recovered groups) along 
#with the transfer of information/awareness among individuals (via a separate communication network). The functionality of the code is explained in detail below:
#a. Function awareRateSIR
#	Inputs:
#	graphContact: A network representing physical contacts through which the disease can be transmitted.
#	graphCom: A network representing communication links used to spread awareness.
#	state: A vector representing each node's state (typically coded as 0 for susceptible, 1 for aware, 2 for infected, and 3 for recovered).
#	awareness: A vector representing the awareness level of each node.
#	Other parameters include:
#	maxstep: Maximum simulation steps.
#	Beta: Base transmission rate of the disease.
#	mu: Recovery rate.
#	fi: Factor that reduces infection rate based on awareness.
#	gama: Rate at which awareness fades over time.
#	alpha: Probability of transmitting awareness.
#	numseed: Number of initially infected nodes.
#	Main Steps in the Function:
#	Seeding the Infection:
#	A given number (numseed) of nodes are randomly selected to become infected.
#	Their state is set to infected (2) and their awareness is initialized to a specific value (here, -1, denoted by the variable komon).
#	Simulation Loop (for up to maxstep steps or until there are no infected nodes):
#	Awareness Update:
#	Nodes with an awareness level below a threshold gradually increase their awareness level.
#	If a susceptible node (state 1) has an awareness of -1, it is updated to the initial awareness level (1).
#	Recovery Process:
#	After a certain number of steps (when idx > |komon|), each infected node has a chance (with probability mu) to recover (state becomes 3) and its awareness is reset to zero.
#	Disease Transmission:
#	For every infected node, its neighbors in the contact network are checked.
#	If a neighbor is susceptible (state 0) or already aware (state 1), the infection rate is determined:
#	For aware nodes, the infection rate is reduced using the formula: infectionrate=(1−fiawareness)×Beta\text{infectionrate} = \left(1 - fi^{\text{awareness}}\right) \times Betainfectionrate=(1−fiawareness)×Beta
#	For susceptible nodes, the full rate Beta is used.
#	Based on a random draw, if the neighbor becomes infected, its state is set to 2 and awareness is initialized to komon.
#	Information Transmission:
#	Nodes that are aware (awareness ≥ 0) spread information to their neighbors in the communication network.
#	With probability alpha, if a neighbor is susceptible, it becomes aware (state is set to 1) with an initial awareness level.
#	Additionally, if the neighbor already has a higher awareness level, its awareness is decreased, indicating the positive impact of timely information.
#	Awareness Fading:
#	For each aware node, there is a chance (with probability gama) that the awareness level increases numerically (representing a fading or reduction in the influence of awareness over time).
#	Output:
#	The function returns the total epidemic size (i.e., the number of nodes that got infected during the simulation).
________________________________________
#b. Function runAwareRateSIR
#	Purpose:
#	This function runs the simulation multiple times (given by iter) and computes the average epidemic size only for simulations where the epidemic grew beyond a certain threshold relative to the total number of nodes.
#	Implementation Details:
#	Two versions are provided:
#	A sequential version that runs simulations in a loop.
#	A parallel version that uses the foreach package for concurrent simulation runs.
#	After running the simulations, it calculates the mean epidemic size for those runs where the epidemic size exceeded a specified threshold (for example, 0.001 times the total number of nodes).



##############################################################
awareRateSIR <- function ( graphContact, graphCom, state, awareness, maxstep=100,
                           Beta, mu, fi, gama, alpha, numseed){
  epidemicsize = 0
  komon = -1
  vertexcnt = min( vcount( graphContact ), vcount( graphCom ) )
  numinfected = numseed;
  #infectedtmp=c(30,52,40,25)
  for(i in 1:numinfected){
    infected = floor(runif(1)*vertexcnt+1);
   # infected=infectedtmp[i];
    while( length( neighborhood( graphContact, 1, infected[1])[[1]][-1] )  == 0
           | state[infected] ==3 ){
      infected = floor(runif(1)*vertexcnt+1);
    }
    state[infected] = 2;
    awareness [infected] = komon;
    epidemicsize = epidemicsize + 1;
  }
  
  
  for( idx in ( 1:maxstep )){
    if( length(which(state==2)) == 0 )
    {return(epidemicsize)}
    I <- which( state==2 )
    ###################
    I_idx <- which(awareness < -1)
    awareness[I_idx] <- awareness[I_idx] + 1
    I_idx2 <- which(awareness == -1 & state  == 1)
    awareness[I_idx2] <- 1
    ##############
    if( idx>abs(komon) ){# first round
      #recovery
      for( j in 1:length(I)){##
        p = runif(1)
        if(p <= mu){
          state[ I[j] ] = 3
          awareness[I[j]] = 0
        }
      }
    }
    ii <- which( state==2)
    if(length(ii)==0){
      return(epidemicsize)
    }
    #disease transmision
    for( i in 1:length(ii)){
      np <- neighborhood( graphContact, 1, ii[i] )[[1]][-1]
      n <- np[which(state[np]==0 | state[np]==1)]
      if(length(n)>0){
        ptm=proc.time()[3]
        for(k in 1:length(n)){
          if (state[n[k]]== 1){
            si = awareness[n[k]]
            infectionrate = ( 1-( fi^si) )*Beta
          }
          if (state[n[k]]== 0 ){
            infectionrate = Beta
          }
          p2 <- runif(1);
          if( p2 < infectionrate ){
            state[n[k]] = 2
            awareness[n[k]] = komon
            epidemicsize = epidemicsize + 1
          }
        }
        tt=proc.time()[3] - ptm
      }
    }
    #information transmision
    aware <- which( awareness >= 0)
    if( length(aware)>0){
      for( k in 1:length(aware)){
        awareIdx <- awareness[aware[k]]
        na <- neighborhood(graphCom, 1, aware[k])[[1]][-1]
        if(length(na)>0){
          for( l in 1:length(na)){
            ptm2=proc.time()[3]
            p3 <- runif(1);
            if( p3 < alpha ){
              if ( state[na[l]] == 0 ){
                state[na[l]] = 1
                awareness[na[l]] = 1
              }
              else if ( awareness[na[l]] > awareIdx )
              {
                awareness[na[l]] <- awareIdx + 1
              }
            }
          }
          t=proc.time()[3] - ptm2
        }
      }
    }
    #information fading
    for( p in 1:length(aware)){
      p4 <- runif(1)
      if(p4 < gama){
        awareness[aware[p]] <- awareness[aware[p]] + 1
      }
    }
  }
  return (epidemicsize)
}

runAwareRateSIR <- function(graphContact,graphCom ,state ,awareness, maxstep, Beta, mu,
                            fi, gama, alpha, iter, seeds){
  
  tp <- rep(0, iter)
  tpm0<- proc.time()[3];
  for ( i in 1:iter){
    tp[i] <- awareRateSIR(graphContact, graphCom,state, awareness, maxstep, Beta, mu,
                          fi, gama, alpha, seeds)
    t <- proc.time()[3]-tp
    cat("the iteration number :", i)
    cat("finish. output is",tp[i])
    cat("time is", t)
    print("\n")
  }
  vcnt = min( vcount(graphContact), vcount(graphCom))
  tmp <- which( tp > (0.01 * vcnt) )
  if( length(tmp) != 0 )
  {return( mean( tp[which( tp > (0.01 * vcnt) ) ]  ) )}
  else
  {return(0)}
}
runAwareRateSIR <- function(graphContact,graphCom ,state ,awareness, maxstep, Beta, mu,
                            fi, gama, alpha, iter, numseed = 2){
  
  # tp <- rep(0, iter)
  # for ( i in 1:iter){
  #   tp[i] <- awareRateSIR(graphContact, graphCom,state, awareness, maxstep, Beta, mu,
  #                         fi, gama, alpha, numseed)
  # }
  
  tp<-foreach( i=1:iter,  .inorder = FALSE, .export = "awareRateSIR",
               .packages = c("igraph", "foreach") ) %dopar%{
                 tp <- awareRateSIR(graphContact, graphCom,state, awareness, maxstep, Beta, mu,
                                    fi, gama, alpha, numseed)
                 return(tp)
               }
  tp <- unlist(tp)
  vcnt = min( vcount(graphContact), vcount(graphCom))
  tmp <- which( tp > (0.001 * vcnt) )
  if( length(tmp) != 0 )
  {return( mean( tp[which( tp > (0.001 * vcnt) ) ]  ) )}
  else
  {return(0)}
}
