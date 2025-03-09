#This code provides an example of how to use functions defined in separate files (such as InitialImmunization.R, awareRateSIR.R, runawarealgo.R, runClassicalgo.R).
# It first reads data for a physical contact graph (with a power-law structure) from a CSV file, then adds a number of random edges to create a virtual layer. 
#Next, it configures epidemic simulation parameters (infection rate, recovery rate, awareness parameters, etc.) 
#and calls functions such as runClassicalgo to calculate the epidemic size under various immunization strategies. You can also use runawarealgo in a similar way. 
#This example demonstrates how multiple source files can be loaded to build both physical and virtual network layers and ultimately perform immunization and SIR simulations.
######################################################

source('.../InitialImmunization.R')
source('..../awareRateSIR.R')
source('.../runawarealgo.R')
source('.../runClassicalgo.R')
 install.packages("igraph")
 install.packages("xlsx")
 install.packages("foreach")
 install.packages("doParallel")
library(igraph)
library(xlsx)
library(foreach)
library(doParallel)

main <- function(){
  no_cores <- detectCores() - 2
  registerDoParallel(no_cores)
  ###########################

  dat = read.csv(".../powerlaw.csv" , header = FALSE)
  el = as.matrix(dat)
  gcontact = graph.edgelist(el, directed = FALSE)

   gcomtmp = graph(edges = NULL, n=vcount(gcontact), directed = FALSE)
   gcom = union(gcomtmp, gcontact)
 
   addnum = 2 * floor(0.4*vcount(gcontact))
   randomedges <- floor(runif( addnum, 1, vcount(gcontact)+1 ) ) 
   gcom <- add.edges(gcom, randomedges) 
   b <- as_ids(V(gcom))
   tmpnum = ecount(gcom) - ecount(gcontact)
   if( tmpnum <= (addnum / 2)){
    while(TRUE){
       if( tmpnum < (addnum / 2)){
         randomedges <- floor(runif(2, 1, vcount(gcontact)+1 ) )
        gcom <- add.edges(gcom, randomedges)
         tmpnum = ecount(gcom) - ecount(gcontact)
       }
       else
      {
        break;
        }
     }
   }

  # Gcom=  get.edgelist(gcom, names=FALSE)
  # Gcom= as_edgelist(gcom, names = TRUE)

  ##########################initializtion#########################
  #iter=20
  maxstep=100;  iter=10;
  fi=0.8; gama=0.1; alpha=0.3;
  Beta = 0.2; mu = 0.1;
  #Beta = 0.9; mu = 0.2;
  
  measure1 = 5;
  Omega = 0.8;#physical layer constant
  Eta = 0.2;#communication layer constant
  #interval <- 0.10; start <- 0.0; finish <- 0.10;
  #interval <- 0.04; start <- 0.02; finish <- 0.78;
  #interval <- 0.005; start <- 0.005; finish <- 0.1;
  interval <- 0.01; start <- 0.01; finish <- 0.10;
  #interval <- 0.02; start <- 0.10; finish <- 0.2;
  lennum <- ((finish - start ) / interval);
  lennum <- lennum + 1;
  res <- rep(0, lennum );
  res2 <- res
  res3 <- res
  cnt <- min ( vcount(gcontact), vcount(gcom));
  bool = TRUE;
  ################################################################
  #create two vector containing constant values of Eta and Omega values
  #OmegaArray = c( 0.8, 0.8, 0.8, 0.8)
  #EtaArray = c(  0.2, 0.2, 0.2, 0.2)
  #fiArray = c(0.8);
  # alphaArray = c ( 0.3, 0.5, 0.7 , 0.9 )
  alphaArray = c ( 0.3)
  #gamaArray = c ( 0.1, 0.1, 0.1, 0.1 )
  #alphaArray = c (0.5,0.7,0.9)
  #numDiag = length(EtaArray)
  numDiag = length(alphaArray)


  # for ( i in 1:numDiag){
    #Omega <- OmegaArray[i]
    #Eta <- EtaArray[i]
    #fi <- fiArray[i]
    #gama <- gamaArray[i]
    # alpha <- alphaArray[i]

    # print(alpha)

    # newlist <- runawarealgo(gcontact, gcom, cnt,interval, start, finish,
    #                          iter, maxstep,
    #                         Beta, mu, fi, gama, alpha, bool)
    # #newlist <-list("bool"=bool, "result"=res)
    # bool = newlist$bool
    # res = newlist$result
  # }
  
  
  
   
  for ( d in 1:5){
    measure1 = d; 
  # bool = FALSE;
  # 
  #   if (d==1)
  #     bool=TRUE;

  res2 <- runClassicalgo(gcontact, gcom, cnt, interval, start, finish,
                         measure1, iter, maxstep, Beta,
                         mu, fi, gama, alpha)
  }
  # alpha=0.3;
  # print(alpha)
  
  # res2 <- runClassicalgo(gcontact, gcom, cnt, interval, start, finish,
  #                        measure1, iter, maxstep, Beta,
  #                        mu, fi, gama, alpha)
  # 
  # for(i in 1:numDiag){
  #   diff = (res2 - res) / cnt
  #   print("difference is")
  #   print(diff)
  #   print( max(diff))
  #   print(min(diff))  
  #   }

  
  print(res)
  # print(res2)
  
}
