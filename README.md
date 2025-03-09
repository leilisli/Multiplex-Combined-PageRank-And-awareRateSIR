# Multiplex-Combined-PageRank-And-awareRateSIR
'MCPR.R' is a R function to compute the "Multiplex Combined PageRank" and 'awareRateSIR.R' is  a simulation of the SIR disease spread model  (Susceptible, Infected, Recovered groups) along with the transfer of information/awareness among individuals (via a separate communication network).
This program is distributed ny the authors in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
Programming Language Used: R	
Database Used: The Powerlaw network is a scale-free network with a power-law distribution and an exponent parameter γ=2.5\gamma = 2.5γ=2.5. In real-world networks, the power-law exponent typically ranges between 2 and 3 and is often close to 2.4. Therefore, in many studies aiming to simulate real-world conditions, γ=2.5\gamma = 2.5γ=2.5 is assumed. The physical interaction layer is constructed with these characteristics. To build the virtual interaction network, for 0.4 of the nodes in the network, we randomly add edges between random nodes.
The address of this dataset is as follows: Codes\dataset\powerlaw.csv
Installation Notes: After installing R, you must first install and load the igraph library in R install.packages("igraph") /  library(igraph)
To use any of these files, you should use the source command: source(directory)
Make sure to update the path inside read.csv according to where the data file is located.
