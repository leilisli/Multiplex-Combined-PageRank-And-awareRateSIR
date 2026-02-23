# Multiplex-Combined-PageRank-And-awareRateSIR

**MCPR.R** is a function in R to calculate the "Multiplex Combined PageRank" and 'awareRateSIR.R' is a simulation of the SIR model (Susceptible, Infected, Recovered), along with the transfer of information or awareness amongst individuals through a separate communication network.

This program is distributed by the author in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

##Programming Language Used: 
**R**

## Database Used:

### 1. Synthetic Power-Law Network
The Powerlaw network is a scale-free network with a power-law distribution and an exponent parameter γ = 2.5. In real-world networks, the power-law exponent typically ranges between 2 and 3 and is often close to 2.4. Therefore, in many studies aiming to simulate real-world conditions, γ = 2.5 is assumed. The physical interaction layer is constructed with these characteristics. To build the virtual interaction network, for 0.4 of the nodes in the network, we randomly add edges between random nodes.
The address of this dataset is as follows:  
DataSet/powerlaw.csv

### 2. Copenhagen Network Dataset (Real-World Data)
In addition to the synthetic Powerlaw network, this project also uses the Copenhagen Network dataset, which contains real-world human interaction and communication data collected in a university environment. This dataset enables realistic modeling of epidemic spreading and awareness diffusion processes.

The address of the Copenhagen dataset is as follows:  
https://networks.skewed.de/net/copenhagen

## Notebook File
**copenhagen_network_study.ipynb**
This Jupyter Notebook provides additional analysis and exploration of the Copenhagen Network dataset. It can be used to:

- Reproduce experiments  
- Visualize network properties  
- Investigate epidemic and awareness spreading dynamics on real-world data 


## Installation Notes 
1. Install **R**
2. Install and load the required package:

```r
install.packages("igraph")
library(igraph)

To use any of these files, you should use the source command:  
source(directory)

Make sure to update the path inside read.csv according to where the data file is located.
