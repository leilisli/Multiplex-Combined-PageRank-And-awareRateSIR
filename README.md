# Multiplex-Combined-PageRank-And-awareRateSIR

## Overview

This repository provides R implementations for analyzing epidemic spreading and awareness diffusion processes in multiplex networks.

The project includes:

- **`MCPR.R`**: A function for computing the *Multiplex Combined PageRank (MCPR)* centrality measure.
- **`awareRateSIR.R`**: A simulation of the *Susceptible–Infected–Recovered (SIR)* epidemic model, incorporating awareness diffusion through a separate communication layer.

This framework is designed for research in complex networks, multilayer systems, epidemic modeling, and information spreading dynamics.

---
## Disclaimer
This program is distributed by the author in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

## Programming Language Used: 
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
```


To use any of these files, you should use the source command:  
source(directory)

Make sure to update the path inside read.csv according to where the data file is located.
