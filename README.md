# README #

### The program for finding the optimal path in  complex free energy landscape (for example oxygen pathways in protein). ###

### The workflow: ###
1. A starting point for wavefront (e.g. an active site of a protein) is chosen and Eikonal equation is solved using finite-differences method. This procedure calculates 3D map of wavefront arrival times.
2. Area where wavefont exits the protein is located and used as a starting point for gradient descent algorithm to find the optimal path back to the 
source of the wave.  
3. The cost of the path and the free energy profile is computed by projecting the free energy map on the path.

###Input:###
 1. 3D slowness map of the media in DX format. 
 2. User configurable parameters: "eikonal.in".

### Running the program:###
Run the executable "eikonal" from the folder with the configuration file "eikonal.in"

###Output:###
1. Wavefront arrival times: times.dx
2. The optimal path: path.xyz
3. Free energy profile along the path: pmf.dat

### References: ###
1. T. Zaraiskaya, S. Vassiliev, D. Bruce, Discovering oxygen channel topology in photosystem II using implicit ligand sampling and wavefront propagation, J. Comput. Sci. 5 (2014) 549–555. doi:10.1016/j.jocs.2013.10.001.
2. S. Vassiliev, T. Zaraiskaya, D. Bruce, Molecular dynamics simulations reveal highly permeable oxygen exit channels shared with water uptake channels in photosystem II, Biochim. Biophys. Acta BBA - Bioenerg. 1827 (2013) 1148–1155. doi:10.1016/j.bbabio.2013.06.008.