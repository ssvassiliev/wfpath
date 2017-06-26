# README #

### The program for finding optimal paths in  complex free energy landscapes. ###

### The workflow: ###
1. A starting point for wavefront (e.g. an active site of a protein) is chosen and Eikonal equation is solved using finite-differences method. This procedure calculates 3D map of wavefront arrival times.
2. Area where wavefont exits the protein is located and used as a starting point for gradient descent algorithm to find the optimal path back to the 
source of the wave.  
3. The cost of the path and the free energy profile is computed by projecting the free energy map on the path.

###Input:###
 1. 3D slowness map of the media in .dx format. 
 2. User configurable parameters: "eikonal.in".

### Running the program:###
Run the executable "eikonal" from the folder with the configuration file "eikonal.in"

###Output:###
1. Wavefront arrival times: times.dx
2. The optimal path: path.xyz
3. Free energy profile along the path: pmf.dat