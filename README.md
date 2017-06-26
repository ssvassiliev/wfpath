# README #

### Programs for finding optimal paths in complex 3D free energy fields ###


### Workflow ###

1. A stating point (e.g. active site of a protein) is chosen and Eikonal equation is solved using finite-differences method. This procedure calculates 3D map of wavefront arrival times.
2. Area where wavefont exits the protein is located and used as a starting point for gradient descent algorithm to find the optimal path back to the source of the wave.  
3. The cost of the path and the free energy profile is computed by projecting the free energy map on the path.