This folder includes some notebooks and a python script used in the Magritte paper in Astronomy and Computing.
The python script generates the illustrations for the adaptive Ng-acceleration section.
The jupyter notebooks are used for the illustrations in the grid remeshing section.

We expect the notebooks to be run in the following order:
- import_and_reduce_phantom
- compute_errors
- Plot_relative_differences

To show the speed of the remeshing algorithm, a Phantom model was remeshed and run using both GMSH and the recursive implementation.
This also includes a remesher based on Haar wavelets, but that implementation did not perform adequately, this did not make it into the paper.
