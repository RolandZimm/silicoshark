# silicoshark

Here, we provide the information required for the usage of the general/shark tooth model
We provide the following files:

1) toothmodel.f90   ! the sourcecode of the tooth model
2) toothmodel.e    ! a compiled version of the model
3) make_mesh_for_gnuplot.f90   ! the sourcecode for a program that transforms the primary output of the simulations into plottable and visualizable data.
4) MMFG.e    !  a compiled version of the latter
5) multicusp.ini   ! initial conditions for the toot model

Some file might be compressed (.gz), and can be extracted in the command line by gzip -d *.gz

The model (1) is derived from a previously published version of ToothMaker:
Salazar-Ciudad I, Jernvall J, A computational model of teeth and the developmental origins of morphological variation. Nature 464, 583 (2010).
http://dead.cthulhu.fi/ToothMaker
In order to be used, the sourcecode has to be compiled. We used the freely available gfortran program for this purpose
The most simple respective compilation command (to be executed in the command line) is: 

gfortran toothmodel.f90 -o toothmodel.e

To run (2), the executable program, you have to provide an input file (containing the developmental parameters that define a specific in silico tooth) as well as the total number of iterations.
For instance, using (5), such a command is:

./toothmodel.e multicusp.ini 12000

This command will, notwithstanding potential compilation errors or inviable parameter combinations, give rise to several output files (containing information about cell positions and neighborhoods at every 1000 iterations), e.g. multicusp.ini.out1000, multiscusp.ini.out2000 etc.
Since these output files cannot be plotted directly, we also provide a short program (3) and its precompiled version (4) that will translate simulation outputs into plottable meshes.
To generate (4) from (3) : gfortran make_mesh_for_gnuplot.f90 -o MMFG.e
To execute the program: ./MMFG.e multiscusp.ini.out10000, which will generate multiscusp.ini.out10000_triangle
This file will contain the following information, arranged in columns:
triangle number : cell number : x : y : z : gene product concentrations. Column 6 represents [Fgf], column 7 [Bmp] and column 10 [Shh].
This can be plotted, for instance using gnuplot. If installed, use the 3rd and 5th column to plot the projection view used throughout the manuscript. The "with line" argument allows to connect edges between neighboring cells which will emulate an impression of surface.

The following files contain x:y-points delineating outlines of 3-cuspid and 5-cuspid target shapes, as used in the paper:
6)   pentac_target.txt
7)   tric_target.txt

These files are extracted from CT-scanned sharkteeth previously part of the following thesis:
Berio F, Multiscale variation of 3D tooth forms in selachians and developmental and evolutionary inferences: Odyssey of a scyliorhinid tooth. Populations and Evolution [q-bio.PE]. Université de Lyon (2021).
