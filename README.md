# Introduction
This repository holds programs written to benchmark the performance of the
algorithms in the [Interface Reconstruction Library (IRL)](https://gitlab.com/robertchiodi/interfacereconstructionlibrary). For now, this mainly consists of a comparison between the half-edge based polyhedron intersection routine in IRL to the open-source packages [VOFTools](http://www.dimf.upct.es/personal/lrj/voftools.html) and [R3D](https://github.com/devonmpowell/r3d).

# Building the Benchmarks
Building the benchmarks is done using CMake (minimum version 3.11). To build the benchmarks, you will first need to update the `config.cmake` file with information for your own machine. This mostly involves providing the paths to IRL ([available here](https://gitlab.com/robertchiodi/interfacereconstructionlibrary)), VOFTools ([available here](http://www.dimf.upct.es/personal/lrj/voftools.html)), and R3D ([available here](https://github.com/devonmpowell/r3d)), with each package already compiled and installed. Please note, VOFTools mut be built with the make target `fortran_3d`.

If you would like to change the compilers or compiler flags, this can be achieved towards the end of `config.cmake`. The default supplied compilers are from the GNU family `g++`,`gcc`, and `gfortran`.

Once `config.cmake` is properly modified, the benchmarks can be built by executing from the root directory the command
`mkdir build && cd build && cmake -C ../config.cmake .. && make`
which will perform an out of source build and place the executable `timing_comp` in the root directory.

# Running the Benchmarks
The executable `timing_comp` expects four command-line arguments (as integers) to be supplied to it. They are (in this order):

 1. The type of results to be generated, chosen by an integer in the range [0,2]

	0. Create a sample file with randomly generated planes that can be plotted with the python script `sample_example.py`.
	1. Run the randomly generated sets of plane intersecting polyhedron tests, which will write out the files `irl_timing.txt`, `voftools_timing.txt`, and `r3d_timing.txt`. The results can be aggregated and presented in a more understandable form by running the python script `postprocess.py` after performing the tests.
	2. Run the volume distribution tests for IRL, which will write the file `distribution_timing.txt` with the average number of cells entered, conservation error, and total time.

2. The number of trials to run (must be >=1000)
3. The max number of planes to test for in the plane intersecting polyhedron tests. Has no effect for case options 0 or 2. (must be >=1)
4. Whether to produce section timings (0), total timings (1), or both (2). Note: This only has an effect if the first input on CLI is 1
