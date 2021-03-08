# INSTALLATION NOTES FOR FMMLIB3D

3/7/21

### Requirements:

The basic requirements are a fortran compiler such as gfortran.

If you want to build mex interfaces to octave, you will need octave and its
development libraries.

David Bindel's [MWrap](http://www.cs.cornell.edu/~bindel/sw/mwrap/) is shipped in this package (in the `contrib` directory), including the executable
(`bin/mwrap`). If you want to recompile mwrap, you will also need bison and flex.

In ubuntu linux the complete set of requirements can be installed via

`sudo apt-get install bison flex octave liboctave-dev`

In fedora/centos/EL linux you instead need

`sudo yum install bison flex octave octave-devel`



### Compilation on linux:

`make test`

for a single-threaded test, or

`make test-openmp`

When requested to `ENTER n`, a good number to enter is 10000 (do not format this
as 1e4, since it can only read integer formats).

Errors in outputs close to machine precision should be reported. 
Warnings about floating-point exceptions are normal and to be ignored.


To prevent segfaults in fortran executables with openmp, you will want to do the following in your shell:

```
export OMP_STACKSIZE=4096M     # omp default stack size is too small
ulimit -c 0                    # no core dumps
ulimit -s unlimited            # unlimited stack
```

This is because large static allocations are used in the fortran test drivers.

Please type `make` to see a list of other make options.


**Note** if using gfortran v 10 or above: Since we use passing of size-1 arrays
as pointers, GCC10+ raises errors. You will need to add
```
FFLAGS+=-std=legacy
```
in the relevant sections of `src/Makefile` and `examples/*.make`
which turns these into mere warnings.


### Octave notes

These routines are identical to the MATLAB ones.
To test the shipped mex binary, use eg
`cd matlab; octave test_hfmm3dpart_direct.m`.
To rebuild the binary use `make mex-octave`, then test as above.

