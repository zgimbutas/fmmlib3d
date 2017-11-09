# INSTALLATION NOTES FOR FMMLIB3D

3/20/17

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


To prevent segfaults in fortran executables with openmp, you will want to do the following in your shell:

```
export OMP_STACKSIZE=4096M     # omp default stack size is too small
ulimit -c 0                    # no core dumps
ulimit -s unlimited            # unlimited stack
```

This is because large static allocations are used in the fortran test drivers.

Please type `make` to see a list of other make options.
