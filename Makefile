# master makefile 

default:
	@echo Specify what to make:
	@echo " test - compile and run a simple test"
	@echo " test-openmp - compile and run a simple OpenMP test"
	@echo " lib - fortran libraries"
	@echo " doc - compile documentation"
	@echo " mwrap - build matlab wrapper generator"
	@echo " mex - build octave/matlab mex files"
	@echo " test-mex-octave - test octave mex files"
	@echo " test-mex-matlab - test matlab mex files"
	@echo " test-mwrap-mex-octave: build and test mwrap wrapper generator" 
	@echo " all - (almost) all of the above"
	@echo " clean - remove object files"

SRC = ./src
DOC = ./doc
MATLAB = ./matlab
EXAMPLES = ./examples
MWRAP = ./contrib/mwrap-0.33.3

test: 
	cd $(EXAMPLES); $(MAKE) clean hfmm3dpart

test-openmp: 
	cd $(EXAMPLES); $(MAKE) clean hfmm3dpart-openmp

lib: 
	cd $(SRC); $(MAKE) 

doc: 
	cd $(DOC); $(MAKE) 

mex: mex-octave mex-matlab

mex-octave:
	cd $(MATLAB); $(MAKE) linux-octave-64bit

test-mex-octave: 
	cd $(MATLAB); octave test_hfmm3dpart_direct.m

mex-matlab: 
	cd $(MATLAB); $(MAKE) linux-matlab-64bit

test-mex-matlab: 
	cd $(MATLAB); matlab -nodesktop -nojvm < test_hfmm3dpart_direct.m

test-mwrap-mex-octave: mwrap 
	cd $(MATLAB); $(MAKE) MWRAP=../bin linux-octave-64bit

mwrap: 
	cd $(MWRAP); $(MAKE) 
	cp -f $(MWRAP)/mwrap ./bin

all: lib test doc mwrap mex test-mex-octave


clean:
	cd $(SRC); $(MAKE) clean
	cd $(MATLAB); $(MAKE) clean
	cd $(EXAMPLES); $(MAKE) clean
	cd $(MWRAP); $(MAKE) clean
	cd $(DOC); $(MAKE) clean

distclean:
	cd $(SRC); $(MAKE) distclean
	cd $(MATLAB); $(MAKE) distclean
	cd $(EXAMPLES); $(MAKE) distclean
	cd $(MWRAP); $(MAKE) realclean
	cd $(DOC); $(MAKE) distclean
	rm -f ./bin/mwrap

