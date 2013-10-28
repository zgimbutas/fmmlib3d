
hfmm3dpart:
	make -f hfmm3dpart.make -j4

lfmm3dpart:
	make -f lfmm3dpart.make -j4

hfmm3dtria:
	make -f hfmm3dtria.make -j4

lfmm3dtria:
	make -f lfmm3dtria.make -j4

hfmm3dpart-openmp:
	make -f hfmm3dpart.make HOST=linux-gfortran-openmp -j4

cfmm3dpart-openmp:
	make -f lfmm3dpart.make HOST=linux-gfortran-openmp -j4

hfmm3dtria-openmp:
	make -f hfmm3dtria.make HOST=linux-gfortran-openmp -j4

cfmm3dtria-openmp:
	make -f lfmm3dtria.make HOST=linux-gfortran-openmp -j4

clean: 
	rm -f *.o
 
distclean: 
	rm -f int2 fort.* *.o

