
make -f makefile.mwrap -j4 TARGET=matlab-windows-w64-openmp clean
make -f makefile.mwrap -j4 TARGET=matlab-windows-w64-openmp 

make -f makefile.mwrap -j4 TARGET=matlab-linux-a64-openmp clean
make -f makefile.mwrap -j4 TARGET=matlab-linux-a64-openmp  

make -f makefile.mwrap -j4 TARGET=octave-linux-openmp clean
make -f makefile.mwrap -j4 TARGET=octave-linux-openmp 

make -f makefile.mwrap clean distclean

