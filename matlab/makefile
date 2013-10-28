
all: linux clean

linux: linux-octave-64bit linux-matlab-64bit

macosx: macosx-matlab-64bit

windows: windows-matlab-32bit windows-matlab-64bit

linux-octave-64bit:
	$(MAKE) -f makefile.mwrap TARGET=octave-linux-openmp clean
	$(MAKE) -f makefile.mwrap TARGET=octave-linux-openmp -j

linux-matlab-64bit:
	$(MAKE) -f makefile.mwrap TARGET=matlab-linux-a64-openmp clean
	$(MAKE) -f makefile.mwrap TARGET=matlab-linux-a64-openmp -j

macosx-matlab-32bit:
	$(MAKE) -f makefile.mwrap TARGET=matlab-maci-openmp clean
	$(MAKE) -f makefile.mwrap TARGET=matlab-maci-openmp -j

macosx-matlab-64bit:
	$(MAKE) -f makefile.mwrap TARGET=matlab-maci64-openmp clean
	$(MAKE) -f makefile.mwrap TARGET=matlab-maci64-openmp -j

windows-matlab-32bit:
	$(MAKE) -f makefile.mwrap TARGET=matlab-windows-w32 clean
	$(MAKE) -f makefile.mwrap TARGET=matlab-windows-w32 -j
	$(MAKE) -f makefile.mwrap TARGET=matlab-windows-w32-openmp clean 
	$(MAKE) -f makefile.mwrap TARGET=matlab-windows-w32-openmp -j

windows-matlab-64bit:
	$(MAKE) -f makefile.mwrap TARGET=matlab-windows-w64 clean
	$(MAKE) -f makefile.mwrap TARGET=matlab-windows-w64 -j
	$(MAKE) -f makefile.mwrap TARGET=matlab-windows-w64-openmp clean 
	$(MAKE) -f makefile.mwrap TARGET=matlab-windows-w64-openmp -j

clean:
	rm -f *.o

distclean:
	rm -f *.o

