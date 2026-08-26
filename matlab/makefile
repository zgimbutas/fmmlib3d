all: linux clean

linux: linux-octave-64bit linux-matlab-64bit

macosx: macosx-matlab-64bit macosx-matlab-arm64

windows: windows-matlab-64bit

linux-octave-64bit:
	$(MAKE) -f makefile.mwrap TARGET=octave-linux-openmp clean
	$(MAKE) -f makefile.mwrap TARGET=octave-linux-openmp -j

linux-matlab-64bit:
	$(MAKE) -f makefile.mwrap TARGET=matlab-linux-a64-openmp clean
	$(MAKE) -f makefile.mwrap TARGET=matlab-linux-a64-openmp -j

macosx-octave-64bit:
	$(MAKE) -f makefile.mwrap TARGET=octave-macos-openmp clean
	$(MAKE) -f makefile.mwrap TARGET=octave-macos-openmp -j

macosx-matlab-64bit:
	$(MAKE) -f makefile.mwrap TARGET=matlab-maci64-openmp clean
	$(MAKE) -f makefile.mwrap TARGET=matlab-maci64-openmp -j

macosx-matlab-arm64:
	$(MAKE) -f makefile.mwrap TARGET=matlab-maca64-openmp clean
	$(MAKE) -f makefile.mwrap TARGET=matlab-maca64-openmp -j

windows-matlab-64bit:
	$(MAKE) -f makefile.mwrap TARGET=matlab-windows-w64-cross-openmp clean 
	$(MAKE) -f makefile.mwrap TARGET=matlab-windows-w64-cross-openmp -j

clean:
	rm -f *.o

distclean:
	rm -f *.o

