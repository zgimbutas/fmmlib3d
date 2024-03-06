include ../../../make.inc
MW=../../../mwrap

mex: wrap 
	$(MEX) -I../src femex.cc ../src/*.o
	mv femex.* *.m ../mwfem

wrap: 
	$(MW) -mex femex -c femex.cc \
		-mb assembler.mw mesh.mw elements.mw 

clean:
	rm -f *~
