include ../../make.inc
MW=../../mwrap

all:
	$(MW) -mex fbmex -m foobar.m foobar.mw
	$(MW) -mex fbmex -c fbmex.c foobar.mw
	$(MEX) fbmex.c

clean:
	rm -f foobar.m fbmex.c fbmex.mex* *.o* *~
