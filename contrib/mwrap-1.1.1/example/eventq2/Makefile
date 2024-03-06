include ../../make.inc
MW=../../mwrap

mex:
	$(MW) -mex eventq2mex -c eventq2mex.cc -mb eventq2.mw
	$(MEX) eventq2mex.cc

clean:
	rm -f EventQ_*.m
	rm -f eventq2mex.* *.o *~ 
