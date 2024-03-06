include ../../make.inc
MW=../../mwrap

all: pmex cmex hmex

pmex:
	$(MW) -mex eventqpmex -c eventqpmex.cc -mb eventq_plain.mw
	$(MEX) eventqpmex.cc

hmex:
	$(MW) -mex eventqhmex -c eventqhmex.cc -mb eventq_handle.mw
	$(MEX) eventqhmex.cc

cmex:
	mkdir -p \@eventq
	$(MW) -mex eventqcmex -c eventqcmex.cc -mb eventq_class.mw
	$(MEX) eventqcmex.cc

clean:
	rm -rf \@eventq/
	rm -f EventQ_*.m eventqh.m
	rm -f eventqcmex.* eventqpmex.* eventqhmex.* *.o *~ 
