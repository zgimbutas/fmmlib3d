include make.inc

bin:
	(cd src; make)

test: 
	(cd testing; make)

demo:
	(cd example; make)

clean:
	rm -f mwrap
	(cd src; make clean)
	(cd example; make clean)
	(cd testing; make clean)

realclean: clean
	(cd src; make realclean)
