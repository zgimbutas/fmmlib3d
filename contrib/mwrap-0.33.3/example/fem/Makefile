all:
	mkdir -p mwfem
	(cd src; make)
	(cd interface; make)

clean:
	(cd src; make clean)
	(cd interface; make clean)
	rm -rf mwfem

realclean: clean
	(cd src; make realclean)
