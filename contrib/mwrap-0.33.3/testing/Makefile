include ../make.inc

all: test_transfers test_cpp_complex $(TESTC99COMPLEX) test_syntax \
	test_typecheck test_catch test_fortran1 test_fortran2 \
	test_redirect test_include

test_transfers: 
	../mwrap -mex test_transfersmex \
		-c test_transfersmex.cc \
		-m test_transfers.m test_transfers.mw
	$(MEX) test_transfersmex.cc

test_cpp_complex: 
	../mwrap -cppcomplex \
		-mex test_cpp_complexmex \
		-c test_cpp_complexmex.cc \
		-m test_cpp_complex.m test_cpp_complex.mw
	$(MEX) test_cpp_complexmex.cc

test_c99_complex: 
	../mwrap -c99complex \
		-mex test_c99_complexmex \
		-c test_c99_complexmex.c \
		-m test_c99_complex.m test_c99_complex.mw
	$(MEX) test_c99_complexmex.c

test_syntax:
	- ../mwrap -cppcomplex test_syntax.mw 2> test_syntax.log
	diff test_syntax.log test_syntax.ref

test_typecheck:
	- ../mwrap -cppcomplex test_typecheck.mw 2> test_typecheck.log 
	diff test_typecheck.log test_typecheck.ref

test_catch:
	../mwrap -mex test_catchmex -catch \
		-c test_catchmex.cc \
		-m test_catch.m \
		test_catch.mw
	$(MEX) test_catchmex.cc

test_fortran1:
	../mwrap -mex test_fortran1mex \
		-c test_fortran1mex.cc \
		-m test_fortran1.m \
		test_fortran1.mw
	$(MEX) test_fortran1mex.cc

test_fortran2:
	../mwrap -mex test_fortran2mex \
		-c test_fortran2mex.c \
		-m test_fortran2.m \
		test_fortran2.mw
	$(MEX) test_fortran2mex.c

test_redirect:
	../mwrap -mb test_redirect.mw 

test_include: 
	../mwrap -mex test_includemex\
		-c test_includemex.cc \
		-m test_include.m test_include.mw
	$(MEX) test_includemex.cc

clean:
	rm -f *~ *.mex* *.o* test_typecheck.log
	rm -f test_fortran1.m test_fortran2.m test_transfers.m test_catch.m
	rm -f test_cpp_complex.m test_c99_complex.m
	rm -f test_transfersmex.cc test_catchmex.cc test_fortranmex.cc
	rm -f test_cpp_complexmex.cc test_c99_complexmex.c 
	rm -f test_redirect.m test_redirect1.m
