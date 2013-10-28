include ../make.inc


# === Primary targets ===

../mwrap: mwrap.o lex.yy.o mwrap-ast.o mwrap-typecheck.o \
	mwrap-mgen.o mwrap-cgen.o mwrap-ast.h
	$(CXX) -o ../mwrap mwrap.o mwrap-ast.o \
		mwrap-typecheck.o mwrap-mgen.o mwrap-cgen.o \
		lex.yy.o

mwrap.o: mwrap.cc lex.yy.c mwrap-ast.h
	$(CXX) -c mwrap.cc

mwrap.cc: mwrap.y
	$(BISON) -d -v mwrap.y -o mwrap.cc

lex.yy.o: lex.yy.c
	$(CC) -c lex.yy.c

lex.yy.c: mwrap.l
	$(FLEX) mwrap.l

mwrap-ast.o: mwrap-ast.cc mwrap-ast.h 
	$(CXX) -c -g mwrap-ast.cc

mwrap-typecheck.o: mwrap-typecheck.cc mwrap-ast.h
	$(CXX) -c -g mwrap-typecheck.cc

mwrap-cgen.o: mwrap-cgen.cc mwrap-ast.h mwrap-support.h
	$(CXX) -c -g mwrap-cgen.cc

mwrap-mgen.o: mwrap-mgen.cc mwrap-ast.h
	$(CXX) -c -g mwrap-mgen.cc

mwrap-support.h: mwrap-support.c stringify
	./stringify mex_header < mwrap-support.c > mwrap-support.h

stringify: stringify.c
	gcc -o stringify stringify.c


# === Clean-up targets ===

clean:
	rm -f mwrap.output
	rm -f *.o *~
	rm -f stringify

realclean: clean
	rm -f lex.yy.c mwrap.cc mwrap.hh mwrap-support.h mwrap.pdf

