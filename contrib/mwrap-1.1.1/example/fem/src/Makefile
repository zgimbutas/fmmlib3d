# Targets for object files (we compile with MEX to avoid headaches
# about -fPIC and other flag compatibility issues)

include ../../../make.inc

OBJS= \
	assembler.o mesh.o \
	gauss2by2.o quad2d1.o \
	scalar1d.o scalar2d.o elastic2d1.o

all: $(OBJS)

.cc.o:
	$(MEX) -c $*.cc

elastic2d1.cc: 
	matexpr elastic2d.cc > elastic2d1.cc

quad2d1.cc:
	matexpr quad2d.cc > quad2d1.cc

clean:
	rm -f *.o *~

realclean:
	rm -f *.o elastic2d1.cc quad2d1.cc
