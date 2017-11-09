# Helmholtz and Laplace particle FMMs in R^2

PROJECT_LIB=fmm3d


###HOST=windows
HOST=linux-gfortran
###HOST=macosx


ifeq ($(HOST),linux)

OBJSUF=o
MODSUF=mod
FC=f77 -c 
FFLAGS=-fast
FLINK=f77 -o $(PROJECT)

else

ifeq ($(HOST),macosx)

OBJSUF=o
MODSUF=mod
FC=gfortran -c 
FFLAGS=-O2
FLINK=gfortran -o $(PROJECT)

else

ifeq ($(HOST),linux-gfortran)

OBJSUF=o
MODSUF=mod
FC=gfortran -c 
FFLAGS=-O2
FLINK=gfortran -o $(PROJECT)

else

ifeq ($(HOST),linux-gfortran-openmp)

# buggy compiler, proceed anyway, use gfortran > 4.4.0
OBJSUF=o
MODSUF=mod
FC=gfortran -c 
FFLAGS=-O3 --openmp
FLINK=gfortran -o $(PROJECT) --openmp
### export OMP_NUM_THREADS=4
### export OMP_STACKSIZE=1024M

else

ifeq ($(HOST),linux-intel-openmp)

OBJSUF=o
MODSUF=mod
FC=ief77 -c
FFLAGS= -xT -O3 -ip -openmp 
FLINK=ief77 -o $(PROJECT) -xT -O3 -ip -openmp -static
WITH_SECOND=1

else 

ifeq ($(HOST),linux-intel)

OBJSUF=o
MODSUF=mod
FC=ief77 -c
FFLAGS=-xT -O3 -ip
FLINK=ief77 -o $(PROJECT) -xT -O3 -ip
WITH_SECOND=1

else 

ifeq ($(HOST),linux-intel-profile)

OBJSUF=o
MODSUF=mod
FC=ief77 -c
FFLAGS=-xT -O3 -ip -p
FLINK=ief77 -o $(PROJECT) -xT -O3 -ip -p 
WITH_SECOND=1

else 

ifeq ($(HOST),linux-intel64)

OBJSUF=o
MODSUF=mod
FC=ief77 -c
FFLAGS=-xT -O3 -ip -mcmodel medium -i-dynamic 
FLINK=ief77 -o $(PROJECT) -xT -O3 -ip -mcmodel medium -i-dynamic 
###export LD_LIBRARY_PATH=/opt/intel/fce/10.1.008/lib
WITH_SECOND=1

else 

ifeq ($(HOST),linux-lahey64)

OBJSUF=o
MODSUF=mod
FC=lf77_x64 -c
FFLAGS=-fast --model medium
FLINK=lf77_x64 -o $(PROJECT) -fast --model medium /opt/lf6481/lib64/libelf.so.0
###FLINK=lf77_x64 -o $(PROJECT) -fast --model medium 
WITH_SECOND=1

else 
     
ifeq ($(HOST),linux-fort77)

OBJSUF=o
MODSUF=mod
FC=fort77 -c
FFLAGS=-fast 
FLINK=fort77 -o $(PROJECT) -static -g
WITH_SECOND=1

else 

ifeq ($(HOST),linux-profile)

OBJSUF=o
MODSUF=mod
FC=fort77-i386-m32 -c
FFLAGS=-fast -profile -pg
FLINK=fort77-i386-m32 -o $(PROJECT) -static -pg -g
WITH_SECOND=1

else 

OBJSUF=obj
MODSUF=mod
FC=lf95 -c
FFLAGS=-O1
FLINK=lf95 -out $(PROJECT)

endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif


.PHONY: $(PROJECT) clean list

.f.$(OBJSUF):
	$(FC) $(FFLAGS) $<

.f.$(MODSUF):
	$(FC) $(FFLAGS) $<

.SUFFIXES: $(MODSUF) .$(OBJSUF) .f .c

# SOURCE FILE LIST
#
vpath %.f .

FMODS = 

FSRCS = cdjseval3d.f helmrouts3d.f laprouts3d.f projections.f		\
        trilib.f d3mtreeplot.f hfmm3dpart.f legeexps.f rotproj.f	\
        triquadflatlib.f d3tstrcr.f hfmm3drouts.f lfmm3dpart.f		\
        rotviarecur3.f yrecursion.f dfft.f hfmm3dtria.f lfmm3drouts.f	\
        h3dterms.f l3dterms.f lfmm3dtria.f triagauc.f h3dtrans.f	\
        l3dtrans.f prini.f triahquad.f h3dtrirouts.f l3dtrirouts.f	\
        prinm.f triasymq.f

ifeq ($(WITH_SECOND),1) 
FSRCS += second-r8.f
endif

#
# object files list
MODS    =  $(FMODS:.f=.$(MODSUF)) 
OBJS    =  $(FMODS:.f=.$(OBJSUF)) $(FSRCS:.f=.$(OBJSUF)) 
#

lib$(PROJECT_LIB).a:   $(MODS)   $(OBJS)
	ar rvs lib$(PROJECT_LIB).a $(OBJS)

#
clean: 
	rm -f $(OBJS)
# 
distclean: 
	rm -f $(OBJS) int2 *.a fort.*
# 
list: $(FSRCS)
	echo $^
#
pack: $(FSRCS)
	cat $^ > _tmp_.pack.f
#
