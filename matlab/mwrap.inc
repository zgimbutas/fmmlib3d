#
#  Generate Matlab/Octave MEX files
#
#  gfortran > 4.4.0 is strongly recommended
#


###TARGET=octave-linux
###TARGET=octave-linux-openmp

###TARGET=octave-linux-intel
###TARGET=octave-linux-intel-openmp

###TARGET=matlab-linux-a64
###TARGET=matlab-linux-a64-openmp

###TARGET=matlab-linux-intel-a64
###TARGET=matlab-linux-a64-intel-openmp

###TARGET=matlab-maci
###TARGET=matlab-maci-openmp
###TARGET=matlab-maci64
###TARGET=matlab-maci64-openmp

###TARGET=matlab-windows-w32
###TARGET=matlab-windows-w32-openmp
###TARGET=matlab-windows-w64
###TARGET=matlab-windows-w64-openmp



ifeq ($(TARGET), octave-linux)

FFLAGS=-g -O2 -fPIC

# Use the following for GNU Octave + gfortran 4.4.x
TESTC99COMPLEX=test_c99_complex
FC=    gfortran
CC=    gcc
CXX=   g++
MEX=   FFLAGS="-g -O2 -fPIC" mkoctfile --mex
MEX_FC=   $(FC)
MEX_CC=   $(CXX)
MW=mwrap
MWFLAGS=-c99complex -mex	
###MEX_LIBS=-lm
MEX_LIBS=-Wl,--no-undefined -lm -ldl -lgfortran
MEX_LINK=

MEX_FLAGS=
OBJSUF=o
MODSUF=mod

endif


ifeq ($(TARGET), octave-linux-openmp)

FFLAGS=-g -O3 --openmp -fPIC

# Use the following for GNU Octave + gfortran 4.4.x
TESTC99COMPLEX=test_c99_complex
FC=    gfortran
CC=    gcc
CXX=   g++
MEX=   FFLAGS="-g -O3 --openmp -fPIC" mkoctfile --mex 
MEX_FC=   $(FC)
MEX_CC=   $(CXX)
MW=mwrap
MWFLAGS=-c99complex -mex	
###MEX_LIBS=-lgomp -lm
MEX_LIBS=-Wl,--no-undefined -lgomp -lm -ldl -lgfortran
MEX_LINK=

MEX_FLAGS=
OBJSUF=o
MODSUF=mod

endif


ifeq ($(TARGET), octave-linux-intel)

FFLAGS=-g -xT -O3 -ip -fPIC
intel_prefix="/opt/intel/fce/10.1.008/lib"

# Use the following for GNU Octave + gfortran 4.4.x
TESTC99COMPLEX=test_c99_complex
FC=    ief77
CC=    gcc
CXX=   g++
MEX=   FFLAGS="-g -xT -O3 -ip -fPIC" mkoctfile --mex 
MEX_FC=   $(FC)
MEX_CC=   $(CXX)
MW=mwrap
MWFLAGS=-c99complex -mex	
###MEX_LIBS=-lgomp -lf2c -lm
MEX_LIBS=-Wl,--no-undefined -lm -ldl -Wl,--rpath,$(intel_libs)  -L$(intel_libs)  -lifcore -lifcoremt -lsvml -lirc -limf -liomp5 -ldl
MEX_LINK=

MEX_FLAGS=
OBJSUF=o
MODSUF=mod

endif


ifeq ($(TARGET), octave-linux-intel-openmp)

FFLAGS=-g -xT -O3 -ip -openmp -fPIC
intel_prefix="/opt/intel/fce/10.1.008/lib"

# Use the following for GNU Octave + gfortran 4.4.x
TESTC99COMPLEX=test_c99_complex
FC=    ief77
CC=    gcc
CXX=   g++
MEX=   FFLAGS="-g -xT -O3 -ip -openmp -fPIC" mkoctfile --mex 
MEX_FC=   $(FC)
MEX_CC=   $(CXX)
MW=mwrap
MWFLAGS=-c99complex -mex	
###MEX_LIBS=-lgomp -lf2c -lm
MEX_LIBS=-Wl,--no-undefined -lm -ldl -Wl,--rpath,$(intel_libs)  -L$(intel_libs)  -lifcore -lifcoremt -lsvml -lirc -limf -liomp5 -ldl
MEX_LINK=

MEX_FLAGS=
OBJSUF=o
MODSUF=mod

endif



ifeq ($(TARGET), matlab-linux-a64-intel)

FFLAGS=-g -xT -O3 -ip -fPIC
intel_prefix="/opt/intel/fce/10.1.008/lib"

# Use the following for Matlab + gfortran 4.4.x
TESTC99COMPLEX=test_c99_complex
FC=    ief77
CC=    gcc
CXX=   g++
# Uncomment this line for new-style classdef support
# OOFLAG=-DR2008OO
###MEX=   mex $(OOFLAG)
# Use the following for 64-bit MEX
MEX=   mex -largeArrayDims $(OOFLAG)
MEX_FC=   $(FC)
MEX_CC=   $(CXX)
MW=mwrap
MWFLAGS=-c99complex -mex	
MEX_LIBS=-Wl,--rpath,$(intel_libs) -L$(intel_libs)  -lifcore -lsvml -limf -lintlc -lm
MEX_LINK=-output $(MEX_GATEWAY).mexa64

MEX_FLAGS=
OBJSUF=o
MODSUF=mod

endif



ifeq ($(TARGET), matlab-linux-a64-intel-openmp)

FFLAGS=-g -xT -O3 -ip --openmp -fPIC
intel_prefix="/opt/intel/fce/10.1.008/lib"

# Use the following for Matlab + gfortran 4.4.x
TESTC99COMPLEX=test_c99_complex
FC=    ief77
CC=    gcc
CXX=   g++
# Uncomment this line for new-style classdef support
# OOFLAG=-DR2008OO
###MEX=   mex $(OOFLAG)
# Use the following for 64-bit MEX
MEX=   mex -largeArrayDims $(OOFLAG)
MEX_FC=   $(FC)
MEX_CC=   $(CXX)
MW=mwrap
MWFLAGS=-c99complex -mex	
MEX_LIBS=-Wl,--rpath,$(intel_libs) -L$(intel_libs) -lifcoremt -lsvml -limf -liomp5 -lintlc -lm
MEX_LINK=-output $(MEX_GATEWAY).mexa64

MEX_FLAGS=
OBJSUF=o
MODSUF=mod

endif



ifeq ($(TARGET), matlab-maci)

FFLAGS=-g -O2 -fPIC

# Use the following for Matlab + gfortran 4.4.x
TESTC99COMPLEX=test_c99_complex
FC=    gfortran
CC=    gcc
CXX=   g++
# Uncomment this line for new-style classdef support
# OOFLAG=-DR2008OO
MEX=   mex $(OOFLAG) CC=gcc LD=gcc ARCH=maci
# Use the following for 64-bit MEX
###MEX=   mex -largeArrayDims $(OOFLAG) ARCH=maci
MEX_FC=   $(FC)
MEX_CC=   $(CXX)
MW=mwrap
MWFLAGS=-c99complex -mex	
MEX_LIBS=-lgfortran -lm
MEX_LINK=

MEX_FLAGS=
OBJSUF=o
MODSUF=mod

endif


ifeq ($(TARGET), matlab-maci-openmp)

FFLAGS=-g -O2 -fPIC --openmp

# Use the following for Matlab + gfortran 4.4.x
TESTC99COMPLEX=test_c99_complex
FC=    gfortran
CC=    gcc
CXX=   g++
# Uncomment this line for new-style classdef support
# OOFLAG=-DR2008OO
MEX=   mex $(OOFLAG) CC=gcc LD=gcc ARCH=maci
# Use the following for 64-bit MEX
###MEX=   mex -largeArrayDims $(OOFLAG) ARCH=maci
MEX_FC=   $(FC)
MEX_CC=   $(CXX)
MW=mwrap
MWFLAGS=-c99complex -mex	
MEX_LIBS=-lgfortran -lgomp -lm
MEX_LINK=

MEX_FLAGS=
OBJSUF=o
MODSUF=mod

endif


ifeq ($(TARGET), matlab-maci64)

FFLAGS=-g -O2 -fPIC

# Use the following for Matlab + gfortran 4.4.x
TESTC99COMPLEX=test_c99_complex
FC=    gfortran
CC=    gcc
CXX=   g++
# Uncomment this line for new-style classdef support
# OOFLAG=-DR2008OO
MEX=   mex $(OOFLAG) CC=gcc LD=gcc ARCH=maci64
# Use the following for 64-bit MEX
###MEX=   mex -largeArrayDims $(OOFLAG) ARCH=maci64
MEX_FC=   $(FC)
MEX_CC=   $(CXX)
MW=mwrap
MWFLAGS=-c99complex -mex        
MEX_LIBS=-lgfortran -lm
MEX_LINK=

MEX_FLAGS=
OBJSUF=o
MODSUF=mod

endif


ifeq ($(TARGET), matlab-maci64-openmp)

FFLAGS=-g -O2 -fPIC --openmp

# Use the following for Matlab + gfortran 4.4.x
TESTC99COMPLEX=test_c99_complex
FC=    gfortran
CC=    gcc
CXX=   g++
# Uncomment this line for new-style classdef support
# OOFLAG=-DR2008OO
MEX=   mex $(OOFLAG) CC=gcc LD=gcc ARCH=maci64
# Use the following for 64-bit MEX
###MEX=   mex -largeArrayDims $(OOFLAG) ARCH=maci64
MEX_FC=   $(FC)
MEX_CC=   $(CXX)
MW=mwrap
MWFLAGS=-c99complex -mex        
MEX_LIBS=-lgfortran -lgomp -lm
MEX_LINK=

MEX_FLAGS=
OBJSUF=o
MODSUF=mod

endif


ifeq ($(TARGET), matlab-linux-a64)

FFLAGS=-g -O2 -fPIC

# Use the following for Matlab + gfortran 4.4.x
TESTC99COMPLEX=test_c99_complex
FC=    gfortran
CC=    gcc
CXX=   g++
# Uncomment this line for new-style classdef support
# OOFLAG=-DR2008OO
###MEX=   mex $(OOFLAG)
# Use the following for 64-bit MEX
MEX=   mex -largeArrayDims $(OOFLAG)
MEX_FC=   $(FC)
MEX_CC=   $(CXX)
MW=mwrap
MWFLAGS=-c99complex -mex	
MEX_LIBS=-lgfortran -lm
MEX_LINK=-output $(MEX_GATEWAY)_single.mexa64

MEX_FLAGS=
OBJSUF=o
MODSUF=mod

endif



ifeq ($(TARGET), matlab-linux-a64-openmp)

FFLAGS=-g -O3 --openmp -fPIC

# Use the following for Matlab + gfortran 4.4.x
TESTC99COMPLEX=test_c99_complex
FC=    gfortran
CC=    gcc
CXX=   g++
# Uncomment this line for new-style classdef support
# OOFLAG=-DR2008OO
###MEX=   mex $(OOFLAG)
# Use the following for 64-bit MEX
MEX=   mex -largeArrayDims $(OOFLAG)
MEX_FC=   $(FC)
MEX_CC=   $(CXX)
MW=mwrap
MWFLAGS=-c99complex -mex	
MEX_LIBS=-lgomp -lgfortran -lm
MEX_LINK=-output $(MEX_GATEWAY).mexa64

MEX_FLAGS=
OBJSUF=o
MODSUF=mod

endif



ifeq ($(TARGET), matlab-windows-w32)

# Use the following for Matlab + gfortran cross compiler.
TESTC99COMPLEX=test_c99_complex

CROSS_FC=$(HOME)/lib/cross_win32/bin/i686-w64-mingw32-gfortran -g -O2 -ffast-math -msse2 -shared -fexceptions
CROSS_CC=$(HOME)/lib/cross_win32/bin/i686-w64-mingw32-gcc -g -O2 -ffast-math -msse2 -shared -fexceptions
CROSS_LIBS=-I$(HOME)/lib/matlab-r2006a/extern/include -L$(HOME)/lib/matlab-r2006a/bin/win32 -lmex -lmat -lmx -lm -lgfortran -static

FC=    $(CROSS_FC)
CC=    $(CROSS_CC)
CXX=   $(CROSS_CC)
MEX=   $(CXX)
MEX_FC=   $(FC)
MEX_CC=   $(CXX)
MW=mwrap
MWFLAGS=-c99complex -mex	
MEX_LIBS=$(CROSS_LIBS) -Wl,--no-undefined 
MEX_LINK=-o $(MEX_GATEWAY)_single.mexw32

MEX_FLAGS=
OBJSUF=o
MODSUF=mod

endif


ifeq ($(TARGET), matlab-windows-w32-openmp)

# Use the following for Matlab + gfortran cross compiler.
TESTC99COMPLEX=test_c99_complex

CROSS_FC=$(HOME)/lib/cross_win32/bin/i686-w64-mingw32-gfortran -g -O3 -shared -fexceptions
CROSS_CC=$(HOME)/lib/cross_win32/bin/i686-w64-mingw32-gcc -g -O3 -shared -fexceptions
CROSS_LIBS=-I$(HOME)/lib/matlab-r2006a/extern/include -L$(HOME)/lib/matlab-r2006a/bin/win32 -lmex -lmat -lmx -lm -lgfortran -static

FC=    $(CROSS_FC)  --openmp
CC=    $(CROSS_CC)  --openmp
CXX=   $(CROSS_CC)  --openmp
MEX=   $(CXX)
MEX_FC=   $(FC)
MEX_CC=   $(CXX)
MW=mwrap
MWFLAGS=-c99complex -mex
MEX_LIBS=$(CROSS_LIBS) --openmp -lgomp -lpthread -Wl,--no-undefined 
MEX_LINK=-o $(MEX_GATEWAY).mexw32

MEX_FLAGS=
OBJSUF=o
MODSUF=mod

endif


ifeq ($(TARGET), matlab-windows-w64)

# Use the following for Matlab + gfortran cross compiler.
TESTC99COMPLEX=test_c99_complex

CROSS_FC=$(HOME)/lib/cross_win64/bin/x86_64-w64-mingw32-gfortran -g -O2 -ffast-math -msse2 -shared -fexceptions
CROSS_CC=$(HOME)/lib/cross_win64/bin/x86_64-w64-mingw32-gcc -g -O2 -ffast-math -msse2 -shared -fexceptions
CROSS_LIBS=-I$(HOME)/lib/matlab-r2010b/extern/include -L$(HOME)/lib/matlab-r2010b/bin/win64 -lmex -lmat -lmx -lm -lgfortran -static

FC=    $(CROSS_FC)
CC=    $(CROSS_CC)
CXX=   $(CROSS_CC)
MEX=   $(CXX)
MEX_FC=   $(FC)
MEX_CC=   $(CXX)
MW=mwrap
MWFLAGS=-c99complex -mex	
MEX_LIBS=$(CROSS_LIBS) -Wl,--no-undefined 
MEX_LINK=-o $(MEX_GATEWAY)_single.mexw64

MEX_FLAGS=
OBJSUF=o
MODSUF=mod

endif


ifeq ($(TARGET), matlab-windows-w64-openmp)

# Use the following for Matlab + gfortran cross compiler.
TESTC99COMPLEX=test_c99_complex

CROSS_FC=$(HOME)/lib/cross_win64/bin/x86_64-w64-mingw32-gfortran -g -O3 -shared -fexceptions
CROSS_CC=$(HOME)/lib/cross_win64/bin/x86_64-w64-mingw32-gcc -g -O3 -shared -fexceptions
CROSS_LIBS=-I$(HOME)/lib/matlab-r2010b/extern/include -L$(HOME)/lib/matlab-r2010b/bin/win64 -lmex -lmat -lmx -lm -lgfortran -static

FC=    $(CROSS_FC)  --openmp
CC=    $(CROSS_CC)  --openmp
CXX=   $(CROSS_CC)  --openmp
MEX=   $(CXX)
MEX_FC=   $(FC)
MEX_CC=   $(CXX)
MW=mwrap
MWFLAGS=-c99complex -mex
MEX_LIBS=$(CROSS_LIBS) --openmp -lgomp -lpthread -Wl,--no-undefined 
MEX_LINK=-o $(MEX_GATEWAY).mexw64

MEX_FLAGS=
OBJSUF=o
MODSUF=mod

endif



