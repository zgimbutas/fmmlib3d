# === Compilers ===

# (You don't need flex and bison unless you're planning to modify the
#  grammar -- the default installation comes with the relevant files
#  already in place.)

# Uncomment this line for new-style classdef support
# OOFLAG=-DR2008OO

# Uncomment this line for C99 complex support
TESTC99COMPLEX=test_c99_complex

CC=    gcc
CXX=   g++
MEX=   mex $(OOFLAG)

# Use the following for 64-bit MEX
# MEX=   mex -largeArrayDims $(OOFLAG)

# Use the following for GNU Octave.
MEX=   mkoctfile --mex

FLEX=  flex
BISON= bison
