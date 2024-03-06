# Helmholtz and Laplace FMM library in R^3.

Date: March 6, 2024

Version 1.2.4

The FMMLIB3D suite permits the evaluation of potential fields governed
by the Laplace, Helmholtz equations in free space.

FMMLIB3D provides subroutines for particle (point) sources as well as
constant layer potential densities on triangles. The codes are easy to
use and reasonably well optimized for performance. A rudimentary
manual is provided in the FMM3D/doc directory.

FMMLIB3D contains both Fortran source code and versions compiled for
MATLAB under Mac OS X (64 bit), Windows (64 bit), and Linux (64 bit).


### License

```
Copyright (C) 2010-2012: Leslie Greengard and Zydrunas Gimbutas
Contact: greengard@cims.nyu.edu

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 

2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holders nor the names of its
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

### Contents

```
src/ - Fortran source code
examples/ - Fortran testing drivers and makefiles
matlab/ - matlab scripts and mex files 
contrib/mwrap-1.1.1/ - mwrap source code
```

To test the library, please type `make test`. 


### Fortran

Particle FMM routines.

```
hfmm3dpartself - Helmholtz particle FMM in R^3.
lfmm3dpartself - Laplace particle FMM in R^3.

hfmm3dparttarg - Helmholtz particle target FMM in R^3.
lfmm3dparttarg - Laplace particle target FMM in R^3.
```

Triangle FMM routines (constant densities on flat triangles).

```
hfmm3dtriaself - Helmholtz particle FMM in R^3.
lfmm3dtriaself - Laplace particle FMM in R^3.

hfmm3dtriatarg - Helmholtz particle target FMM in R^3.
lfmm3dtriatarg - Laplace particle target FMM in R^3.
```

Direct evaluation routines (constant densities on flat triangles).

```
h2dtriadirect - Helmholtz triangle interactions in R^3.
l2dtriadirect - Laplace triangle interactions in R^3.
```

Direct evaluation routines.

```
h2dpartdirect - Helmholtz particle interactions in R^3.
l2dpartdirect - Laplace particle interactions in R^3.
```


### Matlab

```
% Helmholtz and Laplace FMMs in R^3.
%
% Triangle FMM routines (constant densities on flat triangles).
%   hfmm3dtria      - Helmholtz triangle FMM in R^3. 
%   lfmm3dtria      - Laplace triangle FMM in R^3.
%
% Particle FMM routines.
%   hfmm3dpart      - Helmholtz particle FMM in R^3.
%   lfmm3dpart      - Laplace particle FMM in R^3.
%
% Direct evaluation routines (constant densities on flat triangles).
%   h3dtriadirect  - Helmholtz triangle interactions in R^3.
%   l3dtriadirect  - Laplace triangle interactions in R^3.
%
% Direct evaluation routines (particles).
%   h3dpartdirect  - Helmholtz particle interactions in R^3.
%   l3dpartdirect  - Laplace particle interactions in R^3.
%
% Triangulations.
%   atriread - Retrieve Cart3d triangulation from a file. (flat)
%   atriwrite - Store Cart3d triangulation to a file. (flat)
%   atriproc - Process triangulations in Cart3d format. (flat)
%   atrirefine - Refine Cart3d triangulation. (flat)
%   atriplot - Plot Cart3d triangulation. (flat)
%
% Triangle FMM postprocessing routines (constant densities on flat triangles).
%   hfmm3dtriampf    - Helmholtz triangle FMM in R^3, targets only.
%
% Tree generation routines.
%   d3tstrcr - construct the logical structure for a fully adaptive FMM in R^3.
%   d3tstrcrem  - include empty boxes, min and max level restriction.
%   d3tgetb     - retrieve box information.
%   d3tgetl     - retrieve list information.
%
% Testing and debugging routines.
%   test_lfmm3dpart_direct - test Laplace particle FMM and direct routines.
%   test_hfmm3dpart_direct - test Helmholtz particle FMM and direct routines.
%   test_lfmm3dtria_direct - test Laplace triangle FMM and direct routines.
%   test_hfmm3dtria_direct - test Helmholtz triangle FMM and direct routines.
%
% Internal utility functions.
%   fmm3dprini   - initialize simple printing routines.
%
```

### Acknowledgments

This work was supported in part by the Department of Energy and in
part by the Air Force Office of Scientific Research under MURI grant
FA9550-06-1-0337 and NSSEFF Program Award FA9550-10-0180, in part by
the NSF under grant DMS09-34733, and in part by Meyer Sound
Laboratories, Inc.

