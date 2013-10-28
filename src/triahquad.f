cc Copyright (C) 2009-2010: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
c
c
c    $Date: 2011-07-15 09:56:04 -0400 (Fri, 15 Jul 2011) $
c    $Revision: 2240 $
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  
c       Triangle quadrature routines
c
c       This subroutine returns values of integrals 
c       of constant density potentials on flat triangles
c
c       1/r and the derivatives of 1/r up to 3-th order
c       r and the derivatives of r up to 4-th order
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains a set of subroutines for the handling of 
c       of integrals of constant density potentials on flat triangles.
c       It contains 6 user-callable subroutines
c       
c    triahquad - returns values of constant density potentials on flat
c       triangles: used to derive the values of potentials of 1/r and
c       the derivatives of 1/r up to 3-th order, r and the derivatives r
c       up to 4-th order. See triahfun for the detailed description of
c       kernels, and triartable, triabtable for the mapping of different
c       types into the corresponding r and 1/r value and derivative tables.
c       The evaluation point is assumed to be in the near field.
c
c    triartable - returns values of constant density potentials on
c       flat triangles: 1/r and the derivatives of 1/r up to 3-th order
c
c    triabtable - returns values of constant density potentials on
c       flat triangles: r and the derivatives of r up to 4-th order
c
c    triahfun - returns values of kernels to be integrated, 
c       source in (xy) plane
c
c    triaevalp - returns values of constant density potentials on flat
c       triangles: used to derive the values of potentials of 1/r and
c       the derivatives of 1/r up to 3-th order, r and the derivatives r
c       up to 4-th order. See triaefun for the detailed description of
c       kernels, and triartable, triabtable for the mapping of different
c       types into the corresponding r and 1/r value and derivative
c       tables.  The evaluation point is assumed to be in the far field,
c       and all integrals are evaluated via a quadrature formula with
c       user provided nodes. This is much faster than triahquad for
c       small number of quadratures nodes assuming that that the target
c       point is separated by at least 3 or 4 source triangle diameters
c       from the source centroid. The source triangle is assumed to be
c       in R^3.
c
c    triaefun - returns values of kernels to be integrated, source in R^3
c
c
        subroutine triahquad
     $     (itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        implicit none
        integer  itype,iquad
        real *8 x0,y0,z0,hval
        real *8 vert1(2),vert2(2),vert3(2)
        real *8 vert0(2),w(12),v1(2),v2(2),v3(2),verth(2)
        real *8 nx,ny
        real *8 a1,a2,cosa,sina,ds,dx,dy,h,rval
c
c       This subroutine returns values of constant density potentials on
c       flat triangles: used to derive the values of potentials of 1/r
c       and the derivatives of 1/r up to 3-th order, r and the
c       derivatives r up to 4-th order. See triahfun for the detailed
c       description of kernels, and triartable, triabtable for the
c       mapping of different types into the corresponding r and 1/r
c       value and derivative tables.
c
c       The evaluation point is assumed to be in the near field.
c       
c       All integrals are evaluated at the point (x0,y0,z0)
c
c       
c       Input parameters: 
c
c       itype - type of the integral to be evaluated, see below 
c              must be in range (1..51)
c
c       iquad - this should be the sign of z0
c              iquad =+1, if z0>0
c              iquad = 0, if z0=0
c              iquad =-1, if z0<0
c
c       vert1(2), vert2(2), vert3(2) - vertices of triangle
c
c       x0,y0,y0 - target evaluation point
c
c       Output parameters:
c
c       hval - integral value (real*8)
c       
c
c
c       The triangle T is defined in R^3 and its vertices are
c
c       vert1(2), vert2(2), vert3(2) in R^2
c
c       
c            /\vert3
c           /  \
c          /    \
c         /   T  \
c        /________\
c       vert1     vert2
c
c       L = boundary of triangle T
c
c       hval is the defined as 
c
c       \int_T  f_{itype} (x,y,z) dT
c
c       itype=1:   \int_T  1/r dT
c       itype=2:   \int_T  (x-x0)/r^3 dT
c       itype=3:   \int_T  (y-y0)/r^3 dT
c       itype=4:   \int_T  z0/r^3 dT
c       itype=5:   \int_T  3*(x-x0)*z0/r^5 dT
c       itype=6:   \int_T  3*(y-y0)*z0/r^5 dT
c       itype=7:   \int_T  -1/r^3+3*z0*z0/r^5 dT
c
c       itype=8:   \int_T  (x-x0)/r dT
c       itype=9:   \int_T  (y-y0)/r dT
c
c       itype=10:   \int_T  (x-x0)*z0/r^3 dT
c       itype=11:   \int_T  (y-y0)*z0/r^3 dT
c
c       itype=12:   \int_T  -1/r^3+3*(x-x0)*(x-x0)/r^5 dT
c       itype=13:   \int_T  3*(x-x0)*(y-y0)/r^5 dT
c       itype=14:   \int_T  -1/r^3+3*(y-y0)*(y-y0)/r^5 dT
c
c       itype=15:   \int_T  3*z0/r^5-15*(x-x0)*(x-x0)*z0/r^7 dT
c       itype=16:   \int_T  -15*(x-x0)*(y-y0)*z0/r^7 dT
c       itype=17:   \int_T  3*z0/r^5-15*(y-y0)*(y-y0)*z0/r^7 dT
c
c       itype=18:   \int_T  3*(x-x0)/r^5-15*(x-x0)*z0*z0/r^7 dT
c       itype=19:   \int_T  3*(y-y0)/r^5-15*(y-y0)*z0*z0/r^7 dT
c       itype=20:   \int_T  9*z0/r^5-15*z0*z0*z0/r^7 dT
c
c       itype=21:   \int_T  9*(x-x0)/r^5-15*(x-x0)*(x-x0)*(x-x0)/r^7 dT
c       itype=22:   \int_T  3*(y-y0)/r^5-15*(x-x0)*(x-x0)*(y-y0)/r^7 dT
c       itype=23:   \int_T  3*(x-x0)/r^5-15*(x-x0)*(y-y0)*(y-y0)/r^7 dT
c       itype=24:   \int_T  9*(y-y0)/r^5-15*(y-y0)*(y-y0)*(y-y0)/r^7 dT
c
c       itype=25:   \int_T  -1/r+(x-x0)*(x-x0)/r^3 dT
c       itype=26:   \int_T  (x-x0)*(y-y0)/r^3 dT
c       itype=27:   \int_T  -1/r+(y-y0)*(y-y0)/r^3 dT
c
c       itype=28:   \int_T  -z0/r^3+3*(x-x0)*(x-x0)*z0/r^5 dT
c       itype=29:   \int_T          3*(x-x0)*(y-y0)*z0/r^5 dT
c       itype=30:   \int_T  -z0/r^3+3*(y-y0)*(y-y0)*z0/r^5 dT
c
c       itype=31:   \int_T  (x-x0)*(-3/r^3+3*(x-x0)*(x-x0)/r^5) dT
c       itype=32:   \int_T  (y-y0)*(-1/r^3+3*(x-x0)*(x-x0)/r^5) dT
c       itype=33:   \int_T  (x-x0)*(-1/r^3+3*(y-y0)*(y-y0)/r^5) dT
c       itype=34:   \int_T  (y-y0)*(-3/r^3+3*(y-y0)*(y-y0)/r^5) dT
c
c       itype=35:   \int_T  (x-x0)*(-1/r^3+3*z0*z0/r^5) dT
c       itype=36:   \int_T  (y-y0)*(-1/r^3+3*z0*z0/r^5) dT
c
c       itype=37:   \int_T  r dT  
c
c       itype=38:   -3/r^3+18*(x-x0)^2/r^5-15*(x-x0)^4/r^7
c       itype=39:   9*(x-x0)*(y-y0)/r^5-15*(x-x0)^3*(y-y0)/r^7
c       itype=40:   -1/r^3+3*((x-x0)^2+(y-y0)^2)/r^5-15*(x-x0)^2*(y-y0)^2/r^7
c       itype=41:   9*(x-x0)*(y-y0)/r^5-15*(x-x0)*(y-y0)^3/r^7
c       itype=42:   -3/r^3+18*(y-y0)^2/r^5-15*(y-y0)^4/r^7
c
c       itype=43:   9*(x-x0)*z0/r^5-15*(x-x0)^3*z0/r^7
c       itype=44:   3*(y-y0)*z0/r^5-15*(x-x0)^2*(y-y0)*z0/r^7
c       itype=45:   3*(x-x0)*z0/r^5-15*(x-x0)*(y-y0)^2*z0/r^7
c       itype=46:   9*(y-y0)*z0/r^5-15*(y-y0)^3*z0/r^7
c
c       itype=47:   -1/r^3+3*((x-x0)^2+z0^2)/r^5-15*(x-x0)^2*z0^2/r^7
c       itype=48:   3*(x-x0)*(y-y0)/r^5-15*(x-x0)*(y-y0)*z0^2/r^7
c       itype=49:   -1/r^3+3*((y-y0)^2+z0^2)/r^5-15*(y-y0)^2*z0^2/r^7
c
c       itype=50:   9*(x-x0)*z0/r^5-15*(x-x0)*z0^3/r^7
c       itype=51:   9*(y-y0)*z0/r^5-15*(y-y0)*z0^3/r^7
c
c       where r=1/sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2), with z=0
c
c       iquad is the sign of z0
c
c       iquad =+1, if z0>0
c       iquad = 0, if z0=0
c       iquad =-1, if z0<0
c
c
c         diff(r,z)       -> z0/r
c         diff(r,z,z)     -> 1/r-z0^2/r^3
c         diff(r,z,z,z)   -> -3*z0/r^3+3*z0^3/r^5
c         diff(r,z,z,z,z) -> -3/r^3+18*z0^2/r^5-15*z0^4/r^7
c
c
c  2      diff(1/r,x)=-(x-x0)/r^3   
c  3      diff(1/r,y)=-(y-y0)/r^3
c  4      diff(1/r,z)=z0/r^3
c
c  5      diff(1/r,z,x)=3*(x-x0)*z0/r^5
c  6      diff(1/r,z,y)=3*(y-y0)*z0/r^5
c  7      diff(1/r,z,z)=-1/r^3+3*z0*z0/r^5
c
c 12      diff(1/r,x,x)=-1/r^3+3*(x-x0)*(x-x0)/r^5
c 13      diff(1/r,x,y)=       3*(x-x0)*(y-y0)/r^5
c 14      diff(1/r,y,y)=-1/r^3+3*(y-y0)*(y-y0)/r^5
c
c 15      diff(1/r,x,x,z)=3*z0/r^5-15*(x-x0)*(x-x0)*z0/r^7
c 16      diff(1/r,x,y,z)=        -15*(x-x0)*(y-y0)*z0/r^7
c 17      diff(1/r,y,y,z)=3*z0/r^5-15*(y-y0)*(y-y0)*z0/r^7
c
c 18      diff(1/r,x,z,z)=3*(x-x0)/r^5-15*(x-x0)*z0*z0/r^7
c 19      diff(1/r,y,z,z)=3*(y-y0)/r^5-15*(y-y0)*z0*z0/r^7
c 20      diff(1/r,z,z,z)=9*z0/r^5    -15*z0*z0*z0/r^7
c
c 21      diff(1/r,x,x,x)=9*(x-x0)/r^5-15*(x-x0)*(x-x0)*(x-x0)/r^7
c 22      diff(1/r,x,x,y)=3*(y-y0)/r^5-15*(x-x0)*(y-y0)*(x-x0)/r^7
c 23      diff(1/r,x,y,y)=3*(x-x0)/r^5-15*(y-y0)*(x-x0)*(x-x0)/r^7
c 24      diff(1/r,y,y,y)=9*(y-y0)/r^5-15*(y-y0)*(y-y0)*(y-y0)/r^7
c
c
c
c       Regular integrals
c              
c       (itype .eq. 1 .or. itype .eq. 4 .or. itype .eq. 7 .or. itype .eq. 20)
c
c       are computed by splitting triangle into 3 oriented triangles and
c       calling triarquad_ab routine
c              
c
c       Hilbert like integrals are handled via integration by parts formula,
c       aka the gradient identity \int_T \grad f dT = \int_L f n dL, where
c       n is the exterior normal and L is the boundary of T
c
c       \int_T (x-x0)/sqrt((x-x0)^2+(y-y0)^2+z0^3)^3 dT = 
c       -\int_T d/dx 1/sqrt((x-x0)^2+(y-y0)^2+z0^3) dT =
c       -\int_L 1/sqrt((x-x0)^2+(y-y0)^2+z0^2) n_x dL
c
c       \int_T 3*(x-x0)*z0/sqrt((x-x0)^2+(y-y0)^2+z0^2)^5 dT = 
c       -\int_T d/dx z0/sqrt((x-x0)^2+(y-y0)^2+z0^2)^3 dT =
c       -\int_L z0/sqrt((x-x0)^2+(y-y0)^2+z0^2)^3 n_x dL
c
c       They are evaluated by splitting the boundary of triangle into 3
c       oriented segments and calling trialquad and trialquad_hd routines
c
c       Warning: internally, we compute the INWARD normal to the triangle,
c       an all integrands are adjusted to get the signs right
c
c       Regular integral for r
c              
c       (itype .eq. 37 )
c
c       is computed by splitting triangle into 6 oriented triangles and
c       calling triarquad routine
c              
c              
c
c
c
        hval=0
c
        if( itype .eq. 1 .or. itype .eq. 4 .or.
     $     itype .eq. 7 .or. itype .eq. 20 ) then
c
c       ... 3 triangle subdivision algorithm
c
        vert0(1)=x0
        vert0(2)=y0
c
        hval=0
c
        call triahproj(vert0,vert1,vert2,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c       
        call triarquad_ab(itype,iquad,-a1,a2,h,z0,rval)
        hval=hval+rval
c
        call triahproj(vert0,vert2,vert3,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c
        call triarquad_ab(itype,iquad,-a1,a2,h,z0,rval)
        hval=hval+rval 
c
        call triahproj(vert0,vert3,vert1,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c
        call triarquad_ab(itype,iquad,-a1,a2,h,z0,rval)
        hval=hval+rval 
c
        return
        endif
c
c
c       ... Hilbert integrals
c
        if( itype .eq. 2 .or. itype .eq. 5 .or. itype .eq. 8
     $     .or. itype .eq. 10 .or. itype .eq. 18 
     $     .or. itype .eq. 35 .or. itype .eq. 50 ) then
c
        vert0(1)=x0
        vert0(2)=y0
c
        hval=0
c
        call triahproj(vert0,vert1,vert2,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c       
        call trialquad(itype,iquad,-a1,a2,h,z0,rval)
        dx=vert2(1)-vert1(1)
        dy=vert2(2)-vert1(2)
        ds=sqrt(dx**2+dy**2)        
        nx=-dy/ds
        ny= dx/ds
        rval=rval*nx
        hval=hval+rval
c
        call triahproj(vert0,vert2,vert3,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c
        call trialquad(itype,iquad,-a1,a2,h,z0,rval)
        dx=vert3(1)-vert2(1)
        dy=vert3(2)-vert2(2)
        ds=sqrt(dx**2+dy**2)
        nx=-dy/ds
        ny= dx/ds
        rval=rval*nx
        hval=hval+rval 
c
        call triahproj(vert0,vert3,vert1,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c
        call trialquad(itype,iquad,-a1,a2,h,z0,rval)
        dx=vert1(1)-vert3(1)
        dy=vert1(2)-vert3(2)
        ds=sqrt(dx**2+dy**2)        
        nx=-dy/ds
        ny= dx/ds
        rval=rval*nx
        hval=hval+rval 
c
        return
        endif
c
c
        if( itype .eq. 3 .or. itype .eq. 6 .or. itype .eq. 9
     $     .or. itype .eq. 11 .or. itype .eq. 19 
     $     .or. itype .eq. 36 .or. itype .eq. 51 ) then
c
        vert0(1)=x0
        vert0(2)=y0
c
        hval=0
c
        call triahproj(vert0,vert1,vert2,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c       
        call trialquad(itype,iquad,-a1,a2,h,z0,rval)
        dx=vert2(1)-vert1(1)
        dy=vert2(2)-vert1(2)
        ds=sqrt(dx**2+dy**2)        
        nx=-dy/ds
        ny= dx/ds
        rval=rval*ny
        hval=hval+rval
c
        call triahproj(vert0,vert2,vert3,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c
        call trialquad(itype,iquad,-a1,a2,h,z0,rval)
        dx=vert3(1)-vert2(1)
        dy=vert3(2)-vert2(2)
        ds=sqrt(dx**2+dy**2)        
        nx=-dy/ds
        ny= dx/ds
        rval=rval*ny
        hval=hval+rval
c
        call triahproj(vert0,vert3,vert1,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c
        call trialquad(itype,iquad,-a1,a2,h,z0,rval)
        dx=vert1(1)-vert3(1)
        dy=vert1(2)-vert3(2)
        ds=sqrt(dx**2+dy**2)        
        nx=-dy/ds
        ny= dx/ds
        rval=rval*ny
        hval=hval+rval
c
        return
        endif
c
c
c       ... tangential derivatives of Hilbert integrals
c
        if( itype .eq. 12 .or. itype .eq. 15 
     $     .or. itype .eq. 21 .or. itype .eq. 25 
     $     .or. itype .eq. 28 .or. itype .eq. 31 .or. itype .eq. 38
     $     .or. itype .eq. 40
     $     .or. itype .eq. 43 .or. itype .eq. 47 ) then
c
        vert0(1)=x0
        vert0(2)=y0
c
        hval=0
c
        call triahproj(vert0,vert1,vert2,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c       
        dx=vert2(1)-vert1(1)
        dy=vert2(2)-vert1(2)
        ds=sqrt(dx**2+dy**2)        
        nx=-dy/ds
        ny= dx/ds
        cosa=-ny
        sina=+nx
        call trialquad_hd(itype,iquad,cosa,sina,-a1,a2,h,z0,rval)
        rval=rval*nx
        hval=hval+rval
c
        call triahproj(vert0,vert2,vert3,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c
        dx=vert3(1)-vert2(1)
        dy=vert3(2)-vert2(2)
        ds=sqrt(dx**2+dy**2)
        nx=-dy/ds
        ny= dx/ds
        cosa=-ny
        sina=+nx
        call trialquad_hd(itype,iquad,cosa,sina,-a1,a2,h,z0,rval)
        rval=rval*nx
        hval=hval+rval 
c
        call triahproj(vert0,vert3,vert1,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c
        dx=vert1(1)-vert3(1)
        dy=vert1(2)-vert3(2)
        ds=sqrt(dx**2+dy**2)        
        nx=-dy/ds
        ny= dx/ds
        cosa=-ny
        sina=+nx
        call trialquad_hd(itype,iquad,cosa,sina,-a1,a2,h,z0,rval)
        rval=rval*nx
        hval=hval+rval 
c
        return
        endif
c
c
        if( itype .eq. 13 .or. itype .eq. 16 
     $     .or. itype .eq. 22 .or. itype .eq. 26 
     $     .or. itype .eq. 29 .or. itype .eq. 32 .or. itype .eq. 39
     $     .or. itype .eq. 44 .or. itype .eq. 48 ) then
c
        vert0(1)=x0
        vert0(2)=y0
c
        hval=0
c
        call triahproj(vert0,vert1,vert2,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c       
        dx=vert2(1)-vert1(1)
        dy=vert2(2)-vert1(2)
        ds=sqrt(dx**2+dy**2)        
        nx=-dy/ds
        ny= dx/ds
        cosa=-ny
        sina=+nx
        call trialquad_hd(itype,iquad,cosa,sina,-a1,a2,h,z0,rval)
        rval=rval*ny
        hval=hval+rval
c
        call triahproj(vert0,vert2,vert3,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c
        dx=vert3(1)-vert2(1)
        dy=vert3(2)-vert2(2)
        ds=sqrt(dx**2+dy**2)
        nx=-dy/ds
        ny= dx/ds
        cosa=-ny
        sina=+nx
        call trialquad_hd(itype,iquad,cosa,sina,-a1,a2,h,z0,rval)
        rval=rval*ny
        hval=hval+rval 
c
        call triahproj(vert0,vert3,vert1,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c
        dx=vert1(1)-vert3(1)
        dy=vert1(2)-vert3(2)
        ds=sqrt(dx**2+dy**2)        
        nx=-dy/ds
        ny= dx/ds
        cosa=-ny
        sina=+nx
        call trialquad_hd(itype,iquad,cosa,sina,-a1,a2,h,z0,rval)
        rval=rval*ny
        hval=hval+rval 
c
        return
        endif
c
c
        if( itype .eq. 23 .or. itype .eq. 33 .or. itype .eq. 41 
     $     .or. itype .eq. 45 ) then
c
        vert0(1)=x0
        vert0(2)=y0
c
        hval=0
c
        call triahproj(vert0,vert1,vert2,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c       
        dx=vert2(1)-vert1(1)
        dy=vert2(2)-vert1(2)
        ds=sqrt(dx**2+dy**2)        
        nx=-dy/ds
        ny= dx/ds
        cosa=+nx
        sina=+ny
        call trialquad_hd(itype,iquad,cosa,sina,-a1,a2,h,z0,rval)
        rval=rval*nx
        hval=hval+rval
c
        call triahproj(vert0,vert2,vert3,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c
        dx=vert3(1)-vert2(1)
        dy=vert3(2)-vert2(2)
        ds=sqrt(dx**2+dy**2)
        nx=-dy/ds
        ny= dx/ds
        cosa=+nx
        sina=+ny
        call trialquad_hd(itype,iquad,cosa,sina,-a1,a2,h,z0,rval)
        rval=rval*nx
        hval=hval+rval 
c
        call triahproj(vert0,vert3,vert1,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c
        dx=vert1(1)-vert3(1)
        dy=vert1(2)-vert3(2)
        ds=sqrt(dx**2+dy**2)        
        nx=-dy/ds
        ny= dx/ds
        cosa=+nx
        sina=+ny
        call trialquad_hd(itype,iquad,cosa,sina,-a1,a2,h,z0,rval)
        rval=rval*nx
        hval=hval+rval 
c
        endif
c
        if( itype .eq. 14 .or. itype .eq. 17 
     $     .or. itype .eq. 24 .or. itype .eq. 27 
     $     .or. itype .eq. 30 .or. itype .eq. 34 .or. itype .eq. 42
     $     .or. itype .eq. 46 .or. itype .eq. 49 ) then
c
        vert0(1)=x0
        vert0(2)=y0
c
        hval=0
c
        call triahproj(vert0,vert1,vert2,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c       
        dx=vert2(1)-vert1(1)
        dy=vert2(2)-vert1(2)
        ds=sqrt(dx**2+dy**2)        
        nx=-dy/ds
        ny= dx/ds
        cosa=+nx
        sina=+ny
        call trialquad_hd(itype,iquad,cosa,sina,-a1,a2,h,z0,rval)
        rval=rval*ny
        hval=hval+rval
c
        call triahproj(vert0,vert2,vert3,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c
        dx=vert3(1)-vert2(1)
        dy=vert3(2)-vert2(2)
        ds=sqrt(dx**2+dy**2)
        nx=-dy/ds
        ny= dx/ds
        cosa=+nx
        sina=+ny
        call trialquad_hd(itype,iquad,cosa,sina,-a1,a2,h,z0,rval)
        rval=rval*ny
        hval=hval+rval 
c
        call triahproj(vert0,vert3,vert1,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c
        dx=vert1(1)-vert3(1)
        dy=vert1(2)-vert3(2)
        ds=sqrt(dx**2+dy**2)        
        nx=-dy/ds
        ny= dx/ds
        cosa=+nx
        sina=+ny
        call trialquad_hd(itype,iquad,cosa,sina,-a1,a2,h,z0,rval)
        rval=rval*ny
        hval=hval+rval 
c
        return
        endif
c
c
        if( itype .eq. 37 ) then
c       
c       ... 6 triangle subdivision algorithm 
c
        vert0(1)=x0
        vert0(2)=y0
c
        hval=0
c
        call triahproj(vert0,vert1,vert2,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c
        call triarquad(itype,iquad,h,a1,z0,rval)
        hval=hval+rval
        call triarquad(itype,iquad,h,a2,z0,rval)
        hval=hval+rval
c
        call triahproj(vert0,vert2,vert3,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c
        call triarquad(itype,iquad,h,a1,z0,rval)
        hval=hval+rval
        call triarquad(itype,iquad,h,a2,z0,rval)
        hval=hval+rval
c
        call triahproj(vert0,vert3,vert1,h,a1,a2,verth)
ccc        write(*,*) a1,a2,h
c
        call triarquad(itype,iquad,h,a1,z0,rval)
        hval=hval+rval
        call triarquad(itype,iquad,h,a2,z0,rval)
        hval=hval+rval
c
        return
        endif
c
c
ccc        write(20,*) 'itype=', itype, 'hval=', hval
c
        return
        end
c
c
c
c
c
        subroutine triartable
     $     (ix,iy,iz,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        implicit none
        integer ix,iy,iz,iquad,itype
        real *8 vert1(2),vert2(2),vert3(2)
        real *8 vert0(2),w(12),v1(2),v2(2),v3(2),verth(2)
        real *8 x0,y0,z0,hval
c
c       This subroutine returns values of constant density potentials on
c       flat triangles: 1/r and the derivatives of 1/r up to 3-th order
c
c       Input parameters: 
c
c       ix,iy,iz - the order of derivative with respect of x,y,z, respectively
c
c       iquad - this should be the sign of z0
c              iquad =+1, if z0>0
c              iquad = 0, if z0=0
c              iquad =-1, if z0<0
c
c       vert1(2), vert2(2), vert3(2) - vertices of triangle
c
c       x0,y0,y0 - target evaluation point
c
c       Output parameters:
c
c       hval - integral value (real*8)
c       
c
c
        if( ix+iy+iz .eq. 0 ) then
        itype = 1
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        return
        endif
c
        if( ix+iy+iz .eq. 1 ) then
        if( iz .eq. 0 .and. ix .eq. 1 ) itype = 2
        if( iz .eq. 0 .and. iy .eq. 1 ) itype = 3
        if( iz .eq. 1 ) itype = 4
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        if( iz .eq. 0 .and. ix .eq. 1 ) hval=-hval
        if( iz .eq. 0 .and. iy .eq. 1 ) hval=-hval
        return
        endif
c
c
        if( ix+iy+iz .eq. 2 ) then
        if( iz .eq. 0 .and. ix .eq. 2 ) itype = 12
        if( iz .eq. 0 .and. ix .eq. 1 .and. iy .eq. 1 ) itype = 13
        if( iz .eq. 0 .and. iy .eq. 2 ) itype = 14
        if( iz .eq. 1 .and. ix .eq. 1 ) itype = 5
        if( iz .eq. 1 .and. iy .eq. 1 ) itype = 6
        if( iz .eq. 2 ) itype = 7
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        if( ix .eq. 1 .and. iz .eq. 1 ) hval=-hval
        if( iy .eq. 1 .and. iz .eq. 1 ) hval=-hval
        return
        endif
c
c
        if( ix+iy+iz .eq. 3 ) then
        if( iz .eq. 0 .and. ix .eq. 3 ) then
        itype = 21
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        endif
        if( iz .eq. 0 .and. ix .eq. 2 .and. iy .eq. 1 ) then
        itype = 22
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        endif
        if( iz .eq. 0 .and. ix .eq. 1 .and. iy .eq. 2 ) then
        itype = 23
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        endif
        if( iz .eq. 0 .and. iy .eq. 3 ) then
        itype = 24
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        endif
        if( iz .eq. 1 .and. ix .eq. 2 ) then
        itype = 15
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        hval=-hval
        endif
        if( iz .eq. 1 .and. ix .eq. 1 .and. iy .eq. 1 ) then
        itype = 16
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        hval=-hval
        endif
        if( iz .eq. 1 .and. iy .eq. 2 ) then
        itype = 17
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        hval=-hval
        endif
        if( iz .eq. 2 .and. ix .eq. 1 ) then
        itype = 18
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        endif
        if( iz .eq. 2 .and. iy .eq. 1 ) then
        itype = 19
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        endif
        if( iz .eq. 3 ) then
        itype = 20
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        hval=-hval
        endif
        return
        endif
c
c
        return
        end
c
c
c
c
c
        subroutine triabtable
     $     (ix,iy,iz,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        implicit none
        integer ix,iy,iz,iquad,itype
        real *8 vert1(2),vert2(2),vert3(2)
        real *8 x0,y0,z0,hval,rval
        real *8 vert0(2),w(12),v1(2),v2(2),v3(2),verth(2)
c
c       This subroutine returns values of constant density potentials on
c       flat triangles: r and the derivatives of r up to 4-th order
c
c       Input parameters: 
c
c       ix,iy,iz - the order of derivative with respect of x,y,z, respectively
c
c       iquad - this should be the sign of z0
c              iquad =+1, if z0>0
c              iquad = 0, if z0=0
c              iquad =-1, if z0<0
c
c       vert1(2), vert2(2), vert3(2) - vertices of triangle
c
c       x0,y0,y0 - target evaluation point
c
c       Output parameters:
c
c       hval - integral value (real*8)
c       
c
        if( ix+iy+iz .eq. 0 ) then
c
        itype=37
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
c
        return
        endif
c
c
        if( ix+iy+iz .eq. 1 ) then
        if( iz .eq. 0 .and. ix .eq. 1 ) itype = 8
        if( iz .eq. 0 .and. iy .eq. 1 ) itype = 9
        if( iz .eq. 1 ) itype = 1
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        if( iz .eq. 1 ) hval=-z0*hval
        return
        endif
c
        if( ix+iy+iz .eq. 2 ) then
        if( iz .eq. 0 .and. ix .eq. 2 ) then
        itype = 25
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        hval=-hval
        endif
        if( iz .eq. 0 .and. ix .eq. 1 .and. iy .eq. 1 ) then
        itype = 26
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        hval=-hval
        endif
        if( iz .eq. 0 .and. iy .eq. 2 ) then
        itype = 27
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        hval=-hval
        endif
        if( iz .eq. 1 .and. ix .eq. 1 ) then
        itype = 2
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        hval=hval*z0
        endif
        if( iz .eq. 1 .and. iy .eq. 1 ) then
        itype = 3
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        hval=hval*z0
        endif
        if( iz .eq. 2 ) then
        itype = 4
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        hval=-hval*z0
        itype = 1
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,rval)
        hval=hval+rval
        endif
        return
        endif
c
c
        if( ix+iy+iz .eq. 3 ) then
        if( iz .eq. 0 .and. ix .eq. 3 ) then
        itype = 31
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        endif
        if( iz .eq. 0 .and. ix .eq. 2 .and. iy .eq. 1 ) then
        itype = 32
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        endif
        if( iz .eq. 0 .and. ix .eq. 1 .and. iy .eq. 2 ) then
        itype = 33
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        endif
        if( iz .eq. 0 .and. iy .eq. 3 ) then
        itype = 34
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        endif
        if( iz .eq. 1 .and. ix .eq. 2 ) then
        itype = 28
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        hval=-hval
        endif
        if( iz .eq. 1 .and. ix .eq. 1 .and. iy .eq. 1 ) then
        itype = 29
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        hval=-hval
        endif
        if( iz .eq. 1 .and. iy .eq. 2 ) then
        itype = 30
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        hval=-hval
        endif
        if( iz .eq. 2 .and. ix .eq. 1 ) then
        itype = 35
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        endif
        if( iz .eq. 2 .and. iy .eq. 1 ) then
        itype = 36
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        endif
        if( iz .eq. 3 ) then
        itype = 7
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        hval=-hval*z0
        itype = 4
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,rval)
        hval=hval+2*rval
        endif
        return
        endif
c
c
        if( ix+iy+iz .eq. 4 ) then
        if( iz .eq. 0 .and. ix .eq. 4 ) then
        itype = 38
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        endif
        if( iz .eq. 0 .and. ix .eq. 3 .and. iy .eq. 1 ) then
        itype = 39
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        endif
        if( iz .eq. 0 .and. ix .eq. 2 .and. iy .eq. 2 ) then
        itype = 40
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        endif
        if( iz .eq. 0 .and. ix .eq. 1 .and. iy .eq. 3 ) then
        itype = 41
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        endif
        if( iz .eq. 0 .and. iy .eq. 4 ) then
        itype = 42
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        endif
        if( iz .eq. 1 .and. ix .eq. 3 ) then
        itype = 43
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        hval=-hval
        endif
        if( iz .eq. 1 .and. ix .eq. 2 .and. iy .eq. 1 ) then
        itype = 44
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        hval=-hval
        endif
        if( iz .eq. 1 .and. ix .eq. 1 .and. iy .eq. 2 ) then
        itype = 45
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        hval=-hval
        endif
        if( iz .eq. 1 .and. iy .eq. 3 ) then
        itype = 46
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        hval=-hval
        endif
        if( iz .eq. 2 .and. ix .eq. 2 ) then
        itype = 47
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        endif
        if( iz .eq. 2 .and. ix .eq. 1 .and. iy .eq. 1 ) then
        itype = 48
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        endif
        if( iz .eq. 2 .and. iy .eq. 2 ) then
        itype = 49
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        endif
        if( iz .eq. 3 .and. ix .eq. 1 ) then
        itype = 50
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        hval=-hval
        endif
        if( iz .eq. 3 .and. iy .eq. 1 ) then
        itype = 51
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        hval=-hval
        endif
        if( iz .eq. 4 ) then
        itype = 20
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,hval)
        hval=hval*z0
        itype = 7
        call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,rval)
        hval=hval+3*rval
        endif
        return
        endif
c
c
        return
        end
c
c
c
c
c
        subroutine triahfun(itype,x,y,x0,y0,z0,val)
        implicit none
        integer itype
        real *8 x,y,x0,y0,z0,val
        real *8 r
c
c       This subroutine returns the values of kernels to be integrated.
c       Integration domain is assumed to be in (x,y) plane with z=0
c
c       Input parameters: 
c
c       itype - type of the integral to be evaluated, see below 
c              must be in range (1..51)
c
c       x,y - source evaluation point, in (x,y) plane
c
c       x0,y0,y0 - target evaluation point
c
c       Output parameters:
c
c       val - kernel value (real*8)
c
c
c       
        val=0
c
        r=sqrt((x-x0)**2+(y-y0)**2+z0**2)
c
        if( itype .eq. 1 ) val=1/r
        if( itype .eq. 2 ) val=(x-x0)/r**3
        if( itype .eq. 3 ) val=(y-y0)/r**3
        if( itype .eq. 4 ) val=z0/r**3
        if( itype .eq. 5 ) val=3*(x-x0)*z0/r**5
        if( itype .eq. 6 ) val=3*(y-y0)*z0/r**5
        if( itype .eq. 7 ) val=-1/r**3+3*z0*z0/r**5
c
        if( itype .eq. 8 ) val=(x-x0)/r
        if( itype .eq. 9 ) val=(y-y0)/r
c
        if( itype .eq. 10 ) val=(x-x0)*z0/r**3
        if( itype .eq. 11 ) val=(y-y0)*z0/r**3
c
        if( itype .eq. 12 ) val=-1/r**3+3*(x-x0)*(x-x0)/r**5
        if( itype .eq. 13 ) val=        3*(x-x0)*(y-y0)/r**5
        if( itype .eq. 14 ) val=-1/r**3+3*(y-y0)*(y-y0)/r**5
c
        if( itype .eq. 15 ) val=3*z0/r**5-15*(x-x0)*(x-x0)*z0/r**7
        if( itype .eq. 16 ) val=         -15*(x-x0)*(y-y0)*z0/r**7
        if( itype .eq. 17 ) val=3*z0/r**5-15*(y-y0)*(y-y0)*z0/r**7
c
        if( itype .eq. 18 ) val=3*(x-x0)/r**5-15*(x-x0)*z0*z0/r**7
        if( itype .eq. 19 ) val=3*(y-y0)/r**5-15*(y-y0)*z0*z0/r**7
        if( itype .eq. 20 ) val=9*z0/r**5    -15*z0*z0*z0/r**7
c
        if( itype .eq. 21 ) 
     $     val=9*(x-x0)/r**5-15*(x-x0)*(x-x0)*(x-x0)/r**7
        if( itype .eq. 22 ) 
     $     val=3*(y-y0)/r**5-15*(x-x0)*(x-x0)*(y-y0)/r**7
        if( itype .eq. 23 ) 
     $     val=3*(x-x0)/r**5-15*(x-x0)*(y-y0)*(y-y0)/r**7
        if( itype .eq. 24 ) 
     $     val=9*(y-y0)/r**5-15*(y-y0)*(y-y0)*(y-y0)/r**7
c
        if( itype .eq. 25 ) val=-1/r+(x-x0)*(x-x0)/r**3
        if( itype .eq. 26 ) val=(x-x0)*(y-y0)/r**3
        if( itype .eq. 27 ) val=-1/r+(y-y0)*(y-y0)/r**3
c
        if( itype .eq. 28 ) val=-z0/r**3+3*(x-x0)*(x-x0)*z0/r**5
        if( itype .eq. 29 ) val=         3*(x-x0)*(y-y0)*z0/r**5
        if( itype .eq. 30 ) val=-z0/r**3+3*(y-y0)*(y-y0)*z0/r**5
c
        if( itype .eq. 31 ) val=(x-x0)*(-3/r**3+3*(x-x0)*(x-x0)/r**5)
        if( itype .eq. 32 ) val=(y-y0)*(-1/r**3+3*(x-x0)*(x-x0)/r**5)
        if( itype .eq. 33 ) val=(x-x0)*(-1/r**3+3*(y-y0)*(y-y0)/r**5)
        if( itype .eq. 34 ) val=(y-y0)*(-3/r**3+3*(y-y0)*(y-y0)/r**5)
c
        if( itype .eq. 35 ) val=(x-x0)*(-1/r**3+3*z0*z0/r**5)
        if( itype .eq. 36 ) val=(y-y0)*(-1/r**3+3*z0*z0/r**5)
c
        if( itype .eq. 37 ) val=r

        if( itype .eq. 38 )
     $     val= -3/r**3+18*(x-x0)**2/r**5-15*(x-x0)**4/r**7
        if( itype .eq. 39 )
     $     val= 9*(x-x0)*(y-y0)/r**5-15*(x-x0)**3*(y-y0)/r**7
        if( itype .eq. 40 )
     $     val= -1/r**3+3*((x-x0)**2+(y-y0)**2)/r**5
     $          -15*(x-x0)**2*(y-y0)**2/r**7
        if( itype .eq. 41 )
     $     val= 9*(x-x0)*(y-y0)/r**5-15*(x-x0)*(y-y0)**3/r**7
        if( itype .eq. 42 )
     $     val= -3/r**3+18*(y-y0)**2/r**5-15*(y-y0)**4/r**7

        if( itype .eq. 43 )
     $     val= (9*(x-x0)*z0/r**5-15*(x-x0)**3*z0/r**7)
        if( itype .eq. 44 )
     $     val= (3*(y-y0)*z0/r**5-15*(x-x0)**2*(y-y0)*z0/r**7)
        if( itype .eq. 45 )
     $     val= (3*(x-x0)*z0/r**5-15*(x-x0)*(y-y0)**2*z0/r**7)
        if( itype .eq. 46 )
     $     val= (9*(y-y0)*z0/r**5-15*(y-y0)**3*z0/r**7)

        if( itype .eq. 47 )
     $     val= -1/r**3+3*((x-x0)**2+z0**2)/r**5
     $     -15*(x-x0)**2*z0**2/r**7
        if( itype .eq. 48 )
     $     val= 3*(x-x0)*(y-y0)/r**5-15*(x-x0)*(y-y0)*z0**2/r**7
        if( itype .eq. 49 )
     $     val= -1/r**3+3*((y-y0)**2+z0**2)/r**5
     $     -15*(y-y0)**2*z0**2/r**7

        if( itype .eq. 50 )
     $     val= (9*(x-x0)*z0/r**5-15*(x-x0)*z0**3/r**7)
        if( itype .eq. 51 )
     $     val= (9*(y-y0)*z0/r**5-15*(y-y0)*z0**3/r**7)

        return
        end
c 
c 
c 
c
c
        subroutine triaevalp
     $     (itype,rnodes,whts,nnodes,x0,y0,z0,rval)
        implicit none
        integer itype,nnodes,i
        real *8 rnodes(3,nnodes),whts(nnodes)
        real *8 x0,y0,z0,rval
        real *8 d,dx,dy,dz,r,val,x,y,z
c
c       This subtroutine returns values of constant density potentials
c       on flat triangles: used to derive the values of potentials of
c       1/r and the derivatives of 1/r up to 3-th order, r and the
c       derivatives r up to 4-th order. See triaefun for the detailed
c       description of kernels, and triartable, triabtable for the
c       mapping of different types into the corresponding r and 1/r
c       value and derivative tables.  The source triangle is assumed to
c       be in R^3.
c
c       The evaluation point is assumed to be in the far field.
c
c       All integrals are evaluated via a quadrature formula with user
c       provided nodes. This is much faster than triahquad for small
c       number of quadratures nodes assuming that that the target point
c       is separated by at least 3 or 4 source triangle diameters from
c       the source centroid.
c
c       The loops for itype need to be unrolled in order for this to be
c       efficient, also, note the sign flip in dz=(z0-z)...
c
c       Input parameters: 
c
c       itype - type of the integral to be evaluated, see below 
c              must be in range (1..51), see also triahquad
c
c       rnodes(3,nnodes) - user provided quadrature nodes inside the
c       source triangle in R^3
c
c       whts(nnodes) - quadrature weights
c       nnodes - number of quadrature nodes
c
c       x0,y0,y0 - target evaluation point
c
c       Output parameters:
c
c       rval - integral value (real*8)
c       
c
c
c       ... unroll loops
c
        if( itype .eq. 1 ) then
        d = 0
        do i=1,nnodes 
          dx=rnodes(1,i)-x0
          dy=rnodes(2,i)-y0
          dz=z0-rnodes(3,i)
          r=sqrt(dx**2+dy**2+dz**2)          
          val=1/r
c          call triaefun(itype,x,y,z,x0,y0,z0,val)
          d = d + whts(i)*val
        enddo
        rval = d
        return
        endif

        if( itype .ge. 2 .and. itype .le. 4 ) then
        d = 0
        do i=1,nnodes 
          dx=rnodes(1,i)-x0
          dy=rnodes(2,i)-y0
          dz=z0-rnodes(3,i)
          r=sqrt(dx**2+dy**2+dz**2)          
          if( itype .eq. 2 ) val=dx/r**3
          if( itype .eq. 3 ) val=dy/r**3
          if( itype .eq. 4 ) val=dz/r**3
c          call triaefun(itype,x,y,z,x0,y0,z0,val)
          d = d + whts(i)*val
        enddo
        rval = d
        return
        endif

        if( itype .ge. 5 .and. itype .le. 7 ) then
        d = 0
        do i=1,nnodes 
          dx=rnodes(1,i)-x0
          dy=rnodes(2,i)-y0
          dz=z0-rnodes(3,i)
          r=sqrt(dx**2+dy**2+dz**2)          
          if( itype .eq. 5 ) val=3*dx*dz/r**5
          if( itype .eq. 6 ) val=3*dy*dz/r**5
          if( itype .eq. 7 ) val=-1/r**3+3*dz*dz/r**5
c          call triaefun(itype,x,y,z,x0,y0,z0,val)
          d = d + whts(i)*val
        enddo
        rval = d
        return
        endif
c
c
c        if( itype .eq. 7 ) then
c        d = 0
c        do i=1,nnodes 
c          dx=rnodes(1,i)-x0
c          dy=rnodes(2,i)-y0
c          dz=z0-rnodes(3,i)
c          r=sqrt(dx**2+dy**2+dz**2)          
c          val=(-r**2+3*dz*dz)/r**5
c          call triaefun(itype,x,y,z,x0,y0,z0,val)
c          d = d + whts(i)*val
c        enddo
c        rval = d
c        return
c        endif
c
c
c       very slow due to function overhead ...
c
        if( itype .ge. 8 .and. itype .le. 51 ) then
        d = 0
        do i=1,nnodes 
          x=rnodes(1,i)
          y=rnodes(2,i)
          z=rnodes(3,i)
          call triaefun(itype,x,y,z,x0,y0,z0,val)
          d = d + whts(i)*val
        enddo
        rval = d
        return
        endif

        return
        end
c
c
c
c
        subroutine triaefun(itype,x,y,z,x0,y0,z0,val)
        implicit none
        integer itype
        real *8 x,y,z,x0,y0,z0,val
        real *8 dx,dy,dz,r
c
c       This subroutine returns the values of kernels to be integrated.
c       Integration domain is assumed to be in R^3: (x,y,z) 
c
c       Input parameters: 
c
c       itype - type of the integral to be evaluated, see below 
c              must be in range (1..51)
c
c       x,y,z - source evaluation point
c
c       x0,y0,y0 - target evaluation point
c
c       Output parameters:
c
c       val - kernel value (real*8)
c
c       Note the sign flip in dz...
c       
        val=0
c
        dx=x-x0
        dy=y-y0
        dz=z0-z
        r=sqrt(dx**2+dy**2+dz**2)
c
        if( itype .eq. 1 ) val=1/r
        if( itype .eq. 2 ) val=dx/r**3
        if( itype .eq. 3 ) val=dy/r**3
        if( itype .eq. 4 ) val=dz/r**3
        if( itype .eq. 5 ) val=3*dx*dz/r**5
        if( itype .eq. 6 ) val=3*dy*dz/r**5
        if( itype .eq. 7 ) val=-1/r**3+3*dz*dz/r**5
c
        if( itype .eq. 8 ) val=dx/r
        if( itype .eq. 9 ) val=dy/r
c
        if( itype .eq. 10 ) val=dx*dz/r**3
        if( itype .eq. 11 ) val=dy*dz/r**3
c
        if( itype .eq. 12 ) val=-1/r**3+3*dx*dx/r**5
        if( itype .eq. 13 ) val=    3*dx*dy/r**5
        if( itype .eq. 14 ) val=-1/r**3+3*dy*dy/r**5
c
        if( itype .eq. 15 ) val=3*dz/r**5-15*dx*dx*dz/r**7
        if( itype .eq. 16 ) val=         -15*dx*dy*dz/r**7
        if( itype .eq. 17 ) val=3*dz/r**5-15*dy*dy*dz/r**7
c
        if( itype .eq. 18 ) val=3*dx/r**5-15*dx*dz*dz/r**7
        if( itype .eq. 19 ) val=3*dy/r**5-15*dy*dz*dz/r**7
        if( itype .eq. 20 ) val=9*dz/r**5-15*dz*dz*dz/r**7
c
        if( itype .eq. 21 ) val=9*dx/r**5-15*dx*dx*dx/r**7
        if( itype .eq. 22 ) val=3*dy/r**5-15*dx*dx*dy/r**7
        if( itype .eq. 23 ) val=3*dx/r**5-15*dx*dy*dy/r**7
        if( itype .eq. 24 ) val=9*dy/r**5-15*dy*dy*dy/r**7
c
        if( itype .eq. 25 ) val=-1/r+dx*dx/r**3
        if( itype .eq. 26 ) val=dx*dy/r**3
        if( itype .eq. 27 ) val=-1/r+dy*dy/r**3
c
        if( itype .eq. 28 ) val=-dz/r**3+3*dx*dx*dz/r**5
        if( itype .eq. 29 ) val=         3*dx*dy*dz/r**5
        if( itype .eq. 30 ) val=-dz/r**3+3*dy*dy*dz/r**5
c
        if( itype .eq. 31 ) val=dx*(-3/r**3+3*dx*dx/r**5)
        if( itype .eq. 32 ) val=dy*(-1/r**3+3*dx*dx/r**5)
        if( itype .eq. 33 ) val=dx*(-1/r**3+3*dy*dy/r**5)
        if( itype .eq. 34 ) val=dy*(-3/r**3+3*dy*dy/r**5)
c
        if( itype .eq. 35 ) val=dx*(-1/r**3+3*dz*dz/r**5)
        if( itype .eq. 36 ) val=dy*(-1/r**3+3*dz*dz/r**5)
c
        if( itype .eq. 37 ) val=r

        if( itype .eq. 38 )
     $     val= -3/r**3+18*dx**2/r**5-15*dx**4/r**7
        if( itype .eq. 39 )
     $     val= 9*dx*dy/r**5-15*dx**3*dy/r**7
        if( itype .eq. 40 )
     $     val= -1/r**3+3*(dx**2+dy**2)/r**5
     $          -15*dx**2*dy**2/r**7
        if( itype .eq. 41 )
     $     val= 9*dx*dy/r**5-15*dx*dy**3/r**7
        if( itype .eq. 42 )
     $     val= -3/r**3+18*dy**2/r**5-15*dy**4/r**7

        if( itype .eq. 43 )
     $     val= (9*dx*dz/r**5-15*dx**3*dz/r**7)
        if( itype .eq. 44 )
     $     val= (3*dy*dz/r**5-15*dx**2*dy*dz/r**7)
        if( itype .eq. 45 )
     $     val= (3*dx*dz/r**5-15*dx*dy**2*dz/r**7)
        if( itype .eq. 46 )
     $     val= (9*dy*dz/r**5-15*dy**3*dz/r**7)

        if( itype .eq. 47 )
     $     val= -1/r**3+3*(dx**2+dz**2)/r**5
     $          -15*dx**2*dz**2/r**7
        if( itype .eq. 48 )
     $     val= 3*dx*dy/r**5-15*dx*dy*dz**2/r**7
        if( itype .eq. 49 )
     $     val= -1/r**3+3*(dy**2+dz**2)/r**5
     $          -15*dy**2*dz**2/r**7

        if( itype .eq. 50 ) val= (9*dx*dz/r**5-15*dx*dz**3/r**7)
        if( itype .eq. 51 ) val= (9*dy*dz/r**5-15*dy*dz**3/r**7)

        return
        end
c 
c 
c 
c
c
        subroutine triahproj(vert0,vert1,vert2,h,a,b,verth)
        implicit none
        real *8 vert0(2),vert1(2),vert2(2),verth(2)
        real *8 a,b,h
        real *8 v0(2),v1(2),v2(2)
        real *8 w(12)
c
c       Find the signed length of height and projections 
c       of vertex vert0 to segment defined by vertices vert1 and vert2
c       
c       Input parameter:s
c
c       vert0, vert1, vert2 - triangle vertices
c       
c                C vert0
c                .   
c               / \   
c              / | \
c             /  |  \
c            /   |   \
c           /    |h   \      
c          /     |     \
c         /   a  |  b   \
c        /_______________\
c       A        H        B
c       vert1   verth     vert2
c                 
c                  
c       
c
c       Output parameters:
c
c       h = length of  CH (height)
c       a = length of  AH 
c       b = length of  HB 
c       verth = projection point on segment (vert1,vert2)
c
c
c
        call trianini(vert0,vert1,vert2,w)
ccc        call trianfor(w,vert1,v1)
        call trianfor(w,vert2,v2)
        call trianfor(w,vert0,v0)
        h=v0(2)
        a=v0(1)
        b=v2(1)-v0(1)
c 
        v0(2)=0
        call trianbak(w,v0,verth)
c
        return
        end
c
c
c
c
c
        subroutine trialquad(itype,iquad,a,b,y,z,rval)
        implicit none
        integer itype,iquad,iftaylor
        real *8 vert1(2),vert2(2)
        real *8 a,b,y,z,rval,c,hval
c
c       This subroutine integrates various Hilbert kernels
c
c
        rval=0
c
        c=sqrt(y**2+z**2)
c
c
c
c       ... tangential derivatives of the single layer
c
        if( itype .eq. 2 .or. itype .eq. 3 ) then
c
c       \int_a^b 1/sqrt(x^2+y^2+z^2) dx
c
        iftaylor = 0
        if( abs(c) .lt. 1d-8*abs(a) ) iftaylor=1
        if( abs(c) .lt. 1d-8*abs(b) ) iftaylor=1
        if( iftaylor .eq. 1 ) then
cc        if( a .gt. 0 .and. b .gt. 0 ) rval = log(b)-log(a)
cc        if( a .lt. 0 .and. b .lt. 0 ) rval = -log(-b)+log(-a)
cc        if( a .gt. 0 .and. b .gt. 0 ) rval = log(b/a) 
cc        if( a .lt. 0 .and. b .lt. 0 ) rval = log(a/b) 
        if( a .gt. 0 .and. b .gt. 0 ) 
     $     rval = +(log(b/a) + c**2*(1/b**2-1/a**2)/4)
        if( a .lt. 0 .and. b .lt. 0 ) 
     $     rval = -(log(b/a) + c**2*(1/b**2-1/a**2)/4)
        if( a .gt. 0 .and. b .lt. 0 ) rval = 0
        if( a .lt. 0 .and. b .gt. 0 ) rval = 0
        else        
ccc        rval = log(b+sqrt(c**2+b**2)) - log(a+sqrt(c**2+a**2))
ccc        rval = log((b+sqrt(c**2+b**2))/(a+sqrt(c**2+a**2))) 
        call triahquad_sing1(a,b,c,rval)
        endif
c
ccc        write(*,*) '1==>', a, b, c, rval, iftaylor
ccc        write(*,*) '2==>', b+sqrt(c**2+b**2)
ccc        write(*,*) '3==>', a+sqrt(c**2+a**2)
c
        endif
c
c
c
c       ... tangential derivatives of the double layer
c
        if( itype .eq. 5 .or. itype .eq. 6 ) then
c
c       \int_a^b z/sqrt(x^2+y^2+z^2)^3 dx
c
        iftaylor = 0
        if( abs(c) .lt. 1d-8*abs(a) ) iftaylor=1
        if( abs(c) .lt. 1d-8*abs(b) ) iftaylor=1
        if( iftaylor .eq. 1 ) then
        if( a .gt. 0 .and. b .gt. 0 ) rval = +z*(1/a**2-1/b**2)/2
        if( a .lt. 0 .and. b .lt. 0 ) rval = -z*(1/a**2-1/b**2)/2
        if( a .gt. 0 .and. b .lt. 0 ) rval = 0
        if( a .lt. 0 .and. b .gt. 0 ) rval = 0
        else
ccc        rval = z*b/c**2/sqrt(c**2+b**2) - z*a/c**2/sqrt(c**2+a**2)
        call triahquad_sing7(a,b,c,rval)
        rval = z*rval
        endif
c
cc        write(*,*) '1==>', a, b, c, rval, iftaylor
cc        write(*,*) '2==>', b+sqrt(c**2+b**2)
cc        write(*,*) '3==>', a+sqrt(c**2+a**2)
c
        endif
c
c
c
c       ... linear moments of the single layer
c
        if( itype .eq. 8 .or. itype .eq. 9 ) then
c
c> int(x*(1/r),x);  
c                                 2    2    2 1/2
c                               (x  + y  + z )
c
c       \int_a^b sqrt(x^2+y^2+z^2) dx
c
c> int(int(x*(1/r),x),x);
c        2    2    2 1/2                2    2    2 1/2   2
c1/2 x (x  + y  + z )    + 1/2 ln(x + (x  + y  + z )   ) y
c
c                    2    2    2 1/2   2
c     + 1/2 ln(x + (x  + y  + z )   ) z
c
        iftaylor = 0
        if( abs(c) .lt. 1d-8*abs(a) ) iftaylor=1
        if( abs(c) .lt. 1d-8*abs(b) ) iftaylor=1
        if( iftaylor .eq. 1 ) then
        if( a .gt. 0 .and. b .gt. 0 ) 
     $     rval = -(b-a)*(a+b)/2 
     $     -0.5d0*c**2 * (log(b/a) + c**2*(1/b**2-1/a**2)/4)
        if( a .lt. 0 .and. b .lt. 0 ) 
     $     rval = +(b-a)*(a+b)/2 
     $     +0.5d0*c**2 * (log(b/a) + c**2*(1/b**2-1/a**2)/4)
        if( a .gt. 0 .and. b .lt. 0 ) rval = 0
        if( a .lt. 0 .and. b .gt. 0 ) rval = 0
        else        
cc        rval = -0.5d0*(b*sqrt(b**2+c**2)+log(b+sqrt(b**2+c**2))*c**2)
cc     $     +0.5d0*(a*sqrt(a**2+c**2)+log(a+sqrt(a**2+c**2))*c**2)

c        rval = 
c     $     -0.5d0*b*sqrt(b**2+c**2)-0.5d0*log(b+sqrt(b**2+c**2))*c**2
c     $     +0.5d0*a*sqrt(a**2+c**2)+0.5d0*log(a+sqrt(a**2+c**2))*c**2

c        rval = 
c     $     -0.5d0*b*sqrt(b**2+c**2)
c     $     +0.5d0*a*sqrt(a**2+c**2)
        call triahquad_sing_aux1(a,b,c,hval)
        rval = -.5d0*hval
        call triahquad_sing1(a,b,c,hval)
        rval = rval -.5d0*c**2*hval 

        endif
c
        endif
c
c
c       ... linear moments of the double layer
c       
        if( itype .eq. 10 .or. itype .eq. 11 ) then
c
c> int(x*(z/r^3),x);           
c                                        z
c                              - -----------------
c                                  2    2    2 1/2
c                                (x  + y  + z )
c
c       \int_a^b z/sqrt(x^2+y^2+z^2) dx
c
c> int(int(x*(z/r^3),x),x);
c                                     2    2    2 1/2
c                         -z ln(x + (x  + y  + z )   )
c
        iftaylor = 0
        if( abs(c) .lt. 1d-8*abs(a) ) iftaylor=1
        if( abs(c) .lt. 1d-8*abs(b) ) iftaylor=1
        if( iftaylor .eq. 1 ) then
        if( a .gt. 0 .and. b .gt. 0 ) 
     $     rval = +z*(log(b/a) + c**2*(1/b**2-1/a**2)/4)
        if( a .lt. 0 .and. b .lt. 0 ) 
     $     rval = -z*(log(b/a) + c**2*(1/b**2-1/a**2)/4)
        if( a .gt. 0 .and. b .lt. 0 ) rval = 0
        if( a .lt. 0 .and. b .gt. 0 ) rval = 0
        else        
ccc        rval = z*(log(b+sqrt(c**2+b**2)) - log(a+sqrt(c**2+a**2)))
        call triahquad_sing1(a,b,c,rval)
        rval = z*rval
        endif
c
        endif
c
c
c       ... tangential derivatives of the quadruple layer
c
        if( itype .eq. 18 .or. itype .eq. 19 ) then
c
c       \int_a^b diff(z/sqrt(x^2+y^2+z^2)^3,z) dx
c
        iftaylor = 0
        if( abs(c) .lt. 1d-8*abs(a) ) iftaylor=1
        if( abs(c) .lt. 1d-8*abs(b) ) iftaylor=1
        if( iftaylor .eq. 1 ) then
        if( a .gt. 0 .and. b .gt. 0 ) 
     $     rval = +.5d0*(1/a**2-1/b**2)
ccc     $            -.75d0*z**2*(1/a**4-1/b**4)
        if( a .lt. 0 .and. b .lt. 0 )
     $     rval = -.5d0*(1/a**2-1/b**2)
ccc     $            +.75d0*z**2*(1/a**4-1/b**4)
        if( a .gt. 0 .and. b .lt. 0 ) rval = 0
        if( a .lt. 0 .and. b .gt. 0 ) rval = 0
        else        
c
c        rval = 
c     $     -b/sqrt(c**2+b**2)/c**2*
c     $     (z**2*(1/(c**2+b**2)+2/c**2)-1) 
c     $     +a/sqrt(c**2+a**2)/c**2*
c     $     (z**2*(1/(c**2+a**2)+2/c**2)-1)
c 
c        rval = 
c     $     +b/sqrt(c**2+b**2)/c**2
c     $     -a/sqrt(c**2+a**2)/c**2
c     $     -b/sqrt(c**2+b**2)/c**2*
c     $     (z**2*(1/(c**2+b**2)+2/c**2))
c     $     +a/sqrt(c**2+a**2)/c**2*
c     $     (z**2*(1/(c**2+a**2)+2/c**2))
c
        call triahquad_sing7(a,b,c,hval)
        rval = hval
        call triahquad_sing20(a,b,c,hval)
        rval = rval - hval*z*z

        endif
c
ccc        write(*,*) 'rval=*', rval, iftaylor, a, b
c
        endif
c
c
c
c       ... linear moments of the normal derivative of double layer
c
        if( itype .eq. 35 .or. itype .eq. 36 ) then
c
c> int(x*(-1/r^3+3*z*z/r^5),x);
c                              2
c                             z                    1
c                    - ----------------- + -----------------
c                        2    2    2 3/2     2    2    2 1/2
c                      (x  + y  + z )      (x  + y  + z )
c
c       \int_a^b 1/sqrt(x^2+y^2+z^2)-z*z/sqrt(x^2+y^2+z^2)^3 dx
c
c> map(simplify,int(int(x*(-1/r^3+3*z*z/r^5),x),x));
c                         2
c                        z  x                        2    2    2 1/2
c           - --------------------------- + ln(x + (x  + y  + z )   )
c               2    2    2    2    2 1/2
c             (y  + z ) (x  + y  + z )
c
        iftaylor = 0
        if( abs(c) .lt. 1d-8*abs(a) ) iftaylor=1
        if( abs(c) .lt. 1d-8*abs(b) ) iftaylor=1
        if( iftaylor .eq. 1 ) then
        if( a .gt. 0 .and. b .gt. 0 ) 
     $     rval = -(log(b/a) + c**2*(1/b**2-1/a**2)/4)
     $            +z**2/2*(1/a**2-1/b**2)
        if( a .lt. 0 .and. b .lt. 0 )
     $     rval = +(log(b/a) + c**2*(1/b**2-1/a**2)/4)
     $            -z**2/2*(1/a**2-1/b**2)
        if( a .gt. 0 .and. b .lt. 0 ) rval = 0
        if( a .lt. 0 .and. b .gt. 0 ) rval = 0
        else
c        x=b
c        r=sqrt(x**2+c**2)
c        fb=log(x+r)-z**2*x/(c**2)/r
c        x=a
c        r=sqrt(x**2+c**2)
c        fa=log(x+r)-z**2*x/(c**2)/r
c        rval = fb-fa
c        rval = -rval

        call triahquad_sing1(a,b,c,hval)
        rval= -hval 
        call triahquad_sing7(a,b,c,hval)
        rval = rval + z**2*hval 

        endif
c
ccc        write(*,*) 'rval=*', rval, a, b, c, iftaylor
c
        endif
c
c
        if( itype .eq. 50 .or. itype .eq. 51 ) then 
c
c
c> simplify(int(9*x*z/r^5-15*x*z^3/r^7),x));
c> int(9*x*z/r^5-15*x*z^3/r^7,x);          
c                             3
c                            z                      z
c                   3 ----------------- - 3 -----------------
c                       2    2    2 5/2       2    2    2 3/2
c                     (x  + y  + z )        (x  + y  + z )
c
c> simplify(int(1/r^3,x)); 
c                                       x
c                          ---------------------------
c                            2    2    2    2    2 1/2
c                          (y  + z ) (x  + y  + z )
c
c> simplify(int(1/r^5,x)); 
c
c                                    2      2      2
c                              x (3 y  + 3 z  + 2 x )
c                       1/3 ----------------------------
c                             2    2 2   2    2    2 3/2
c                           (y  + z )  (x  + y  + z )
c
c
        iftaylor = 0
        if( abs(c) .lt. 1d-8*abs(a) ) iftaylor=1
        if( abs(c) .lt. 1d-8*abs(b) ) iftaylor=1
        if( iftaylor .eq. 1 ) then
        if( a .gt. 0 .and. b .gt. 0 ) 
     $     rval = -1.5d0*z*(1/b**2-1/a**2)
ccc     $            +.75d0*z**3*(1/b**4-1/a**4) 
ccc     $            -1.25d0*z**3*(1/b**6-1/a**6)*c**2
        if( a .lt. 0 .and. b .lt. 0 )
     $     rval = +1.5d0*z*(1/b**2-1/a**2)
ccc     $            -.75d0*z**3*(1/b**4-1/a**4)
ccc     $            +1.25d0*z**3*(1/b**6-1/a**6)*c**2
        if( a .gt. 0 .and. b .lt. 0 ) rval = 0
        if( a .lt. 0 .and. b .gt. 0 ) rval = 0
        else
c        x=b
c        r=sqrt(x**2+c**2)
c        fb=z**3*x*(3*c**2+2*x*x)/(c**2)**2/r**3
c        fb=fb-3*z*x/(c**2)/r
c        x=a
c        r=sqrt(x**2+c**2)
c        fa=z**3*x*(3*c**2+2*x*x)/(c**2)**2/r**3
c        fa=fa-3*z*x/(c**2)/r
c        rval = fb-fa
c        rval = -rval

        call triahquad_sing7(a,b,c,hval)
        rval= 3*z*hval 
        call triahquad_sing20(a,b,c,hval)
        rval = rval - z**3*hval 

        endif
c
ccc        write(*,*) 'rval=*', rval, a, b, c, iftaylor
c
        endif
c
c
c       ... evaluation point is too close to the line, abort
c        
c
c        eps=1d-12
c        if( abs(c) .lt. eps*abs(b) ) return
c        if( abs(c) .lt. eps*abs(a) ) return
c
c
c
        return
        end
c
c
c
c
c  
        subroutine trialquad_hd(itype,iquad,cosa,sina,a,b,y,z,rval)
        implicit none
        integer itype,iquad,iftaylor
        real *8 vert1(2),vert2(2)
        real *8 cosa,sina,a,b,y,z,rval
        real *8 c,fa,fb,hval,r,x
c
c       This subroutine integrates various Hilbert kernels
c
c
        rval=0
c
        c=sqrt(y**2+z**2)
c
c       ... tangential derivatives of the modified linear moments of 
c       hilbert integral (single layer) 
c
        if( itype .eq. 25 .or. itype .eq. 26 .or. itype .eq. 27 ) then
c
c> (int(-1/r+x*x/r^3,x));        
c                          x         
c                - ----------------- 
c                    2    2    2 1/2
c                  (x  + y  + z )
c
c       \int_a^b [x cos(alpha) + y sin(alpha)]/sqrt(x^2+y^2+z^2) dx
c
c> int(x/r,x);
c                                 2    2    2 1/2
c                               (x  + y  + z )
c> int(y/r,x);
c                                     2    2    2 1/2
c                          y ln(x + (x  + y  + z )   )
c
        iftaylor = 0
        if( abs(c) .lt. 1d-8*abs(a) ) iftaylor=1
        if( abs(c) .lt. 1d-8*abs(b) ) iftaylor=1
        if( iftaylor .eq. 1 ) then
        if( a .gt. 0 .and. b .gt. 0 ) 
     $    rval = -((b-a)+0.5d0*(1/b-1/a)*c**2) *cosa
     $     -(log(b/a) + c**2*(1/b**2-1/a**2)/4)*y*sina
        if( a .lt. 0 .and. b .lt. 0 ) 
     $    rval = +((b-a)+0.5d0*(1/b-1/a)*c**2) *cosa
     $     +(log(b/a) + c**2*(1/b**2-1/a**2)/4)*y*sina
        if( a .gt. 0 .and. b .lt. 0 ) rval = 0
        if( a .lt. 0 .and. b .gt. 0 ) rval = 0
        else        
c        x=b
c        r=sqrt(x**2+c**2)
c        fb=r*cosa+y*log(x+r)*sina        
c        x=a
c        r=sqrt(x**2+c**2)
c        fa=r*cosa+y*log(x+r)*sina
c        rval = fb-fa
c        rval = -rval

ccc        rval = (b-a)*(b+a)/(sqrt(b**2+c**2)+sqrt(a**2+c**2))*cosa
        rval = (sqrt(b**2+c**2)-sqrt(a**2+c**2))*cosa 
        call triahquad_sing1(a,b,c,hval)
        rval= rval + hval*y*sina
        rval = -rval

        endif
c
        endif
c
c
c
c       ... tangential derivatives of the hilbert integral (single layer)
c
        if( itype .eq. 12 .or. itype .eq. 13 .or. itype .eq. 14 ) then
c
c> simplify(int(-1/r^3+3*x*x/r^5,x));
c                                        x
c                              - -----------------
c                                  2    2    2 3/2
c                                (x  + y  + z )
c
c       \int_a^b [x cos(alpha) + y sin(alpha)]/sqrt(x^2+y^2+z^2)^3 dx
c
c> simplify(int(diff(1/r,x),x));  
c                                       1
c                               -----------------
c                                 2    2    2 1/2
c                               (x  + y  + z )
c
c> simplify(int(diff(1/r,y),x));  
c                                       x y
c                         - ---------------------------
c                             2    2    2    2    2 1/2
c                           (y  + z ) (x  + y  + z )
c
        iftaylor = 0
        if( abs(c) .lt. 1d-8*abs(a) ) iftaylor=1
        if( abs(c) .lt. 1d-8*abs(b) ) iftaylor=1
        if( iftaylor .eq. 1 ) then
        if( a .gt. 0 .and. b .gt. 0 ) 
     $     rval = cosa*(1/b-1/a)+y*sina/2*(1/b**2-1/a**2)
        if( a .lt. 0 .and. b .lt. 0 ) 
     $     rval = -cosa*(1/b-1/a)-y*sina/2*(1/b**2-1/a**2)
        if( a .gt. 0 .and. b .lt. 0 ) rval = 0
        if( a .lt. 0 .and. b .gt. 0 ) rval = 0
        else        
c        x=b
c        r=sqrt(x**2+c**2)
c        fb=1/r*cosa-y*x/c**2/r*sina
c        x=a
c        r=sqrt(x**2+c**2)
c        fa=1/r*cosa-y*x/c**2/r*sina
c        rval = fb-fa        

        x=b
        r=sqrt(x**2+c**2)
        fb=1/r*cosa
        x=a
        r=sqrt(x**2+c**2)
        fa=1/r*cosa
        rval = fb-fa        
        call triahquad_sing7(a,b,c,hval)
        rval = rval - y*sina*hval

        endif
c
ccc        write(*,*) '===>', sina, cosa, c, rval
c       
        endif
c
c
        if( itype .eq. 15 .or. itype .eq. 16 .or. itype .eq. 17 ) then
c
c       ... tangential derivatives of the hilbert integral (double layer)
c
c       \int_a^b 3*[x cos(alpha) + y sin(alpha)]*z/sqrt(x^2+y^2+z^2)^5 dx
c
c> simplify(int(diff(1/r,x,z),x));
c                                       z
c                              - -----------------
c                                  2    2    2 3/2
c                                (x  + y  + z )
c
c> simplify(int(diff(1/r,y,z),x));
c                                    2      2      2
c                          x z y (3 y  + 3 z  + 2 x )
c                         ----------------------------
c                           2    2 2   2    2    2 3/2
c                         (y  + z )  (x  + y  + z )
c
        iftaylor = 0
        if( abs(c) .lt. 1d-8*abs(a) ) iftaylor=1
        if( abs(c) .lt. 1d-8*abs(b) ) iftaylor=1
        if( iftaylor .eq. 1 ) then
        if( a .gt. 0 .and. b .gt. 0 ) 
     $   rval = -z*cosa*(1/b**3-1/a**3) 
ccc     $          -z*y*sina*.75d0*(1/b**4-1/a**4) 
        if( a .lt. 0 .and. b .lt. 0 ) 
     $   rval = +z*cosa*(1/b**3-1/a**3) 
ccc     $          +z*y*sina*.75d0*(1/b**4-1/a**4) 
        if( a .gt. 0 .and. b .lt. 0 ) rval = 0
        if( a .lt. 0 .and. b .gt. 0 ) rval = 0
        else        
c        x=b
c        r=sqrt(x**2+c**2)
c        fb=-z/r**3*cosa
c     $     +z*y*x*(3*c**2+2*x**2)/c**4/r**3*sina
c        x=a
c        r=sqrt(x**2+c**2)
c        fa=-z/r**3*cosa
c     $     +z*y*x*(3*c**2+2*x**2)/c**4/r**3*sina
c        rval = fb-fa

        x=b
        r=sqrt(x**2+c**2)
        fb=-z/r**3*cosa
        x=a
        r=sqrt(x**2+c**2)
        fa=-z/r**3*cosa
        rval = fb-fa
        call triahquad_sing20(a,b,c,hval)
        rval = rval + z*y*hval*sina

        endif
c
        endif
c
c
c
        if( itype .eq. 21 .or. itype .eq. 22 
     $     .or. itype .eq. 23 .or. itype .eq. 24) then
c
c       +\int_a^b 3*[x cos(alpha) + y sin(alpha)]**2/sqrt(x^2+y^2+z^2)^5 dx
c       -\int_a^b 1/sqrt(x^2+y^2+z^2)^3 dx
c
c> simplify(int(3*(x*cosa)**2/r^5,x));       
c                                    3     2
c                                   x  cosa
c                          ---------------------------
c                            2    2    2    2    2 3/2
c                          (y  + z ) (x  + y  + z )
c
c> simplify(int(3*(x*cosa)*(y*sina)*2/r^5,x));
c                                   y sina cosa
c                             -2 -----------------
c                                  2    2    2 3/2
c                                (x  + y  + z )
c
c> simplify(int(3*(y*sina)**2/r^5,x));        
c                         2     2       2      2      2
c                        y  sina  x (3 y  + 3 z  + 2 x )
c                        -------------------------------
c                           2    2 2   2    2    2 3/2
c                         (y  + z )  (x  + y  + z )
c
c> simplify(int(1/r^3,x));
c                                       x
c                          ---------------------------
c                            2    2    2    2    2 1/2
c                          (y  + z ) (x  + y  + z )
c
        iftaylor = 0
        if( abs(c) .lt. 1d-8*abs(a) ) iftaylor=1
        if( abs(c) .lt. 1d-8*abs(b) ) iftaylor=1
        if( iftaylor .eq. 1 ) then
        if( a .gt. 0 .and. b .gt. 0 ) 
     $   rval = +1.5d0*cosa**2*(1/b**2-1/a**2) 
     $          -0.5d0*(1/b**2-1/a**2)
     $          +2*y*cosa*sina*(1/b**3-1/a**3) 
ccc     $          +0.75d0*y**2*sina**2*(1/b**4-1/a**4) 
        if( a .lt. 0 .and. b .lt. 0 ) 
     $   rval = -1.5d0*cosa**2*(1/b**2-1/a**2) 
     $          +.5d0*(1/b**2-1/a**2)
     $          -2*y*cosa*sina*(1/b**3-1/a**3) 
ccc     $          -.75d0*y**2*sina**2*(1/b**4-1/a**4) 
        if( a .gt. 0 .and. b .lt. 0 ) rval = 0
        if( a .lt. 0 .and. b .gt. 0 ) rval = 0
        else        
c        x=b
c        r=sqrt(x**2+c**2)
c        fb=x**3*cosa**2/c**2/r**3
c        fb=fb-2*y*cosa*sina/r**3
c        fb=fb+y**2*sina**2*x*(3*c**2+2*x**2)/c**4/r**3
c        fb=fb-x/c**2/r
c        x=a
c        r=sqrt(x**2+c**2)
c        fa=x**3*cosa**2/c**2/r**3
c        fa=fa-2*y*cosa*sina/r**3
c        fa=fa+y**2*sina**2*x*(3*c**2+2*x**2)/c**4/r**3
c        fa=fa-x/c**2/r
c        rval = fb-fa
c        rval = -rval

        x=b
        r=sqrt(x**2+c**2)
        fb=-2*y*cosa*sina/r**3
        x=a
        r=sqrt(x**2+c**2)
        fa=-2*y*cosa*sina/r**3
        rval = fb-fa
        call triahquad_sing7(a,b,c,hval)
        rval = rval - hval
        call triahquad_sing20(a,b,c,hval)
        rval = rval + y**2*sina**2*hval
        call triahquad_sing20h(a,b,c,hval)
        rval = rval + cosa**2*hval
        rval = -rval
        
        endif
c
ccc        write(*,*) fa,fb
ccc        pause
ccc        write(12,*) 'rval=*',a,b,c, rval
ccc        write(12,*) iftaylor,itype,sina,cosa
c
        endif
c
c
c       ... tangential derivatives of the modified linear moments of 
c       hilbert integral (double layer) 
c
        if( itype .eq. 28 .or. itype .eq. 29 .or. itype .eq. 30 ) then
c
c> simplify((int(-z/r^3+3*x*x*z/r^5,x)));  
c                                      x z
c                              - -----------------
c                                  2    2    2 3/2
c                                (x  + y  + z )
c
c
c       \int_a^b [x cos(alpha) + y sin(alpha)]*z/sqrt(x^2+y^2+z^2)^3 dx
c
c> simplify(int(diff(1/r,x),x));  
c                                       1
c                               -----------------
c                                 2    2    2 1/2
c                               (x  + y  + z )
c
c> simplify(int(diff(1/r,y),x));  
c                                       x y
c                         - ---------------------------
c                             2    2    2    2    2 1/2
c                           (y  + z ) (x  + y  + z )
c
c
        iftaylor = 0
        if( abs(c) .lt. 1d-8*abs(a) ) iftaylor=1
        if( abs(c) .lt. 1d-8*abs(b) ) iftaylor=1
        if( iftaylor .eq. 1 ) then
        if( a .gt. 0 .and. b .gt. 0 ) 
     $   rval = +z*cosa*(1/b-1/a) +.5d0*z*y*sina*(1/b**2-1/a**2)
        if( a .lt. 0 .and. b .lt. 0 ) 
     $   rval = -z*cosa*(1/b-1/a) -.5d0*z*y*sina*(1/b**2-1/a**2)
        if( a .gt. 0 .and. b .lt. 0 ) rval = 0
        if( a .lt. 0 .and. b .gt. 0 ) rval = 0
        else        
c        x=b
c        r=sqrt(x**2+c**2)
c        fb=z/r*cosa-y*x*z/c**2/r*sina
c        x=a
c        r=sqrt(x**2+c**2)
c        fa=z/r*cosa-y*x*z/c**2/r*sina
c        rval = fb-fa

        x=b
        r=sqrt(x**2+c**2)
        fb=z/r*cosa
        x=a
        r=sqrt(x**2+c**2)
        fa=z/r*cosa
        rval = fb-fa

        call triahquad_sing7(a,b,c,hval)
        rval= rval - z*y*sina*hval         


        endif
c
        endif
c
c
c
        if( itype .eq. 31 .or. itype .eq. 34 .or.
     $     itype .eq. 32 .or. itype .eq. 33) then
c
c> simplify(int(x*(-3/r^3+3*x*x/r^5),x));
c                                      2    2
c                                     y  + z
c                                -----------------
c                                  2    2    2 3/2
c                                (x  + y  + z )
c
c       -\int_a^b [x cos(alpha) + y sin(alpha)]**2/sqrt(x^2+y^2+z^2)^3 dx
c       +\int_a^b 1/sqrt(x^2+y^2+z^2) dx
c
c> int((x*cosa)**2/r^3,x);        
c                        2
c                    cosa  x            2          2    2    2 1/2
c             - ----------------- + cosa  ln(x + (x  + y  + z )   )
c                 2    2    2 1/2
c               (x  + y  + z )
c
c> simplify(int((x*cosa)*(y*sina)*2/r^3,x));
c                                   cosa y sina
c                             -2 -----------------
c                                  2    2    2 1/2
c                                (x  + y  + z )
c
c> simplify(int((y*sina)**2/r^3,x));   
c                                    2     2
c                                   y  sina  x
c                           ---------------------------
c                             2    2    2    2    2 1/2
c                           (y  + z ) (x  + y  + z )
c
c
c> simplify(int(1/r,x));  
c                                    2    2    2 1/2
c                           ln(x + (x  + y  + z )   )
c
        iftaylor = 0
        if( abs(c) .lt. 1d-8*abs(a) ) iftaylor=1
        if( abs(c) .lt. 1d-8*abs(b) ) iftaylor=1
        if( iftaylor .eq. 1 ) then
        if( a .gt. 0 .and. b .gt. 0 ) 
     $   rval = 
     $          -2*y*cosa*sina*(1/b-1/a)
     $          -.5d0*y*y*sina**2*(1/b**2-1/a**2)
     $          -sina**2*(log(b/a)+ c**2*(1/b**2-1/a**2)/4)
        if( a .lt. 0 .and. b .lt. 0 ) 
     $   rval = 
     $          +2*y*cosa*sina*(1/b-1/a)
     $          +.5d0*y*y*sina**2*(1/b**2-1/a**2)
     $          +sina**2*(log(b/a)+ c**2*(1/b**2-1/a**2)/4)
        if( a .gt. 0 .and. b .lt. 0 ) rval = 0
        if( a .lt. 0 .and. b .gt. 0 ) rval = 0
        else        
cc        x=b
cc        r=sqrt(x**2+c**2)
cc        fb=cosa**2*(-x/r+log(x+r))
cc        fb=fb-2*y*cosa*sina/r
cc        fb=fb+y**2*sina**2*x/c**2/r
cc        fb=fb-log(x+r)
cc        x=a
cc        r=sqrt(x**2+c**2)
cc        fa=cosa**2*(-x/r+log(x+r))
cc        fa=fa-2*y*cosa*sina/r
cc        fa=fa+y**2*sina**2*x/c**2/r
cc        fa=fa-log(x+r)
cc        rval = fb-fa

        x=b
        r=sqrt(x**2+c**2)
        fb=cosa**2*(-x/r)
        fb=fb-2*y*cosa*sina/r
        x=a
        r=sqrt(x**2+c**2)
        fa=cosa**2*(-x/r)
        fa=fa-2*y*cosa*sina/r
        rval = fb-fa
        call triahquad_sing1(a,b,c,hval)
        rval=rval -hval*sina**2 
        call triahquad_sing7(a,b,c,hval)
        rval=rval +hval*sina**2 *y**2 

        endif
c
        endif
c
c
c
        if( itype .eq. 43 .or. itype .eq. 44 
     $     .or. itype .eq. 45 .or. itype .eq. 46) then
c
c> diff(r,z,x,x);  
c                               2
c                            z x                   z
c                    3 ----------------- - -----------------
c                        2    2    2 5/2     2    2    2 3/2
c                      (x  + y  + z )      (x  + y  + z )
c
c       +\int_a^b 3*[x cos(alpha) + y sin(alpha)]**2/sqrt(x^2+y^2+z^2)^5 dx
c       -\int_a^b 1/sqrt(x^2+y^2+z^2)^3 dx
c
c> simplify(int(3*(x*cosa)**2/r^5,x));       
c                                    3     2
c                                   x  cosa
c                          ---------------------------
c                            2    2    2    2    2 3/2
c                          (y  + z ) (x  + y  + z )
c
c> simplify(int(3*(x*cosa)*(y*sina)*2/r^5,x));
c                                   y sina cosa
c                             -2 -----------------
c                                  2    2    2 3/2
c                                (x  + y  + z )
c
c> simplify(int(3*(y*sina)**2/r^5,x));        
c                         2     2       2      2      2
c                        y  sina  x (3 y  + 3 z  + 2 x )
c                        -------------------------------
c                           2    2 2   2    2    2 3/2
c                         (y  + z )  (x  + y  + z )
c
c> simplify(int(1/r^3,x));
c                                       x
c                          ---------------------------
c                            2    2    2    2    2 1/2
c                          (y  + z ) (x  + y  + z )
c
        iftaylor = 0
        if( abs(c) .lt. 1d-8*abs(a) ) iftaylor=1
        if( abs(c) .lt. 1d-8*abs(b) ) iftaylor=1
        if( iftaylor .eq. 1 ) then
        if( a .gt. 0 .and. b .gt. 0 ) then
         rval = +1.5d0*cosa**2*(1/b**2-1/a**2)
     $          +2*y*cosa*sina*(1/b**3-1/a**3)
ccc     $          +.75d0*y*y*sina**2*(1/b**4-1/a**4)
     $          -.5d0*(1/b**2-1/a**2)
        rval = z*rval
        endif        
        if( a .lt. 0 .and. b .lt. 0 ) then
         rval = -1.5d0*cosa**2*(1/b**2-1/a**2)
     $          -2*y*cosa*sina*(1/b**3-1/a**3)
ccc     $          -.75d0*y*y*sina**2*(1/b**4-1/a**4)
     $          +.5d0*(1/b**2-1/a**2)
        rval = z*rval
        endif        
        if( a .gt. 0 .and. b .lt. 0 ) rval = 0
        if( a .lt. 0 .and. b .gt. 0 ) rval = 0
        else        
c        x=b
c        r=sqrt(x**2+c**2)
c        fb=x**3*cosa**2/c**2/r**3
c        fb=fb-2*y*cosa*sina/r**3
c        fb=fb+y**2*sina**2*x/c**2/r**3
c        fb=fb+y**2*sina**2*x*2/c**4/r
c        fb=fb-x/c**2/r
c        x=a
c        r=sqrt(x**2+c**2)
c        fa=x**3*cosa**2/c**2/r**3
c        fa=fa-2*y*cosa*sina/r**3
c        fa=fa+y**2*sina**2*x/c**2/r**3
c        fa=fa+y**2*sina**2*x*2/c**4/r
c        fa=fa-x/c**2/r
c        rval = fb-fa
c        rval = -rval*z

        x=b
        r=sqrt(x**2+c**2)
        fb=-2*y*cosa*sina/r**3
        x=a
        r=sqrt(x**2+c**2)
        fa=-2*y*cosa*sina/r**3
        rval = fb-fa
        call triahquad_sing7(a,b,c,hval)
        rval = rval - hval
        call triahquad_sing20(a,b,c,hval)
        rval = rval + hval* y**2 * sina**2
        call triahquad_sing20h(a,b,c,hval)
        rval= rval + cosa**2*hval         
        rval = -rval*z


        endif
c
        endif
c
c
        if( itype .eq. 47 .or. itype .eq. 48 .or. itype .eq. 49 ) then
c
c> diff(r,z,z,x);
c                             2
c                            z  x                  x
c                    3 ----------------- - -----------------
c                        2    2    2 5/2     2    2    2 3/2
c                      (x  + y  + z )      (x  + y  + z )
c
c
c       \int_a^b [x cos(alpha) + y sin(alpha)]/sqrt(x^2+y^2+z^2)^3 dx
c
c> simplify(int(diff(1/r,x),x));  
c                                       1
c                               -----------------
c                                 2    2    2 1/2
c                               (x  + y  + z )
c
c> simplify(int(diff(1/r,y),x));  
c                                       x y
c                         - ---------------------------
c                             2    2    2    2    2 1/2
c                           (y  + z ) (x  + y  + z )
c
c       \int_a^b [x cos(alpha) + y sin(alpha)]*3*z*z/sqrt(x^2+y^2+z^2)^5 dx
c
c> int(3*z*z*x/r^5,x);
c                                        2
c                                       z
c                              - -----------------
c                                  2    2    2 3/2
c                                (x  + y  + z )
c
c> int(3*z*z*y/r^5,x);
c   2   /                   x
c3 z  y |4/3 -------------------------------
c       |        2      2    2    2    2 3/2
c       \    (4 y  + 4 z ) (x  + y  + z )
c
c                           x                \
c     + 32/3 --------------------------------|
c                2      2 2   2    2    2 1/2|
c            (4 y  + 4 z )  (x  + y  + z )   /
c
        iftaylor = 0
        if( abs(c) .lt. 1d-8*abs(a) ) iftaylor=1
        if( abs(c) .lt. 1d-8*abs(b) ) iftaylor=1
        if( iftaylor .eq. 1 ) then
        if( a .gt. 0 .and. b .gt. 0 ) 
     $   rval = +cosa*(1/b-1/a)
     $          +.5d0*y*sina*(1/b**2-1/a**2)
     $          -cosa*z**2*(1/b**3-1/a**3)
ccc     $          -.75d0*z**2*sina*y*(1/b**4-1/a**4)
        if( a .lt. 0 .and. b .lt. 0 ) 
     $   rval = -cosa*(1/b-1/a)
     $          -.5d0*y*sina*(1/b**2-1/a**2)
     $          +cosa*z**2*(1/b**3-1/a**3)
ccc     $          +.75d0*z**2*sina*y*(1/b**4-1/a**4)
        if( a .gt. 0 .and. b .lt. 0 ) rval = 0
        if( a .lt. 0 .and. b .gt. 0 ) rval = 0
        else        
c        x=b
c        r=sqrt(x**2+c**2)
c        fb=1/r*cosa-y*x/c**2/r*sina
c        fb=fb-cosa*z**2/r**3
c        fb=fb+sina*z**2*y*x*(1/(c**2)/r**3+2/(c**2)**2/r)
c        x=a
c        r=sqrt(x**2+c**2)
c        fa=1/r*cosa-y*x/c**2/r*sina
c        fa=fa-cosa*z**2/r**3
c        fa=fa+sina*z**2*y*x*(1/(c**2)/r**3+2/(c**2)**2/r)
c        rval = fb-fa

        x=b
        r=sqrt(x**2+c**2)
        fb=1/r*cosa
        fb=fb-cosa*z**2/r**3
        x=a
        r=sqrt(x**2+c**2)
        fa=1/r*cosa
        fa=fa-cosa*z**2/r**3
        rval = fb-fa
c
        call triahquad_sing7(a,b,c,hval)
        rval= rval - y*sina*hval         
        call triahquad_sing20(a,b,c,hval)
        rval = rval + z**2*y*sina*hval 
c
        endif
c
        endif
c
c
        if( itype .eq. 38 .or. itype .eq. 39 
     $     .or. itype .eq. 41 .or. itype .eq. 42 ) then
c
c> diff(r,x,x,x);
c                              3
c                             x                     x
c                    3 ----------------- - 3 -----------------
c                        2    2    2 5/2       2    2    2 3/2
c                      (x  + y  + z )        (x  + y  + z )
c
c
c       \int_a^b [x cos(alpha) + y sin(alpha)]/sqrt(x^2+y^2+z^2)^3 dx
c
c> simplify(int(diff(1/r,x),x));  
c                                       1
c                               -----------------
c                                 2    2    2 1/2
c                               (x  + y  + z )
c
c> simplify(int(diff(1/r,y),x));  
c                                       x y
c                         - ---------------------------
c                             2    2    2    2    2 1/2
c                           (y  + z ) (x  + y  + z )
c
c       \int_a^b 3*[x cos(alpha) + y sin(alpha)]**3/sqrt(x^2+y^2+z^2)^3 dx
c
c> simplify(int(3*x*x*x/r^5,x));
c                                  2      2      2
c                               3 x  + 2 y  + 2 z
c                             - ------------------
c                                 2    2    2 3/2
c                               (x  + y  + z )
c
c> simplify(int(3*x*x*y/r^5,x));
c                                      3
c                                     x  y
c                          ---------------------------
c                            2    2    2    2    2 3/2
c                          (y  + z ) (x  + y  + z )
c
c> int(3*x*y*y/r^5,x);
c                                        2
c                                       y
c                              - -----------------
c                                  2    2    2 3/2
c                                (x  + y  + z )
c
c> simplify(int(3*y*y*y/r^5,x));
c                           3       2      2      2
c                          y  x (3 y  + 3 z  + 2 x )
c                         ----------------------------
c                           2    2 2   2    2    2 3/2
c                         (y  + z )  (x  + y  + z )
c
        iftaylor = 0
        if( abs(c) .lt. 1d-8*abs(a) ) iftaylor=1
        if( abs(c) .lt. 1d-8*abs(b) ) iftaylor=1
        if( iftaylor .eq. 1 ) then
        if( a .gt. 0 .and. b .gt. 0 ) 
     $   rval = +3*cosa*(1/b-1/a)
     $          +1.5d0*y*sina*(1/b**2-1/a**2)
     $          -3*cosa**3*(1/b-1/a)
     $          -4.5d0*cosa**2*sina*y*(1/b**2-1/a**2)
     $          -3*cosa*sina**2*y**2*(1/b**3-1/a**3)
ccc     $          -.75d0*sina**3*y**3*(1/b**4-1/a**4)
        if( a .lt. 0 .and. b .lt. 0 ) 
     $   rval = -3*cosa*(1/b-1/a)
     $          -1.5d0*y*sina*(1/b**2-1/a**2)
     $          +3*cosa**3*(1/b-1/a)
     $          +4.5d0*cosa**2*sina*y*(1/b**2-1/a**2)
     $          +3*cosa*sina**2*y**2*(1/b**3-1/a**3)
ccc     $          +.75d0*sina**3*y**3*(1/b**4-1/a**4)
        if( a .gt. 0 .and. b .lt. 0 ) rval = 0
        if( a .lt. 0 .and. b .gt. 0 ) rval = 0
        else        
c        x=b
c        r=sqrt(x**2+c**2)
c        fb=3/r*cosa-3*y*x/c**2/r*sina
c        fb=fb-cosa**3*(3*x*x+2*c**2)/r**3
c        fb=fb+3*cosa**2*sina*y*x**3/(c**2)/r**3
c        fb=fb-3*cosa*sina**2*y*y/r**3
c        fb=fb+sina**3*y**3*x*(3*c**2+2*x*x)/(c**2)**2/r**3
c        x=a
c        r=sqrt(x**2+c**2)
c        fa=3/r*cosa-3*y*x/c**2/r*sina
c        fa=fa-cosa**3*(3*x*x+2*c**2)/r**3
c        fa=fa+3*cosa**2*sina*y*x**3/(c**2)/r**3
c        fa=fa-3*cosa*sina**2*y*y/r**3
c        fa=fa+sina**3*y**3*x*(3*c**2+2*x*x)/(c**2)**2/r**3
c        rval = fb-fa

        x=b
        r=sqrt(x**2+c**2)
        fb=3/r*cosa
        fb=fb-cosa**3*(3*x*x+2*c**2)/r**3
        fb=fb-3*cosa*sina**2*y*y/r**3
        x=a
        r=sqrt(x**2+c**2)
        fa=3/r*cosa
        fa=fa-cosa**3*(3*x*x+2*c**2)/r**3
        fa=fa-3*cosa*sina**2*y*y/r**3
        rval = fb-fa

        call triahquad_sing7(a,b,c,hval)
        rval= rval - 3*y*sina*hval         
        call triahquad_sing20(a,b,c,hval)
        rval= rval + y**3*sina**3*hval         
        call triahquad_sing20h(a,b,c,hval)
        rval= rval + 3*cosa**2*sina*y*hval         

        endif
c
        endif
c
c
        if( itype .eq. 40 ) then
c
c> diff(r,x,y,y);
c                                2
c                             x y                   x
c                    3 ----------------- -   -----------------
c                        2    2    2 5/2       2    2    2 3/2
c                      (x  + y  + z )        (x  + y  + z )
c
c
c       \int_a^b [x cos(alpha) + y sin(alpha)]/sqrt(x^2+y^2+z^2)^3 dx
c
c> simplify(int(diff(1/r,x),x));  
c                                       1
c                               -----------------
c                                 2    2    2 1/2
c                               (x  + y  + z )
c
c> simplify(int(diff(1/r,y),x));  
c                                       x y
c                         - ---------------------------
c                             2    2    2    2    2 1/2
c                           (y  + z ) (x  + y  + z )
c
c       \int_a^b 3*[x cos(alpha) + y sin(alpha)]*[-x sin(alpha)+y cos(alpha)]^2
c              /sqrt(x^2+y^2+z^2)^3 dx
c
c> (x *ca +y *sa )*(-x *sa+y *ca)^2;
c                                                     2
c                         (x ca + y sa) (-x sa + y ca)
c
c> expand(%);
c   3      2      2   2            3  2       3  2      2   2         3      2
c  x  ca sa  - 2 x  ca  sa y + x ca  y  + y sa  x  - 2 y  sa  x ca + y  sa ca
c
c> simplify(int(3*x*x*x/r^5,x));
c                                  2      2      2
c                               3 x  + 2 y  + 2 z
c                             - ------------------
c                                 2    2    2 3/2
c                               (x  + y  + z )
c
c> simplify(int(3*x*x*y/r^5,x));
c                                      3
c                                     x  y
c                          ---------------------------
c                            2    2    2    2    2 3/2
c                          (y  + z ) (x  + y  + z )
c
c> int(3*x*y*y/r^5,x);
c                                        2
c                                       y
c                              - -----------------
c                                  2    2    2 3/2
c                                (x  + y  + z )
c
c> simplify(int(3*y*y*y/r^5,x));
c                           3       2      2      2
c                          y  x (3 y  + 3 z  + 2 x )
c                         ----------------------------
c                           2    2 2   2    2    2 3/2
c                         (y  + z )  (x  + y  + z )
c
        iftaylor = 0
        if( abs(c) .lt. 1d-8*abs(a) ) iftaylor=1
        if( abs(c) .lt. 1d-8*abs(b) ) iftaylor=1
        if( iftaylor .eq. 1 ) then
        if( a .gt. 0 .and. b .gt. 0 ) 
     $   rval = +cosa*(1/b-1/a)
     $          +0.5d0*y*sina*(1/b**2-1/a**2)
     $          -3*cosa*sina**2*(1/b-1/a)
     $          -1.5d0*(-2*cosa**2*sina+sina**3)*y*(1/b**2-1/a**2)
     $          -(cosa**3-2*cosa**sina**2)*y*y*(1/b**3-1/a**3)
ccc     $          -.75d0*cosa**2*sina*y**3*(1/b**4-1/a**4)
        if( a .lt. 0 .and. b .lt. 0 ) 
     $   rval = -cosa*(1/b-1/a)
     $          -0.5d0*y*sina*(1/b**2-1/a**2)
     $          +3*cosa*sina**2*(1/b-1/a)
     $          +1.5d0*(-2*cosa**2*sina+sina**3)*y*(1/b**2-1/a**2)
     $          +(cosa**3-2*cosa**sina**2)*y*y*(1/b**3-1/a**3)
ccc     $          +.75d0*cosa**2*sina*y**3*(1/b**4-1/a**4)
        if( a .gt. 0 .and. b .lt. 0 ) rval = 0
        if( a .lt. 0 .and. b .gt. 0 ) rval = 0
        else        
c        x=b
c        r=sqrt(x**2+c**2)
c        fb=1/r*cosa-y*x/c**2/r*sina
c        fb=fb-cosa*sina**2*(3*x*x+2*c**2)/r**3
c        fb=fb+(-2*cosa**2*sina+sina**3)*y*x**3/(c**2)/r**3
c        fb=fb-(cosa**3-2*cosa*sina**2)*y*y/r**3
c        fb=fb+
c     $     sina*cosa**2*y**3*x*(3*c**2+2*x*x)/(c**2)**2/r**3
c        x=a
c        r=sqrt(x**2+c**2)
c        fa=1/r*cosa-y*x/c**2/r*sina
c        fa=fa-cosa*sina**2*(3*x*x+2*c**2)/r**3
c        fa=fa+(-2*cosa**2*sina+sina**3)*y*x**3/(c**2)/r**3
c        fa=fa-(cosa**3-2*cosa*sina**2)*y*y/r**3
c        fa=fa+
c     $     sina*cosa**2*y**3*x*(3*c**2+2*x*x)/(c**2)**2/r**3
c        rval = fb-fa

        x=b
        r=sqrt(x**2+c**2)
        fb=1/r*cosa
        fb=fb-cosa*sina**2*(3*x*x+2*c**2)/r**3
        fb=fb-(cosa**3-2*cosa*sina**2)*y*y/r**3
        x=a
        r=sqrt(x**2+c**2)
        fa=1/r*cosa
        fa=fa-cosa*sina**2*(3*x*x+2*c**2)/r**3
        fa=fa-(cosa**3-2*cosa*sina**2)*y*y/r**3
        rval = fb-fa

        call triahquad_sing7(a,b,c,hval)
        rval= rval - y*sina*hval 
        call triahquad_sing20(a,b,c,hval)
        rval= rval + sina*cosa**2*y**3*hval 
        call triahquad_sing20h(a,b,c,hval)
        rval= rval + (-2*cosa**2*sina+sina**3)*y*hval         

        endif
c
        endif
c
c
c       ... evaluation point is too close to the line, abort
c        
c        eps=1d-12
c        if( abs(c) .lt. eps*abs(b) ) return
c        if( abs(c) .lt. eps*abs(a) ) return
c
        return
        end
c
c
c
c
c  
        subroutine triarquad(itype,iquad,a,b,z,rval)
        implicit none
        integer itype,iquad
        real *8 a,b,z,rval
        real *8 eps,val1,val2,val3,z0
c
c       This subroutine integrates constant densities on flat triangles 
c
c       The triangle is defined in R^3 and its vertices are
c
c       (0,0,0), (a,0,0), (a,b,0)
c       
c            /|b
c           / |
c          /  |
c         / T |
c        /____|
c       0     a
c
c       All integrals are evaluated at the point (0,0,z)
c
c       \int_T  f_{itype} (x,y,z) dT
c
c       itype=1:   \int_T  1/r dT
c       itype=2:   \int_T  x/r^3 dT
c       itype=3:   \int_T  y/r^3 dT
c       itype=4:   \int_T  z/r^3 dT
c       itype=5:   \int_T  3*x*z/r^5 dT
c       itype=6:   \int_T  3*y*z/r^5 dT
c       itype=7:   \int_T  -1/r^3+3*z*z/r^5 dT
c
c       itype=20:  \int_T  9*z0/r^5-15*z0*z0*z0/r^7 dT
c
c       itype=20:  \int_T  r dT
c
c       where r=1/sqrt(x^2+y^2+z^2).
c
c       iquad is the sign of z
c
c       iquad =+1, if z>0
c       iquad = 0, if z=0
c       iquad =-1, if z<0
c
c
c       diff(1/r,x)=-x/r^3
c       diff(1/r,y)=-y/r^3
c       diff(1/r,z)=-z/r^3
c
c       diff(1/r,z,x)=3*x*z/r^5
c       diff(1/r,z,y)=3*y*z/r^5
c       diff(1/r,z,z)=-1/r^3+3*z*z/r^5
c
c       diff(1/r,z,z,z)=9*z0/r^5    -15*z0*z0*z0/r^7
c
c
c>f := int(int(1/r,y=0..b/a*x),x=0..a);
c               2   2     4     2   2 1/2             2     2
c(ln(b~ a~ + (a~  b~  + a~  + a~  z~ )   ) - 1/2 ln(a~  + z~ ) - ln(a~)) a~
c  
c                               2               2
c                       a~ (I a~  + z~ a~ + I b~ )
c     - 1/2 z~ arctan(-------------------------------)
c                           2   2     4     2   2 1/2
c                     b~ (a~  b~  + a~  + a~  z~ )
c
c                                2               2
c                       a~ (-I a~  + z~ a~ - I b~ )                 a~
c     - 1/2 z~ arctan(-------------------------------) + z~ arctan(----)
c                           2   2     4     2   2 1/2               b~
c                     b~ (a~  b~  + a~  + a~  z~ )
c
c> simplify((u+v)/(1-u*v));
c
c                                      2     2     2 1/2
c                          z~ a~ b~ (b~  + a~  + z~ )
c                     -2 ---------------------------------
c                          2   2     2   2     4     2   2
c                        a~  b~  - b~  z~  + a~  + a~  z~
c
c> 
c> simplify(expand(solve(arctan(s/a)+arctan(s/a)=arctan(t),t)));
c
c                                      s a~
c                                 -2 ---------
c                                       2    2
c                                    -a~  + s
c
c       a = z*b; s = a*sqrt(a**2+b**2+z**2)
c
c
        rval=0
c
c
c       ... triangle has large aspect ratio, abort
c        
        eps=1d-14 
        if( abs(a) .lt. eps*abs(b) ) return
        if( abs(b) .lt. eps*abs(a) ) return
        if( abs(a) .eq. 0 .or. abs(b) .eq. 0 ) return
c
c
        if( itype .eq. 1 ) then 
c
        if( iquad .eq. +1 ) then 
ccc        rval=(log(b+sqrt(a**2+b**2+z**2)))*a-log(sqrt(a**2+z**2))*a
ccc     $     +z*(atan2(a,b)-atan2(a*sqrt(a**2+b**2+z**2),z*b))
        call aloga_safe1(a,z,val1)
        call aloga_safe2(a,b,z,val2)
        call atan2_safe(a*sqrt(a**2+b**2+z**2),z*b,val3)
        rval=val2-val1+z*(atan2(a,b)-val3)
        endif
c
        if( iquad .eq.  0 ) then 
ccc        write(*,*) a,b
ccc        rval=log(b+sqrt(a**2+b**2))*a-log(abs(a))*a
        z0=0.0d0
        call aloga_safe1(a,z0,val1)
        call aloga_safe2(a,b,z0,val2)
        rval=val2-val1
        endif
c
        if( iquad .eq. -1 ) then 
cc        rval=log(b+sqrt(a**2+b**2+z**2))*a-log(sqrt(a**2+z**2))*a
cc     $     +abs(z)*(atan2(a,b)-atan2(a*sqrt(a**2+b**2+z**2),abs(z)*b))
        call aloga_safe1(a,z,val1)
        call aloga_safe2(a,b,z,val2)
        call atan2_safe(a*sqrt(a**2+b**2+z**2),abs(z)*b,val3)
        rval=val2-val1+abs(z)*(atan2(a,b)-val3)
        endif
c
ccc        write(12,*) iquad, a,b,z, rval
c
        endif
c
c> fx := int(int(x/r^3,y=0..b/a*x),x=0..a);
c         2     2 1/2      2     2     2 1/2
cb~ ln((b~  + a~ )    + (b~  + a~  + z~ )   )     b~ ln(z~)
c-------------------------------------------- - --------------
c                  2     2 1/2                     2     2 1/2
c               (b~  + a~ )                     (b~  + a~ )
c
c                   2     2     2 1/2        2       2     2
c     - 1/4 ln(2 (b~  + a~  + z~ )    b~ + a~  + 2 b~  + z~ )
c
c                2       2     2        2     2     2 1/2
c     + 1/4 ln(a~  + 2 b~  + z~  - 2 (b~  + a~  + z~ )    b~)
c
        if( itype .eq. 2 ) then 
c
        if( iquad .eq. +1 ) then
        rval=b*(log(sqrt(a**2+b**2)+sqrt(a**2+b**2+z**2))-log(abs(z)))
        rval=rval/sqrt(a**2+b**2)
        rval=rval-0.25d0*log(+2*sqrt(a**2+b**2+z**2)*b+a**2+2*b**2+z**2)
        rval=rval+0.25d0*log(-2*sqrt(a**2+b**2+z**2)*b+a**2+2*b**2+z**2)
        endif
c
        if( iquad .eq.  0 ) then
c
c       log(r) singularity at this point
c        
        rval=0
        endif
c
        if( iquad .eq. -1 ) then
        rval=b*(log(sqrt(a**2+b**2)+sqrt(a**2+b**2+z**2))-log(abs(z)))
        rval=rval/sqrt(a**2+b**2)
        rval=rval-0.25d0*log(+2*sqrt(a**2+b**2+z**2)*b+a**2+2*b**2+z**2)
        rval=rval+0.25d0*log(-2*sqrt(a**2+b**2+z**2)*b+a**2+2*b**2+z**2)
        endif
c
        endif
c
c> fy := int(int(y/r^3,y=0..b/a*x),x=0..a);
c             2     2 1/2            2     2 1/2       2     2 1/2
cfy := - (-(b~  + a~ )    ln(a~ + (a~  + z~ )   ) + (b~  + a~ )    ln(z~)
c
c                2     2 1/2      2     2     2 1/2                 /
c     + a~ ln((b~  + a~ )    + (b~  + a~  + z~ )   ) - a~ ln(z~))  /
c                                                                 /
c
c       2     2 1/2
c    (b~  + a~ )
c
        if( itype .eq. 3 ) then 
c
        if( iquad .eq. +1 ) then
        rval=sqrt(b**2+a**2)*(log(abs(z))-log(a+sqrt(a**2+z**2)))
     $     +a*(log(sqrt(a**2+b**2)+sqrt(a**2+b**2+z**2))-log(abs(z)))
        rval=-rval/sqrt(a**2+b**2)
        endif
c
        if( iquad .eq.  0 ) then
c
c       log(r) singularity at this point
c        
        rval=0
        endif
c
        if( iquad .eq. -1 ) then
        rval=sqrt(b**2+a**2)*(log(abs(z))-log(a+sqrt(a**2+z**2)))
     $     +a*(log(sqrt(a**2+b**2)+sqrt(a**2+b**2+z**2))-log(abs(z)))
        rval=-rval/sqrt(a**2+b**2)
        endif
c
        endif
c
c> fz := int(int(z/r^3,y=0..b/a*x),x=0..a);
c                           2               2
c                   a~ (I a~  + z~ a~ + I b~ )
cfz := 1/2 arctan(-------------------------------)
c                       2   2     4     2   2 1/2
c                 b~ (a~  b~  + a~  + a~  z~ )
c
c                             2               2
c                    a~ (-I a~  + z~ a~ - I b~ )              a~
c     + 1/2 arctan(-------------------------------) - arctan(----)
c                        2   2     4     2   2 1/2            b~
c                  b~ (a~  b~  + a~  + a~  z~ )
c
        if( itype .eq. 4 ) then 
c
        if( iquad .eq. +1 ) then 
ccc        rval=-atan2(a,b)+atan2(a*sqrt(a**2+b**2+z**2),abs(z)*b)
        call atan2_safe(a*sqrt(a**2+b**2+z**2),abs(z)*b,val3)
        rval=-atan2(a,b)+val3
        endif
c
        if( iquad .eq.  0 ) then 
        rval=0
        endif
c
        if( iquad .eq. -1 ) then 
ccc        rval=-atan2(a,b)+atan2(a*sqrt(a**2+b**2+z**2),abs(z)*b)
        call atan2_safe(a*sqrt(a**2+b**2+z**2),abs(z)*b,val3)
        rval=-atan2(a,b)+val3
        rval=-rval
        endif
c
ccc        write(21,*) 'rval', rval
c
        endif
c
c> dy := int(int(x*z/r^5,y=0..b/a*x),x=0..a);
c                                            3
c                                          a~  b~
c     dx := - 1/3 --------------------------------------------------------
c                                                2   2     4     2   2 1/2
c                 z~ (-a~ + I z~) (a~ + I z~) (b~  a~  + a~  + a~  z~ )
c
c
        if( itype .eq. 5 ) then 
c
        if( iquad .eq. +1 ) then 
        rval=a**2*b/(z*(z**2+a**2)*sqrt(b**2+a**2+z**2)) 
        endif
c
        if( iquad .eq. 0 ) then 
c
c       1/r singularity at this point
c        
        rval=0
        endif
c
        if( iquad .eq. -1 ) then 
        rval=a**2*b/(z*(z**2+a**2)*sqrt(b**2+a**2+z**2)) 
        endif
c
        endif
c
c> dy := int(int(y*z/r^5,y=0..b/a*x),x=0..a);
c                            2     2 1/2      2   2     4     2   2 1/2
c                 a~ (-a~ (a~  + z~ )    + (b~  a~  + a~  + a~  z~ )   )
c       dy := 1/3 ------------------------------------------------------
c                           2   2     4     2   2 1/2    2     2 1/2
c                     z~ (b~  a~  + a~  + a~  z~ )    (a~  + z~ )
c
        if( itype .eq. 6 ) then 
c
        if( iquad .eq. +1 ) then 
        rval=a*(-sqrt(z**2+a**2)+sqrt(b**2+a**2+z**2))/
     $     (z*sqrt(b**2+a**2+z**2)*sqrt(a**2+z**2))
        endif
c
        if( iquad .eq. 0 ) then 
c
c       1/r singularity at this point
c        
        rval=0
        endif
c
        if( iquad .eq. -1 ) then 
        rval=a*(-sqrt(z**2+a**2)+sqrt(b**2+a**2+z**2))/
     $     (z*sqrt(b**2+a**2+z**2)*sqrt(a**2+z**2))
        endif
c
        endif
c
c> simplify(diff(fz,z));
c                                        2     2
c                               b~ a~ (a~  + z~ )
c     ---------------------------------------------------------------------
c         2                 2     2                 2     2     2     2 1/2
c     (-a~  + 2 I a~ z~ + z~ ) (a~  + 2 I a~ z~ - z~ ) (b~  + a~  + z~ )
c
        if( itype .eq. 7 ) then 
c
        if( iquad .eq. +1 ) then 
        rval=b*a/
     $     (sqrt(b**2+a**2+z**2)*(a**2+z**2))
        endif
c
        if( iquad .eq. 0 ) then 
        rval=b*a/
     $     (sqrt(b**2+a**2+z**2)*(a**2+z**2))
        endif
c
        if( iquad .eq. -1 ) then 
        rval=b*a/
     $     (sqrt(b**2+a**2+z**2)*(a**2+z**2))
        endif
c
        endif
c       
c
c> simplify(diff(b*a/(sqrt(b**2+a**2+z**2)*(a**2+z**2)),z));
c                                     2      2      2
c                           b a z (3 a  + 3 z  + 2 b )
c                        - ----------------------------
c                            2    2    2 3/2   2    2 2
c                          (b  + a  + z )    (a  + z )
c
        if( itype .eq. 20 ) then 
c
        if( iquad .eq. +1 ) then 
        rval=-b*a*z*(3*a**2+3*z**2+2*b**2)/
     $     (sqrt(b**2+a**2+z**2)**3*(a**2+z**2)**2)
        endif
c
        if( iquad .eq. 0 ) then 
        rval=-b*a*z*(3*a**2+3*z**2+2*b**2)/
     $     (sqrt(b**2+a**2+z**2)**3*(a**2+z**2)**2)
        endif
c
        if( iquad .eq. -1 ) then 
        rval=-b*a*z*(3*a**2+3*z**2+2*b**2)/
     $     (sqrt(b**2+a**2+z**2)**3*(a**2+z**2)**2)
        endif
c
        endif
c
c
        if( itype .eq. 37 ) then 
c
c        rval=
c     $     (+1/2.0d0*a*z**2+1/6.0d0*a**3)*
c     $     (log(b+sqrt(a**2+b**2+z**2))-log(sqrt(a**2+z**2)))
c     $     +1/6.0d0*a*b*sqrt(a**2+b**2+z**2)
c        call atan2_safe(a*sqrt(a**2+b**2+z**2),abs(z)*b,val3)
c        rval=rval+z**3/3*(atan2(a,b)-val3)

        if( iquad .eq. 1 ) then
        rval=1/6.0d0*a*b*sqrt(a**2+b**2+z**2)
        call aloga_safe1(a,z,val1)
        call aloga_safe2(a,b,z,val2)
        rval=rval+(1/2.0d0*z**2+1/6.0d0*a**2)*(val2-val1)
        call atan2_safe(a*sqrt(a**2+b**2+z**2),abs(z)*b,val3)
        rval=rval+z**3/3*(atan2(a,b)-val3)
        endif

        if( iquad .eq. 0 ) then
        z0=0.0d0
        rval=1/6.0d0*a*b*sqrt(a**2+b**2)
        call aloga_safe1(a,z0,val1)
        call aloga_safe2(a,b,z0,val2)
        rval=rval+(1/6.0d0*a**2)*(val2-val1)
        endif

        if( iquad .eq. -1 ) then
        rval=1/6.0d0*a*b*sqrt(a**2+b**2+z**2)
        call aloga_safe1(a,z,val1)
        call aloga_safe2(a,b,z,val2)
        rval=rval+(1/2.0d0*z**2+1/6.0d0*a**2)*(val2-val1)
        call atan2_safe(a*sqrt(a**2+b**2+z**2),abs(z)*b,val3)
        rval=rval-z**3/3*(atan2(a,b)-val3)
        endif

        endif
c
        return
        end
c
c
c
c
c
        subroutine triarquad_ab(itype,iquad,a,b,y,z,rval)
        implicit none
        integer itype,iquad,iftaylor
        real *8 a,b,y,z,rval
        real *8 vert1(2),vert2(2)
        real *8 c,hval
c
c       This subroutine integrates constant densities on flat triangles 
c
c
c       The triangle is defined in R^3 and its vertices are
c
c       (a,0,0), (b,0,0), (0,y,z)
c       
c             . y
c            / \
c           /   \
c          /     \
c         /   T   \
c        /_________\
c       a           b
c
c       All integrals are evaluated at the point (0,y,z)
c
c       \int_T  f_{itype} (x,y,z) dT
c
c       itype=1:   \int_T  1/r dT
c       itype=4:   \int_T  z/r^3 dT
c       itype=7:   \int_T  -1/r^3+3*z*z/r^5 dT
c       itype=20:  \int_T  9*z0/r^5-15*z0*z0*z0/r^7 dT
c
c       where r=1/sqrt(x^2+y^2+z^2).
c
c       iquad is the sign of z
c
c       iquad =+1, if z>0
c       iquad = 0, if z=0
c       iquad =-1, if z<0
c
c
c       diff(1/r,z)=-z/r^3
c       diff(1/r,z,z)=-1/r^3+3*z*z/r^5
c       diff(1/r,z,z,z)=9*z0/r^5    -15*z0*z0*z0/r^7
c
c
        rval=0
c
        c=sqrt(y**2+z**2)
c
c
        if( itype .eq. 1 ) then
c
        iftaylor = 0
        if( abs(c) .lt. 1d-8*abs(a) ) iftaylor=1
        if( abs(c) .lt. 1d-8*abs(b) ) iftaylor=1
        if( iftaylor .eq. 1 ) then
        if( a .gt. 0 .and. b .gt. 0 ) then
        rval = +(log(b/a) + c**2*(1/b**2-1/a**2)/4)*y
        if( iquad .ne. 0 ) 
     $     rval = rval -abs(z)*(-.5d0*y*z*(1/a**2-1/b**2) 
     $                 -(atan2(y,b)-atan2(y,a)))        
        endif
        if( a .lt. 0 .and. b .lt. 0 ) then
        rval = -(log(b/a) + c**2*(1/b**2-1/a**2)/4)*y
        if( iquad .ne. 0 ) 
     $     rval = rval -abs(z)*(-.5d0*y*z*(1/a**2-1/b**2) 
     $                 -(atan2(y,b)-atan2(y,a)))
        if( a .gt. 0 .and. b .lt. 0 ) rval = 0
        if( a .lt. 0 .and. b .gt. 0 ) rval = 0
        endif
        else        
cc        rval = log(b+sqrt(c**2+b**2)) - log(a+sqrt(c**2+a**2))
        call triahquad_sing1(a,b,c,rval)
        rval=y*rval
c
        if( iquad .ne. 0 ) then
        call triahquad_sing4(a,b,c,y,z,hval)
        rval = rval - abs(z)*hval
        endif
c
        endif
c
        endif
c
c
        if( itype .eq. 4 ) then
c
        if( iquad .eq. 0 ) then
           rval = 0
        else
c
        iftaylor = 0
        if( abs(c) .lt. 1d-8*abs(a) ) iftaylor=1
        if( abs(c) .lt. 1d-8*abs(b) ) iftaylor=1
        if( iftaylor .eq. 1 ) then
        if( a .gt. 0 .and. b .gt. 0 ) 
     $     rval = -.5d0*y*z*(1/a**2-1/b**2) -(atan2(y,b)-atan2(y,a))
        if( a .lt. 0 .and. b .lt. 0 ) 
     $     rval = -.5d0*y*z*(1/a**2-1/b**2) -(atan2(y,b)-atan2(y,a))
        if( a .gt. 0 .and. b .lt. 0 ) rval = 0
        if( a .lt. 0 .and. b .gt. 0 ) rval = 0
        else        
c        call atan2_safe(y*sqrt(b**2+c**2),abs(z)*b,val3)
c        call atan2_safe(y*sqrt(a**2+c**2),abs(z)*a,val4)
c        rval=rval+(val3-val4)
c        call atan2_safe(y,b,val3)
c        call atan2_safe(y,a,val4)
c        rval=rval-(val3-val4)
        call triahquad_sing4(a,b,c,y,z,rval)
        endif
c
        if( iquad .eq. -1) rval=-rval
c
        endif
c
ccc        write(21,*) iftaylor, rval
c
        endif
c
c
        if( itype .eq. 7 ) then
c
        iftaylor = 0
        if( abs(c) .lt. 1d-8*abs(a) ) iftaylor=1
        if( abs(c) .lt. 1d-8*abs(b) ) iftaylor=1
        if( iftaylor .eq. 1 ) then
        if( a .gt. 0 .and. b .gt. 0 ) rval = +y*(1/a**2-1/b**2)/2
        if( a .lt. 0 .and. b .lt. 0 ) rval = -y*(1/a**2-1/b**2)/2
        if( a .gt. 0 .and. b .lt. 0 ) rval = 0
        if( a .lt. 0 .and. b .gt. 0 ) rval = 0
        else
ccc        rval = y*b/c**2/sqrt(c**2+b**2) - y*a/c**2/sqrt(c**2+a**2)
        call triahquad_sing7(a,b,c,rval)
        rval=y*rval
        endif
c
        endif
c
c
        if( itype .eq. 20 ) then
c
        iftaylor = 0
        if( abs(c) .lt. 1d-8*abs(a) ) iftaylor=1
        if( abs(c) .lt. 1d-8*abs(b) ) iftaylor=1
        if( iftaylor .eq. 1 ) then
        if( a .gt. 0 .and. b .gt. 0 ) 
     $    rval = 0 
ccc     $           +.75d0*y*z*(1/b**4-1/a**4)
ccc     $           +1.25d0*y*z*c**2*(1/b**6-1/a**6)
        if( a .lt. 0 .and. b .lt. 0 ) 
     $    rval = 0
ccc     $           -.75d0*y*z*(1/b**4-1/a**4)
ccc     $           -1.25d0*y*z*c**2*(1/b**6-1/a**6)
        if( a .gt. 0 .and. b .lt. 0 ) rval = 0
        if( a .lt. 0 .and. b .gt. 0 ) rval = 0
        else
c        rval1=-b*y*z*(3*c**2+2*b**2)/
c     $     (sqrt(b**2+c**2)**3*(c**2)**2)
c        rval2=-a*y*z*(3*c**2+2*a**2)/
c     $     (sqrt(a**2+c**2)**3*(c**2)**2)
c        rval=rval1-rval2

        call triahquad_sing20(a,b,c,rval)
        rval= -y*z*rval
c
        endif
c
        endif
c
c
        return
        end
c
c
c
c
c
        subroutine triahquad_sing_aux1(a,b,c,rval)
        implicit none
        integer itype
        real *8 a,b,c,rval
c
c       evaluate  b*sqrt(b**2+c**2) - a*sqrt(a**2+c**2)
c
c
        itype = 1
        if( a .gt. 0 .and. b. gt. 0 ) itype = 2
        if( a .lt. 0 .and. b. lt. 0 ) itype = 2

        if( itype .eq. 1 ) then
        rval=b*sqrt(b**2+c**2)-a*sqrt(a**2+c**2)
        endif

        if( itype .eq. 2 ) then
ccc        rval=b*sqrt(b**2+c**2)-a*sqrt(a**2+c**2)
        rval=(b**2*(b**2+c**2)-a**2*(a**2+c**2))/
     $     (b*sqrt(b**2+c**2)+a*sqrt(a**2+c**2))
        endif

        return
        end
c
c
c
c
        subroutine triahquad_sing1(a,b,c,rval)
        implicit none
        integer itype
        real *8 a,b,c,rval
        real *8 val1,val2
c
c       evaluate
c
c       log(b+sqrt(c**2+b**2)) - log(a+sqrt(c**2+a**2))
c       log((b+sqrt(c**2+b**2))/(a+sqrt(c**2+a**2))) 
c
c
        itype = 1
        if( a .gt. 0 .and. b. gt. 0 ) itype = 2
        if( a .lt. 0 .and. b. lt. 0 ) itype = 2

        if( itype .eq. 1 ) then
        call alog_safe3(b,c,val1)
        call alog_safe3(a,c,val2)
        rval=(val1-val2)
        endif

        if( itype .eq. 2 ) then
        call alog_safe3(b,c,val1)
        call alog_safe3(a,c,val2)
        rval=(val1-val2)
        endif

        return
        end
c
c
c
c
        subroutine triahquad_sing4(a,b,c,y,z,rval)
        implicit none
        integer itype
        real *8 a,b,c,y,z,rval
        real *8 val1,val2,val3,val4
c
c       evaluate
c
        itype = 1
        if( a .gt. 0 .and. b. gt. 0 ) itype = 2
        if( a .lt. 0 .and. b. lt. 0 ) itype = 2

        rval = 0

        if( itype .eq. 1 ) then
        call atan2_safe(y*sqrt(b**2+c**2),abs(z)*b,val1)
        call atan2_safe(y*sqrt(a**2+c**2),abs(z)*a,val2)
        rval=rval+(val1-val2)
        call atan2_safe(y,b,val3)
        call atan2_safe(y,a,val4)
        rval=rval-(val3-val4)
        endif

        if( itype .eq. 2 ) then
        call atan2_safe(y*sqrt(b**2+c**2),abs(z)*b,val1)
        call atan2_safe(y*sqrt(a**2+c**2),abs(z)*a,val2)
        rval=rval+(val1-val2)
        call atan2_safe(y,b,val3)
        call atan2_safe(y,a,val4)
        rval=rval-(val3-val4)
        endif
c
        return
        end
c
c
c
c
        subroutine triahquad_sing7(a,b,c,rval)
        implicit none
        integer itype
        real *8 a,b,c,rval
c
c       evaluate   b/c**2/sqrt(c**2+b**2) - a/c**2/sqrt(c**2+a**2)
c
        itype = 1
        if( a .gt. 0 .and. b. gt. 0 ) itype = 2
        if( a .lt. 0 .and. b. lt. 0 ) itype = 2

        if( itype .eq. 1 ) then
        rval = b/c**2/sqrt(c**2+b**2) - a/c**2/sqrt(c**2+a**2)
        endif

        if( itype .eq. 2 ) then
        rval = (b-a)*(b+a)
     $     /sqrt(c**2+b**2)/sqrt(c**2+a**2)
     $     /(b*sqrt(a**2+c**2)+a*sqrt(b**2+c**2))   
        endif
c
        return
        end
c
c
c
c
        subroutine triahquad_sing20(a,b,c,rval)
        implicit none
        integer itype
        real *8 a,b,c,rval
c
c       evaluate   
c
c       +b*(2*c**2+2*b**2+c**2)/
c          (sqrt(b**2+c**2)**3*(c**2)**2)
c       -a*(2*c**2+2*a**2+c**2)/
c          (sqrt(a**2+c**2)**3*(c**2)**2)
c
        itype = 1
        if( a .gt. 0 .and. b. gt. 0 ) itype = 2
        if( a .lt. 0 .and. b. lt. 0 ) itype = 2

        if( itype .eq. 1 ) then
        rval = 
     $     +b*(2*c**2+2*b**2+c**2)/
     $       (sqrt(b**2+c**2)**3*(c**2)**2)
     $     -a*(2*c**2+2*a**2+c**2)/
     $       (sqrt(a**2+c**2)**3*(c**2)**2)
        endif

        if( itype .eq. 2 ) then
c        rval = 
c     $     +2*b/(sqrt(b**2+c**2)*(c**2)**2)
c     $     -2*a/(sqrt(a**2+c**2)*(c**2)**2)
c     $     +b/(sqrt(b**2+c**2)**3*(c**2))
c     $     -a/(sqrt(a**2+c**2)**3*(c**2))
        rval = 
     $     1/c**2*
     $     (2*(b-a)*(b+a)
     $     /sqrt(c**2+b**2)/sqrt(c**2+a**2)
     $     /(b*sqrt(a**2+c**2)+a*sqrt(b**2+c**2))
     $     +b/(sqrt(b**2+c**2)**3)
     $     -a/(sqrt(a**2+c**2)**3))
        endif
c
        return
        end
c
c
c
c
        subroutine triahquad_sing20h(a,b,c,rval)
        implicit none
        integer itype
        real *8 a,b,c,rval
c
c       evaluate   
c
c       +b**3/(sqrt(b**2+c**2)**3*(c**2))
c       -a**3/(sqrt(a**2+c**2)**3*(c**2))
c
        itype = 1
        if( a .gt. 0 .and. b. gt. 0 ) itype = 2
        if( a .lt. 0 .and. b. lt. 0 ) itype = 2

        if( itype .eq. 1 ) then
        rval = 
     $     +b**3/(sqrt(b**2+c**2)**3)
     $     -a**3/(sqrt(a**2+c**2)**3)
        rval = rval/c**2
        endif

        if( itype .eq. 2 ) then
        rval = 
     $     (+b**3*(sqrt(a**2+c**2))
     $      -a**3*(sqrt(b**2+c**2)))
     $     + 
     $     a**2*b**2*(b-a)*(b+a)/
     $     (b*sqrt(a**2+c**2)+a*sqrt(b**2+c**2))
        rval = rval/(sqrt(b**2+c**2)**3)/(sqrt(a**2+c**2)**3)

        endif
c
        return
        end
c
c
c
c
        subroutine aloga_safe1(a,z,val)
        implicit none
        real *8 a,z,val
c
c       evaluate log(sqrt(a**2+z**2))*a
c
ccc        if( abs(sqrt(a**2+z**2)) .lt. 1d-14 ) then
        if( abs(sqrt(a**2+z**2)) .lt. 1d-200 ) then
        val = 0
        else
        val = log(sqrt(a**2+z**2))*a
        endif
c
        return
        end
c
c
c
c
        subroutine aloga_safe2(a,b,z,val)
        implicit none
        real *8 a,b,z,val
c
c       evaluate log(b+sqrt(a**2+b**2+z**2))*a
c
        call alog_safe3(b,sqrt(a**2+z**2),val)
        val=val*a
c
        return
c
ccc        if( abs(b+sqrt(a**2+b**2+z**2)) .lt. 1d-14 ) then
        if( (b+sqrt(a**2+b**2+z**2)) .lt. 1d-200 ) then
        val = 0
        else
ccc        val = log(b+sqrt(a**2+b**2+z**2))*a
        call alog_safe3(b,sqrt(a**2+z**2),val)
        val=val*a
        endif
c
        return
        end
c
c
c
        subroutine alog_safe3(a,c,val)
        implicit none
        real *8 a,c,val
c
c       evaluate log(a+sqrt(a**2+c**2)) without catastrophic cancellation
c
        if ( a .gt. 0 ) val = log(a+sqrt(c**2+a**2))
        if ( a .le. 0 ) then
        if( -a+sqrt(c**2+a**2) .eq. 0 .or. c**2 .eq. 0 ) then
           val = 0
        else
        val = log((c**2))-log(-a+sqrt(c**2+a**2))
ccc           val = log((c**2)/(-a+sqrt(c**2+a**2)))
        endif
        endif
c
        return
        end
c
c
c
        subroutine atan2_safe(y,x,val)
        implicit none
        real *8 y,x,val
c
c       evaluate atan2(y,x)
c
        if( abs(x) .eq. 0 .and. abs(y) .eq. 0 ) then
        val = 0 
        else
        val = atan2(y,x)
        endif
c
        return
        end
c
c
c
c
