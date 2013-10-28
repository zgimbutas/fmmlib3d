cc Copyright (C) 2009-2012: Leslie Greengard and Zydrunas Gimbutas
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c    Translation operators for multipole expansions.
c    The principal user-callable routines are 
c
cc   l3dmpmpquadu:       multipole -> multipole translation
c
cc   l3dmpmpquadu_add:   multipole -> multipole translation,
c                        *incrementing* second expansion.
c
cc   l3dmpmpzshift:      utility function, 
c                        mp-mp shift along z axis.
c
cc   l3dmplocquadu:       multipole -> local translation
c
cc   l3dmplocquadu_add:   multipole -> local translation,
c                        *incrementing* second expansion.
c
cc   l3dmplocquadu_trunc: multipole -> local translation,
c                        allows different expansion dimensions.
c
cc   l3dmplocquadu_add_trunc: multipole -> local translation,
c                        *incrementing* second expansion and
c                        allowing different expansion dimensions.
c
cc   l3dmploczshiftstab:  utility function, 
c                        mp-loc shift along z axis.
c
cc   l3dmplocquadu2_add_trunc: multipole -> local translation,
c                        *incrementing* second expansion and
c                        allowing different expansion dimensions.
c                        *OPTIMIZED VERSION* for small value of nterms.
c
cc   l3dmplocquadu2_trunc: multipole -> local translation,
c                        allowing different expansion dimensions.
c                        *OPTIMIZED VERSION* for small value of nterms.
c
cc   rotprint: utility printing function
c
cc   l3dmploczshiftstab_fast:  utility function, 
c                        mp-loc shift along z axis (OPTIMIZED).
c
cc   l3dloclocquadu:       local -> local translation
c
cc   l3dloclocquadu_add:   local -> local translation,
c                        *incrementing* second expansion.
c
cc   l3dlocloczshift:    utility function, 
c                        loc-loc shift along z axis.
c
C***********************************************************************
C     Multipole -> Multipole routines
C***********************************************************************
C
C***********************************************************************
      subroutine l3dmpmpquadu(sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,mpolen,nterms2,ier)
C***********************************************************************
C
C     memory management wrapper for 
C     subroutine l3dmpmpquad0 (below).
C
C     Usage:
C
C     Shift center of multipole expansion.
C     This is a reasonably fast "point and shoot" version which
C     first rotates the coordinate system, then shifts
C     along the Z-axis, and then rotates back to the original
C     coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     sc1     :  scaling parameter for mpole expansion
C     x0y0z0  :  center of original multiple expansion
C     mpole   :  coefficients of original multiple expansion
C     nterms  :  order of multipole expansion
C     sc2     :  scaling parameter for shifted expansion
C     xnynzn  :  center of shifted expansion
C     nterms2 :  order of shifted expansion
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     mpolen = coefficients of shifted mpole expansion
C     ier   = error return flag
C             CURRENTLY UNUSED
C
C     Notes: Work arrays carved out of w.
C
C           marray   = work array used to hold various intermediate 
C                      rotated expansions.
C           dc       = work array contain the square roots of 
C                      some binomial coefficients.
C           rd1,rd2  = work arrays used to compute rotation matrices
C                      about Y-axis recursively.
C           ephi     = work array 
C           fr      = work array 
C
C***********************************************************************
      implicit none
      integer  nterms,nterms2,ier,l,m,jnew,knew
      integer  ldc,imarray,lmarray,imarray1,lmarray1,iephi,lephi
      integer  ifr,lused
      real *8 x0y0z0(3),xnynzn(3)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpolen(0:nterms2,-nterms2:nterms2)
      complex *16 imag
c
c     local allocated workspace array
c
      real *8, allocatable :: w(:)
      complex *16, allocatable :: cw(:)
c
      data imag/(0.0d0,1.0d0)/
C
      ldc = max(nterms,nterms2)
      imarray = 1
      lmarray = (ldc+1)*(2*ldc+1) 
      imarray1 = imarray+lmarray
      lmarray1 = (ldc+1)*(2*ldc+1) 
      iephi = imarray1+lmarray1
      lephi = (2*ldc+3) 
      lused = iephi + lephi
      allocate (cw(lused))
      allocate (w(2*ldc+3))
c
      call l3dmpmpquad0(sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,mpolen,nterms2,cw(imarray),cw(imarray1),
     2           ldc,cw(iephi),w,ier)
      return
      end
c
c
C***********************************************************************
      subroutine l3dmpmpquadu_add(sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,mpolen,ldc,nterms2,ier)
C***********************************************************************
C
C     memory management wrapper for 
C     subroutine l3dmpmpquad0 (below).
C
C     Usage:
C
C           Shift center of multipole expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts
C           along the Z-axis, and then rotates back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     sc1     : scaling parameter for mpole expansion
C     x0y0z0  : center of original multiple expansion
C     mpole   : coefficients of original multiple expansion
C     nterms  : order of multipole expansion
C     sc2     : scaling parameter for shifted expansion
C     xnynzn  : center of shifted expansion
C     ldc     : dimension parameter of shifted expansion
C     nterms2 : order of shifted expansion
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     mpolen  : coefficients of shifted mpole expansion
C     ier     : error return flag
C                     CURRENTLY UNUSED
C
C     Notes: Work arrays carved out of w.
C
C           marray   = work array used to hold various intermediate 
c                      rotated expansions.
C           dc       = work array contain the square roots of 
C                      some binomial coefficients.
C           rd1,rd2  = work arrays used to compute rotation matrices
C                      about Y-axis recursively.
C           ephi     = work array 
C           fr      = work array 
C
C***********************************************************************
      implicit none
      integer  nterms,nterms2,ier,l,m,jnew,knew,ldc
      real *8 x0y0z0(3),xnynzn(3)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpolen(0:ldc,-ldc:ldc)
      complex *16 imag
c
c     local allocated workspace array
c
      complex *16, allocatable :: mptemp(:,:)
c
      data imag/(0.0d0,1.0d0)/
C
      allocate( mptemp(0:nterms2,-nterms2:nterms2) )

      call l3dmpmpquadu(sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,mptemp,nterms2,ier)

      do l = 0,min(ldc,nterms2)
         do m=-l,l
            mpolen(l,m) = mpolen(l,m)+mptemp(l,m)
         enddo
      enddo

      return
      end
c
c
c
c
C***********************************************************************
      subroutine l3dmpmpquad0(sc1,x0y0z0,mpole,nterms,sc2,
     1           xnynzn,mpolen,nterms2,marray,marray1,ldc,ephi,
     2           fr,ier)
C***********************************************************************
C
C     Usage:
C
C     Shift multipole expansion.
C     This is a reasonably fast "point and shoot" version which
C     first rotates the coordinate system, then doing the shifting
C     along the Z-axis, and then rotating back to the original
C     coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     sc1     : scaling parameter for mpole expansion
C     x0y0z0  : center of original multiple expansion
C     mpole   : coefficients of original multiple expansion
C     nterms  : order of multipole expansion
C     sc2     : scaling parameter for shifted expansion
C     xnynzn  : center of shifted expansion
C     nterms2 : order of shifted expansion
C     marray  : work array
C     marray1 : work array
C     ldc     : dimension parameter for work arrays
C     ephi    : work array
C     fr      : work array
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     mpolen  : coefficients of shifted expansion
C     ier     : error flag - UNUSED.
C
C     Work Arrays:
C
C           marray = work array used to hold various intermediate 
c                    expansions.
C           dc     = work array contain the square roots of 
C                    some binomial coefficients.
C           ldc      determines dimension of dc
c                    must exceed max(nterms,nterms2).
C           rd     = work arrays used to store rotation matrices
C                    about Y-axis.
C           ephi    = work array 
C
C***********************************************************************
C
      implicit none
      integer  nterms, lw, lused, ier, nq, nquad, nquse,ldc,nterms2
      real *8 x0y0z0(3),xnynzn(3)
      real *8 rshift
      real *8 d,theta,ctheta,phi,sc1,sc2,rvec(3)
      real *8 fr(0:nterms+1)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpolen(0:nterms2,-nterms2:nterms2)
      complex *16 marray1(0:ldc,-ldc:ldc)
      complex *16 marray(0:ldc,-ldc:ldc)
c
      complex *16 ephi(-ldc-1:ldc+1),imag
      integer  l,m,jnew,knew
      data imag/(0.0d0,1.0d0)/
C
      rvec(1) = xnynzn(1) - x0y0z0(1)
      rvec(2) = xnynzn(2) - x0y0z0(2)
      rvec(3) = xnynzn(3) - x0y0z0(3)
      call cart2polarl(rvec,d,theta,phi)
c
      ephi(1) = exp(imag*phi)
      ephi(0)=1.0d0
      ephi(-1)=dconjg(ephi(1))
c
c----- create array of powers e^(i*m*phi).
c
      do l = 1,ldc
         ephi(l+1) = ephi(l)*ephi(1)
         ephi(-1-l) = dconjg(ephi(l+1))
      enddo
c
c----- a rotation of THETA radians about the Yprime axis after PHI
c      radians about the z-axis.
c      The PHI rotation is carried out on the fly by multiplying 
c      mpole and ephi inside the following loop. 
c
      do l=0,nterms
         do m=-l,l
            marray1(l,m)=mpole(l,m)*ephi(m)
         enddo
      enddo
      do l=0,nterms2
         do m=-nterms2,nterms2
            marray(l,m)= 0.0d0
         enddo
      enddo
      call rotviarecur3f90(theta,nterms,nterms,nterms,marray1,
     1        ldc,marray,ldc)
c
c
c----- shift the mpole expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      rshift = d
      call l3dmpmpzshift(sc1,marray,ldc,nterms,sc2,mpolen,
     1           nterms2,nterms2,rshift,fr)
c
c
c     Reverse THETA rotation.
c     I.e. rotation of -THETA radians about Yprime axis.
c
      call rotviarecur3f90(-theta,nterms2,nterms2,nterms2,mpolen,
     1        nterms2,marray,ldc)
c
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
         do m=-l,l
            mpolen(l,m)=ephi(-m)*marray(l,m)
         enddo
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine l3dmpmpzshift(scale,mpole,lmp,nterms,scale2,mpolen,
     1      lmpn,nterms2,zshift,fr)
c***********************************************************************
c
c     This subroutine converts a multipole expansion centered at the 
c     origin to a multipole expansion centered at (0,0,zhift).
c     The expansion is rescaled to that of the shifted expansion.
c
c     INPUT:
c
c     scale    : scale parameter for mpole
c     mpole    : coefficients of original multipole exp.
c     lmp      : leading dim of mpole (may be a work array)
c     nterms   : number of terms in the orig. expansion
c     scale2   : scale parameter for new expansion (mpolen)
c     lmpn     : leading dim of shifted (may be work array)
c     nterms2  : number of terms in output expansion
c     zshift   : shifting distance along z-axis
c                              (always assumed positive)
c     fr       : work array
c
c     OUTPUT:
c
c     mpolen  (complex *16)  : coefficients of shifted exp.
c
c-----------------------------------------------------------------------
      implicit none
      integer l0,nmax,nterms,nterms2,nquad,ier,lmp,lmpn,ldc
      integer l,m,jnew,knew
      real *8 d,zshift,scale,scale2,rat
      real *8 fr(0:*)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 mpolen(0:lmpn,-lmpn:lmpn)
      real *8, allocatable :: dc(:,:)
      real *8, allocatable :: carray(:,:)
C
C----- shift along z-axis
C
      nmax = max(nterms,nterms2)
c
      allocate( dc(0:2*nmax,0:2*nmax) )
      allocate( carray(0:2*nmax,0:2*nmax) )
c
      do l = 0,2*nmax
         carray(l,0) = 1.0d0
         dc(l,0) = 1.0d0
      enddo
      do m=1,2*nmax
         carray(m,m) = 1.0d0
         dc(m,m) = 1.0d0
         do l=m+1,2*nmax
	    carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
	    dc(l,m)=sqrt(carray(l,m))
         enddo
      enddo
c
      d = zshift
      fr(0) = 1.0d0
      d = d*scale
      fr(1) = d
      do l=2,nmax+1
         fr(l) = fr(l-1)*d
      enddo
c
      do jnew = 0,nterms2
         l0 = max(0,jnew-nterms)
         do knew = -jnew,jnew
	    mpolen(jnew,knew) = 0.0d0
	    do l=l0,jnew-iabs(knew)
	       mpolen(jnew,knew)=mpolen(jnew,knew)+mpole(jnew-l,knew)*
     1             fr(l)*dc(jnew-knew,l)*dc(jnew+knew,l)*(-1)**l
            enddo
         enddo
      enddo
      do jnew = 1,nterms2
         do knew = -jnew,jnew
	    mpolen(jnew,knew)=mpolen(jnew,knew)*(scale2/scale)**jnew
         enddo
      enddo
      return
      end
C
C
C***********************************************************************
C     Multipole -> Local routines
C***********************************************************************
C
C***********************************************************************
      subroutine l3dmplocquadu(sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,local,nterms2,ier)
C***********************************************************************
C
C     Memory management wrapper for subroutine l3dmplocquad0 (below).
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts along
C           the Z-axis, and then rotates back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     sc1     : scaling parameter for mpole expansion
C     x0y0z0  : center of original multiple expansion
C     mpole   : coefficients of original multiple expansion
C     nterms  : order of multipole expansion
C     sc2     : scaling parameter for local expansion
C     xnynzn  : center of shifted local expansion
C     nterms2 : order of local expansion
C---------------------------------------------------------------------
C     OUTPUT:
C
C     local = coefficients of shifted local expansion
C     ier   = error return flag
C                   CURRENTLY NOT USED
C
C     Work arrays carved out of w.
C
C           marray = work array used to hold various intermediate 
C                    rotated expansions.
C           dc     = work array contain the square roots of 
C                    some binomial coefficients.
C           rd1,rd2  = work arrays used to compute rotation matrices
C                    about Y-axis recursively.
C           ephi    = work array 
C
C***********************************************************************
C
      implicit none
      integer nterms,nterms2,ier,l,m,jnew,knew,lused
      integer  ldc,imarray,lmarray,imarray1,lmarray1,iephi,lephi
      real *8 x0y0z0(3),xnynzn(3)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 local(0:nterms2,-nterms2:nterms2)
      complex *16 imag
c
c     local allocated workspace array
c
      real *8, allocatable :: w(:)
      complex *16, allocatable :: cw(:)
c
      data imag/(0.0d0,1.0d0)/
C
      ldc = max(nterms,nterms2)
      imarray = 1
      lmarray = (ldc+1)*(2*ldc+1) 
      imarray1 = imarray+lmarray
      lmarray1 = (ldc+1)*(2*ldc+1)
      iephi = imarray1+lmarray1
      lephi = (2*ldc+3)
      lused = iephi+lephi
      allocate (cw(lused))
      allocate (w(2*ldc+3))
c
      call l3dmplocquad0(sc1,x0y0z0,mpole,nterms,sc2,xnynzn,
     1         local,nterms2,cw(imarray),cw(imarray1),ldc,
     2         cw(iephi),w,ier)

      return
      end
c
c
C***********************************************************************
      subroutine l3dmplocquadu_add(sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,local,ldc,nterms2,ier)
C***********************************************************************
C
C     Memory management wrapper for subroutine l3dmplocquad0 (below).
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts along
C           the Z-axis, and then rotates back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     sc1     = scaling parameter for mpole expansion
C     x0y0z0 = center of original multiple expansion
C     mpole  = coefficients of original multiple expansion
C     nterms = order of multipole expansion
C     sc2     = scaling parameter for local expansion
C     xnynzn = center of shifted local expansion
C     ldc    = dimension parameter of local expansion
C     nterms2 = order of local expansion
C---------------------------------------------------------------------
C     OUTPUT:
C
C     local = coefficients of shifted local expansion
C     ier   = error return flag
C                   CURRENTLY NOT USED
C
C     Notes: Work arrays carved out of w.
C
C           marray = work array used to hold various intermediate 
C                    rotated expansions.
C           dc     = work array contain the square roots of 
C                    some binomial coefficients.
C           rd1,rd2  = work arrays used to compute rotation matrices
C                    about Y-axis recursively.
C           ephi    = work array 
C
C***********************************************************************
C
      implicit none
      integer nterms,nterms2,ldc,ier,l,m,jnew,knew
      real *8 x0y0z0(3),xnynzn(3)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 local(0:ldc,-ldc:ldc)
      complex *16 imag
c
c     local allocated workspace array
c
      complex *16, allocatable :: mptemp(:,:)
c
      data imag/(0.0d0,1.0d0)/
C
      allocate( mptemp(0:nterms2,-nterms2:nterms2) )

      call l3dmplocquadu(sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,mptemp,nterms2,ier)

      do l = 0,min(ldc,nterms2)
         do m=-l,l
            local(l,m) = local(l,m)+mptemp(l,m)
         enddo
      enddo

      return
      end
c
c
c
c
C***********************************************************************
      subroutine l3dmplocquad0(sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,local,nterms2,marray,marray1,ldc,ephi,
     2           fr,ier)
C***********************************************************************

C     USAGE:
C
C           Convert multipole expansion to a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then doing the shifting
C           along the Z-axis, and then rotating back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     sc1        scaling parameter for mpole expansion
C     x0y0z0     center of original multiple expansion
C     mpole      coefficients of original multiple expansion
C     nterms     order of multipole expansion
C     sc2        scaling parameter for local expansion
C     xnynzn     center of shifted local expansion
C     nterms2    order of local expansion
C     marray     work array
C     marray1    work array
C     ldc        dimension parameter for work arrays
C     ephi       work array
C     fr         work array
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     local      coefficients of shifted local expansion
C---------------------------------------------------------------------
      implicit none
      integer  nterms,ier,l,m,jnew,knew,nterms2,ldc,mp
      real *8 d,theta,ctheta,phi,sc1,sc2
      real *8 x0y0z0(3),xnynzn(3)
      real *8 rvec(3)
      real *8 rshift
      real *8  fr(0:2*ldc)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 marray1(0:ldc,-ldc:ldc)
      complex *16 local(0:nterms2,-nterms2:nterms2)
      complex *16 marray(0:ldc,-ldc:ldc)
      complex *16 ephi(-ldc-1:ldc+1),imag
      data imag/(0.0d0,1.0d0)/
C
      rvec(1) = xnynzn(1) - x0y0z0(1)
      rvec(2) = xnynzn(2) - x0y0z0(2)
      rvec(3) = xnynzn(3) - x0y0z0(3)
      call cart2polarl(rvec,d,theta,phi)
c
      ephi(1) = exp(imag*phi)
      ephi(0)=1.0d0
      ephi(-1)=dconjg(ephi(1))
c
c     create array of powers e^(i*m*phi).
c
      do l = 1,ldc
         ephi(l+1) = ephi(l)*ephi(1)
         ephi(-1-l) = dconjg(ephi(l+1))
      enddo
c
c     a rotation of THETA radians about the Yprime axis after PHI
c     radians about the z-axis.
c     The PHI rotation is carried out on the fly by multiplying 
c     mpole and ephi inside the following loop. 
c
      do l=0,nterms
         do mp=-l,l
            marray1(l,mp)  = mpole(l,mp)*ephi(mp)
         enddo
      enddo
      do l=0,nterms2
         do mp=-nterms2,nterms2
            marray(l,mp)  = 0.0d0
         enddo
      enddo
c
ccc      call rotviarecur3f90(theta,nterms,nterms,nterms,marray1,
      call rotviarecur3f90(theta,nterms,nterms,nterms2,marray1,
     1     ldc,marray,ldc)
c
c
c----- shift the local expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      rshift = d
      call l3dmploczshiftstab(marray,sc1,ldc,nterms,local,
     1      sc2,nterms2,nterms2,rshift,fr)

c
c     reverse THETA rotation. 
c     I.e. rotation of -THETA radians about the Yprime axis.
c
      call rotviarecur3f90(-theta,nterms2,nterms2,nterms2,local,
     1     nterms2,marray,ldc)
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
         do m=-l,l
            local(l,m)=ephi(-m)*marray(l,m)
         enddo
      enddo
      return
      end
c
c
c
C***********************************************************************
      subroutine l3dmplocquadu_trunc(sc1,x0y0z0,mpole,nterms,nterms1,
     1           sc2,xnynzn,local,nterms2,ier)
C***********************************************************************
C
C     Memory management wrapper for subroutine l3dmplocquad0 (below).
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts along
C           the Z-axis, and then rotates back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     sc1      scaling parameter for mpole expansion
C     x0y0z0   center of original multiple expansion
C     mpole    coefficients of original multiple expansion
C     nterms   dimension of original multipole expansion
C     nterms1  truncated order of original multipole expansion
C     sc2      scaling parameter for local expansion
C     xnynzn   center of shifted local expansion
C     nterms2  order of local expansion
C---------------------------------------------------------------------
C     OUTPUT:
C
C     local    coefficients of shifted local expansion
C     ier      error return flag
C              CURRENTLY NOT USED
C
C     Work arrays carved out of w.
C
C           marray = work array used to hold various intermediate 
C                    rotated expansions.
C           dc     = work array contain the square roots of 
C                    some binomial coefficients.
C           rd1,rd2  = work arrays used to compute rotation matrices
C                    about Y-axis recursively.
C           ephi    = work array 
C
C***********************************************************************
C
      implicit none
      integer nterms,nterms1,nterms2,ier,l,m,jnew,knew
      integer  ldc,imarray,lmarray,imarray1,lmarray1,iephi,lephi
      integer  lused
      real *8 x0y0z0(3),xnynzn(3)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 local(0:nterms2,-nterms2:nterms2)
      complex *16 wavek,imag
c
c     local allocated workspace array
c
      real *8, allocatable :: w(:)
      complex *16, allocatable :: cw(:)
c
c
      data imag/(0.0d0,1.0d0)/
C
      ldc = max(nterms,nterms2)
      ldc = max(ldc,nterms1) 
      ldc = ldc+2
      imarray = 1
      lmarray = (ldc+1)*(2*ldc+1)
      imarray1 = imarray+lmarray
      lmarray1 = (ldc+1)*(2*ldc+1)
      iephi = imarray1+lmarray1
      lephi = (2*ldc+3) 
      lused = iephi+ lephi
      allocate (cw(lused))
      allocate (w(2*ldc+3))
c
      call l3dmplocquad_trunc0
     $   (sc1,x0y0z0,mpole,nterms,nterms1,sc2,xnynzn,
     1         local,nterms2,cw(imarray),cw(imarray1),ldc,
     2         cw(iephi),w,ier)
      return
      end
c
c
C***********************************************************************
      subroutine l3dmplocquadu_add_trunc
     $     (sc1,x0y0z0,mpole,nterms,nterms1,
     1           sc2,xnynzn,local,ldc,nterms2,ier)
C***********************************************************************
C
C     Memory management wrapper for subroutine l3dmplocquad0 (below).
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts along
C           the Z-axis, and then rotates back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     sc1       scaling parameter for mpole expansion
C     x0y0z0    center of original multiple expansion
C     mpole     coefficients of original multiple expansion
C     nterms    dimension of multipole expansion
C     nterms1   order of truncated multipole expansion
C     sc2       scaling parameter for local expansion
C     xnynzn    center of shifted local expansion
C     ldc       dimension parameter of local expansion
C     nterms2   order of local expansion
C---------------------------------------------------------------------
C     OUTPUT:
C
C     local = coefficients of shifted local expansion
C     ier   = error return flag
C             CURRENTLY NOT USED
C
C     Work arrays carved out of w.
C
C           marray = work array used to hold various intermediate 
C                    rotated expansions.
C           dc     = work array contain the square roots of 
C                    some binomial coefficients.
C           rd1,rd2  = work arrays used to compute rotation matrices
C                    about Y-axis recursively.
C           ephi    = work array 
C
C***********************************************************************
C
      implicit none
      integer nterms,nterms1,nterms2,ldc,ier,l,m,jnew,knew
      real *8 x0y0z0(3),xnynzn(3)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 local(0:ldc,-ldc:ldc)
      complex *16 wavek,imag
c
c     local allocated workspace array
c
      complex *16, allocatable :: mptemp(:,:)
c
      data imag/(0.0d0,1.0d0)/
C
      allocate( mptemp(0:nterms2,-nterms2:nterms2) )

      call l3dmplocquadu_trunc(sc1,x0y0z0,mpole,nterms,nterms1,
     1           sc2,xnynzn,mptemp,nterms2,ier)

      do l = 0,min(ldc,nterms2)
         do m=-l,l
            local(l,m) = local(l,m)+mptemp(l,m)
         enddo
      enddo

      return
      end
c
c
c
c
C***********************************************************************
      subroutine l3dmplocquad_trunc0(sc1,x0y0z0,mpole,nterms,nterms1,
     1           sc2,xnynzn,local,nterms2,marray,marray1,ldc,ephi,
     2           fr,ier)

C***********************************************************************

C     USAGE:
C
C     Convert multipole expansion to a local expansion.
C     This is a reasonably fast "point and shoot" version which
C     first rotates the coordinate system, then doing the shifting
C     along the Z-axis, and then rotating back to the original
C     coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     sc1      scaling parameter for mpole expansion
C     x0y0z0   center of original multiple expansion
C     mpole    coefficients of original multiple expansion
C     nterms   dimension of multipole expansion
C     nterms1  order of truncated multipole expansion
C     sc2      scaling parameter for local expansion
C     xnynzn   center of shifted local expansion
C     nterms2  order of local expansion
C     marray   work array
C     marray1  work array
C     ldc      dimension parameter for work arrays
C     ephi     work array
C     fr       work array
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     local = coefficients of shifted local expansion
C
C     Work Arrays:
C
C           marray = work array used to hold various intermediate 
c                    expansions.
C           ldc      must exceed max(nterms,nterms2).
C           rd1,rd2  work arrays used to store rotation matrices
C                    about Y-axis.
C           ephi    = work array 
C
C           LOTS MORE
C
C
C---------------------------------------------------------------------
      implicit none
      integer  nterms,nterms1,nterms2,ier,l,m,jnew,knew,ldc,mp
      real *8 d,theta,ctheta,phi,sc1,sc2
      real *8 x0y0z0(3),xnynzn(3)
      real *8 rvec(3)
      real *8 rshift
      real *8  fr(0:2*ldc)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 local(0:nterms2,-nterms2:nterms2)
      complex *16 marray(0:ldc,-ldc:ldc)
      complex *16 marray1(0:nterms1,-nterms1:nterms1)
      complex *16 ephi(-ldc-1:ldc+1),imag
      data imag/(0.0d0,1.0d0)/
C
      rvec(1) = xnynzn(1) - x0y0z0(1)
      rvec(2) = xnynzn(2) - x0y0z0(2)
      rvec(3) = xnynzn(3) - x0y0z0(3)
      call cart2polarl(rvec,d,theta,phi)
      ephi(1) = exp(imag*phi)
c
      ephi(0)=1.0d0
      ephi(-1)=dconjg(ephi(1))
c
c     create array of powers e^(i*m*phi).
c
      do l = 1,ldc
         ephi(l+1) = ephi(l)*ephi(1)
         ephi(-1-l) = dconjg(ephi(l+1))
      enddo
c
c     a rotation of THETA radians about the Yprime axis after PHI
c     radians about the z-axis.
      do l=0,nterms1
         do mp=-l,l
            marray1(l,mp)  = mpole(l,mp)*ephi(mp)
         enddo
      enddo
      do l=0,nterms2
         do mp=-nterms2,nterms2
            marray(l,mp)  = 0.0d0
         enddo
      enddo
c
ccc      call rotviarecur3f90(theta,nterms1,nterms1,nterms1,marray1,
      call rotviarecur3f90(theta,nterms1,nterms1,nterms2,marray1,
     1     nterms1,marray,ldc)
c
c----- shift the local expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      rshift = d
      call l3dmploczshiftstab(marray,sc1,ldc,nterms1,local,
     1      sc2,nterms2,nterms2,rshift,fr)

c
c     reverse THETA rotation. 
c     I.e. rotation of -THETA radians about the Yprime axis.
c
      call rotviarecur3f90(-theta,nterms2,nterms2,nterms2,local,
     1     nterms2,marray,ldc)
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
         do m=-l,l
            local(l,m)=ephi(-m)*marray(l,m)
         enddo
      enddo
      return
      end
c
c
c***********************************************************************
      subroutine l3dmploczshiftstab(mpole,scale,lmp,nterms,local,
     2      scale2,lmpn,nterms2,zshift,fr)
c***********************************************************************
c
c     This subroutine converts a multipole expansion centered at the 
c     origin to a local expansion centered at (0,0,zhift).
c     The expansion is rescaled to that of the local expansion.
c
C---------------------------------------------------------------------
c     INPUT:
c
c     mpole    : coefficients of original multipole exp.
c     scale    : scale parameter for mpole
c     lmp      : leading dim of mpole (may be a work array)
c     nterms   : number of terms in original expansion
c
c     scale2   : scale parameter for local
c     lmpn     : leading dim of local (may be a work array)
c     nterms2  : number of terms in output local exp.
c     zshift   : shifting distance along z-axis
c                             (always assumed positive)
c     fr       : work array (0: max(nterms,nterms2))
c
C---------------------------------------------------------------------
c     OUTPUT:
c
c     local    : coefficients of shifted local exp.
C---------------------------------------------------------------------
      implicit none
      integer nterms,nterms2,nquad,ier
      integer l,lw,m,jnew,knew,kk,ll,msign,nmax,lmpn,lmp
      real *8 scale,scale2
      real *8 zshift,d
      real *8 fr(0:*)
      complex *16 mpole(0:lmp,-lmp:lmp),zk
      complex *16 local(0:lmpn,-lmpn:lmpn)
      real *8, allocatable :: dc(:,:)
      real *8, allocatable :: carray(:,:)
C
C----- shift along z-axis by evaluating field on target sphere and
C     projecting onto spherical harmonics and scaling by j_n(kR).
C
      nmax = max(nterms,nterms2)
c
      allocate( dc(0:2*nmax,0:2*nmax) )
      allocate( carray(0:2*nmax,0:2*nmax) )
c
      do l = 0,2*nmax
         carray(l,0) = 1.0d0
         dc(l,0) = 1.0d0
      enddo
      do m=1,2*nmax
         carray(m,m) = 1.0d0
         dc(m,m) = 1.0d0
         do l=m+1,2*nmax
	    carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
	    dc(l,m)=sqrt(carray(l,m))
         enddo
      enddo
c
      d = 1.0d0/zshift
      fr(0) = d
      d = d/scale
      fr(1) = fr(0)*d
      do l=2,2*nmax
        fr(l) = fr(l-1)*d
      enddo
c
      do jnew = 0,nterms2
         do knew = -jnew,jnew
	    local(jnew,knew) = 0.0d0
	    kk  =iabs(knew)
	    msign = (-1)**(jnew+knew)
	    do l=kk,nterms
	       ll=l+jnew
	       local(jnew,knew)=local(jnew,knew)+mpole(l,knew)*
     1              fr(ll)*dc(ll,l-knew)*dc(ll,l+knew)*msign
            enddo
         enddo
      enddo
      do jnew = 0,nterms2
         do knew = -jnew,jnew
	    local(jnew,knew)=local(jnew,knew)*(scale/scale2)**jnew
         enddo
      enddo
      return
      end
C
c
c
c
c
c
c
c
c
C***********************************************************************
      subroutine l3dmplocquadu2_add_trunc
     $     (sc1,x0y0z0,mpole,nterms,nterms1,
     1           sc2,xnynzn,local,ldc,nterms2,ier,
     $     rotmatf,rotmatb,ldm)
C***********************************************************************
C
C     Memory management wrapper for subroutine l3dmplocquad0 (below).
C
C     Usage:
C
C     Converts multipole expansion to a local expansion.
C     This is a reasonably fast "point and shoot" version which
C     first rotates the coordinate system, then shifts along
C     the Z-axis, and then rotates back to the original
C     coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     sc1       scaling parameter for mpole expansion
C     x0y0z0    center of original multiple expansion
C     mpole     coefficients of original multiple expansion
C     nterms    dimension of multipole expansion
C     nterms1   order of truncated multipole expansion
C     sc2       scaling parameter for local expansion
C     xnynzn    center of shifted local expansion
C     ldc       dimension of local expansion
C     nterms2   order of local expansion
C     rotmatf   precomputed array
C     rotmatb   precomputed array
C     ldm       dimension for rotmatf,rotmatb
C---------------------------------------------------------------------
C     OUTPUT:
C
C     local = coefficients of shifted local expansion
C     ier   = error return flag
C             CURRENTLY NOT USED
C
C     Work arrays carved out of w.
C
C           marray = work array used to hold various intermediate 
C                    rotated expansions.
C           dc     = work array contain the square roots of 
C                    some binomial coefficients.
C           rd1,rd2  = work arrays used to compute rotation matrices
C                    about Y-axis recursively.
C           ephi    = work array 
C
C***********************************************************************
C
      implicit none
      integer nterms,nterms1,nterms2,ier,l,m,jnew,knew,ldc,ldm
      real *8 x0y0z0(3),xnynzn(3)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 local(0:ldc,-ldc:ldc)
      complex *16 wavek,imag
      real *8 rotmatf(0:ldm,0:ldm,-ldm:ldm)
      real *8 rotmatb(0:ldm,0:ldm,-ldm:ldm)
c
c     local allocated workspace array
c
      complex *16, allocatable :: mptemp(:,:)
c
      data imag/(0.0d0,1.0d0)/
C
      allocate( mptemp(0:nterms2,-nterms2:nterms2) )

      call l3dmplocquadu2_trunc(sc1,x0y0z0,mpole,nterms,nterms1,
     1           sc2,xnynzn,mptemp,nterms2,ier,rotmatf,rotmatb,ldm)

      do l = 0,min(ldc,nterms2)
         do m=-l,l
            local(l,m) = local(l,m)+mptemp(l,m)
         enddo
      enddo

      return
      end
c
c
c
c
C***********************************************************************
      subroutine l3dmplocquadu2_trunc(sc1,x0y0z0,mpole,nterms,nterms1,
     1           sc2,xnynzn,local,nterms2,ier,rotmatf,rotmatb,ldm)
C***********************************************************************
C
C     Memory management wrapper for subroutine l3dmplocquad0 (below).
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts along
C           the Z-axis, and then rotates back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     sc1      scaling parameter for mpole expansion
C     x0y0z0   center of original multiple expansion
C     mpole    coefficients of original multiple expansion
C     nterms   dimension of multipole expansion
C     nterms1  order of truncated multipole expansion
C     sc2      scaling parameter for local expansion
C     xnynzn   center of shifted local expansion
C     nterms2  order of local expansion
C     rotmatf  precomputed array
C     rotmatb  precomputed array
C     ldm      dimension for rotmatf,rotmatb
C---------------------------------------------------------------------
C     OUTPUT:
C
C     local = coefficients of shifted local expansion
C     ier   = error return flag
C                   CURRENTLY NOT USED
C
C     Work arrays carved out of w.
C
C           marray = work array used to hold various intermediate 
C                    rotated expansions.
C           dc     = work array contain the square roots of 
C                    some binomial coefficients.
C           rd1,rd2  = work arrays used to compute rotation matrices
C                    about Y-axis recursively.
C           ephi    = work array 
C
C***********************************************************************
C
      implicit none
      integer nterms,nterms1,nterms2,ier,l,m,jnew,knew,ldm
      integer ldc,imarray,lmarray,imarray1,lmarray1,iephi,lephi
      integer ifr,lused
      real *8 x0y0z0(3),xnynzn(3)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 local(0:nterms2,-nterms2:nterms2)
      complex *16 wavek,imag
      real *8 rotmatf(0:ldm,0:ldm,-ldm:ldm)
      real *8 rotmatb(0:ldm,0:ldm,-ldm:ldm)
c
c     local allocated workspace array
c
      real *8, allocatable :: w(:)
c
c
      data imag/(0.0d0,1.0d0)/
C
      ldc = max(nterms,nterms2)
      ldc = max(ldc,nterms1) 
      ldc = ldc+2
      imarray = 1
      lmarray = 2*(ldc+1)*(2*ldc+1) + 3 
      imarray1 = imarray+lmarray
      lmarray1 = 2*(ldc+1)*(2*ldc+1) + 3 
      iephi = imarray1+lmarray1
      lephi = 2*(2*ldc+3) + 3 
      ifr = iephi+lephi
      lused = ifr+ 2*(2*ldc+3) + 3
      allocate (w(lused))
c
      call l3dmplocquad2_trunc0
     $   (sc1,x0y0z0,mpole,nterms,nterms1,sc2,xnynzn,
     1         local,nterms2,w(imarray),w(imarray1),ldc,
     2         w(iephi),w(ifr),ier,rotmatf,rotmatb,ldm)
      return
      end
c
c
c
C***********************************************************************
      subroutine l3dmplocquad2_trunc0(sc1,x0y0z0,mpole,nterms,nterms1,
     1           sc2,xnynzn,local,nterms2,marray,marray1,ldc,ephi,
     2           fr,ier,rotmatf,rotmatb,ldm)

C***********************************************************************

C     USAGE:
C
C     Convert multipole expansion to a local expansion.
C     This is a reasonably fast "point and shoot" version which
C     first rotates the coordinate system, then doing the shifting
C     along the Z-axis, and then rotating back to the original
C     coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     sc1       scaling parameter for mpole expansion
C     x0y0z0    center of original multiple expansion
C     mpole     coefficients of original multiple expansion
C     nterms    dimension of multipole expansion
C     nterms1   order of truncated multipole expansion
C     sc2       scaling parameter for local expansion
C     xnynzn    center of shifted local expansion
C     nterms2   order of local expansion
C     marray    work array
C     marray1   work array
C     ldc       dimension parameter for work arrays
C     ephi      work array
C     fr        work array
C     rotmatf   precomputed array
C     rotmatb   precomputed array
C     ldm       dimension for rotmatf,rotmatb
C---------------------------------------------------------------------
C     OUTPUT:
C
C     local    coefficients of shifted local expansion
C
C     Work Arrays:
C
C           marray = work array used to hold various intermediate 
c                    expansions.
C           ldc      must exceed max(nterms,nterms2).
C           rd1,rd2  work arrays used to store rotation matrices
C                    about Y-axis.
C           ephi    = work array 
C
C           LOTS MORE
C
C
C---------------------------------------------------------------------
      implicit none
      integer  nterms,nterms1,nterms2,ier,l,m,jnew,knew,ldm,ldc,mp
      real *8 d,theta,ctheta,phi,sc1,sc2
      real *8 x0y0z0(3),xnynzn(3)
      real *8 rvec(3)
      real *8 rshift
      real *8  fr(0:2*ldc)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 local(0:nterms2,-nterms2:nterms2)
      complex *16 marray(0:ldc,-ldc:ldc)
      complex *16 marray1(0:nterms1,-nterms1:nterms1)
      complex *16 ephi(-ldc-1:ldc+1),imag
      real *8 rotmatf(0:ldm,0:ldm,-ldm:ldm)
      real *8 rotmatb(0:ldm,0:ldm,-ldm:ldm)
      data imag/(0.0d0,1.0d0)/
C
      rvec(1) = xnynzn(1) - x0y0z0(1)
      rvec(2) = xnynzn(2) - x0y0z0(2)
      rvec(3) = xnynzn(3) - x0y0z0(3)
      call cart2polarl(rvec,d,theta,phi)
      ephi(1) = exp(imag*phi)
c
      ephi(0)=1.0d0
      ephi(-1)=dconjg(ephi(1))
c
c     create array of powers e^(i*m*phi).
c
      do l = 1,ldc
         ephi(l+1) = ephi(l)*ephi(1)
         ephi(-1-l) = dconjg(ephi(l+1))
      enddo
c
c     a rotation of THETA radians about the Yprime axis after PHI
c     radians about the z-axis.
      do l=0,nterms1
         do mp=-l,l
            marray1(l,mp)  = mpole(l,mp)*ephi(mp)
         enddo
      enddo
      do l=0,nterms2
         do mp=-nterms2,nterms2
            marray(l,mp)  = 0.0d0
       enddo
      enddo
c
ccc      call rotviarecur3f90(theta,nterms1,nterms1,nterms1,marray1,
ccc     1     nterms1,marray,ldc)
ccc        call rotviarecur3p_apply(theta,nterms1,nterms1,nterms1,marray1,
        call rotviarecur3p_apply(theta,nterms1,nterms1,nterms2,marray1,
     1     nterms1,marray,ldc,rotmatf,ldm)

c
c----- shift the local expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      rshift = d
      call l3dmploczshiftstab_fast(marray,sc1,ldc,nterms1,local,
     1      sc2,nterms2,nterms2,rshift,fr)

c
c     reverse THETA rotation. 
c     I.e. rotation of -THETA radians about the Yprime axis.
c
ccc      call rotviarecur3f90(-theta,nterms2,nterms2,nterms2,local,
ccc     1     nterms2,marray,ldc)
      call rotviarecur3p_apply(-theta,nterms2,nterms2,nterms2,local,
     1     nterms2,marray,ldc,rotmatb,ldm)
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
         do m=-l,l
            local(l,m)=ephi(-m)*marray(l,m)
         enddo
      enddo
      return
      end
c
c

C*****************************************************************
        subroutine rotprint(nterms,rotmat,ldc)
C*****************************************************************
c
c       Purpose: print the rotation matrix
c
c       nterms: order of multipole expansion
c       rotmat:  real *8 (0:ldc,0:ldc,-ldc:ldc): rotation matrix 
c       ldc   : leading dim for rotation matrix (must exceed nterms)
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       marray   coefficients of rotated expansion.
c
C---------------------------------------------------------------------
	implicit none
	integer ij,m,mp,nterms,ldc
	real *8 rotmat(0:ldc,0:ldc,-ldc:ldc)
        real *8, allocatable :: rd1(:,:)
c
        allocate( rd1(0:ldc,-ldc:ldc) )
c
c       ... print rotation matrix
c
        do ij=0,nterms
c
         do m=-ij,ij
            do mp=0,ij
               rd1(mp,m) = rotmat(ij,mp,m)
            enddo
         enddo

c         if( ij.eq.nterms ) then
         call prinf('rd1=*',ij,1)
         do m=-ij,ij
         call prin2(' *',rd1(0,m),ij+1)
         enddo
c         endif
c
        enddo
c
        return
        end
c
c
c                
c***********************************************************************
      subroutine l3dmploczshiftstab_fast(mpole,scale,lmp,nterms,local,
     2      scale2,lmpn,nterms2,zshift,fr)
c***********************************************************************
c
c     This subroutine converts a multipole expansion centered at the 
c     origin to a local expansion centered at (0,0,zhift).
c     The expansion is rescaled to that of the local expansion.
c
C---------------------------------------------------------------------
c     INPUT:
c
c     mpole    : coefficients of original multipole exp.
c     scale    : scale parameter for mpole
c     lmp      : leading dim of mpole (may be a work array)
c     nterms   : number of terms in original expansion
c
c     scale2   : scale parameter for local
c     lmpn     : leading dim of local (may be a work array)
c     nterms2  : number of terms in output local exp.
c     zshift   : shifting distance along z-axis
c                             (always assumed positive)
c     fr       : work array (0: max(nterms,nterms2))
c
C---------------------------------------------------------------------
c     OUTPUT:
c
c     local    : coefficients of shifted local exp.
c
C---------------------------------------------------------------------
      implicit none
      integer nterms,nterms2,nquad,ier,lmp,lmpn
      integer l,lw,m,jnew,knew,kk,ll,msign,nmax
      real *8 zshift,scale,scale2
      real *8 d
      real *8 fr(0:*)
      complex *16 mpole(0:lmp,-lmp:lmp),zk
      complex *16 local(0:lmpn,-lmpn:lmpn)
      real *8, allocatable :: dc(:,:)
      real *8, allocatable :: carray(:,:)
C
C----- shift along z-axis by evaluating field on target sphere and
C     projecting onto spherical harmonics and scaling by j_n(kR).
C
      nmax = max(nterms,nterms2)
c
      allocate( dc(0:2*nmax,0:2*nmax) )
      allocate( carray(0:2*nmax,0:2*nmax) )
c
      do l = 0,2*nmax
         carray(l,0) = 1.0d0
         dc(l,0) = 1.0d0
      enddo
      do m=1,2*nmax
         carray(m,m) = 1.0d0
         dc(m,m) = 1.0d0
         do l=m+1,2*nmax
	    carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
	    dc(l,m)=sqrt(carray(l,m))
         enddo
      enddo
c
      d = 1.0d0/zshift
      fr(0) = d
      d = d/scale
      fr(1) = fr(0)*d
      do l=2,2*nmax
        fr(l) = fr(l-1)*d
      enddo
c
      d = 1
      do jnew = 0,nterms2
         do knew = -jnew,jnew
	    local(jnew,knew) = 0.0d0
	    kk  =abs(knew)
	    if( mod(abs(jnew+knew),2) .eq. 1 ) then
            msign = -1
            else
            msign = +1
            endif
	    do l=kk,nterms
	       ll=l+jnew
               if( msign .eq. 1 ) then
	       local(jnew,knew)=local(jnew,knew)+mpole(l,knew)*
     1              (fr(ll)*dc(ll,l-knew)*dc(ll,l+knew))
               else
	       local(jnew,knew)=local(jnew,knew)-mpole(l,knew)*
     1              (fr(ll)*dc(ll,l-knew)*dc(ll,l+knew))
               endif
            enddo
         local(jnew,knew)=local(jnew,knew)*d
         enddo
         d = d*scale/scale2
      enddo

      return
      end
C
C
C
C***********************************************************************
C     Local -> Local routines
C***********************************************************************
C
C
C***********************************************************************
      subroutine l3dloclocquadu(sc1,x0y0z0,locold,nterms,
     1           sc2,xnynzn,local,nterms2,ier)
C***********************************************************************
C
C     memory management wrapper for 
C     subroutine l3dloclocquad0 (below).
C
C     Usage:
C
C     Shift center of a local expansion.
C     This is a reasonably fast "point and shoot" version which
C     first rotates the coordinate system, then shifts along
C     the Z-axis, and then rotates back to the original
C     coordinates.
C
C     INPUT:
C
C     sc1       scaling parameter for locold expansion
C     x0y0z0    center of original expansion
C     locold    coefficients of original expansion
C     nterms    order of original expansion
C     sc2       scaling parameter for local expansion
C     xnynzn    center of shifted expansion
C     nterms2   order of shifted expansion
C
C     OUTPUT:
C
C     local = coefficients of shifted expansion
C     ier   = error return flag
C             CURRENTLY NOT USED
C
C***********************************************************************
      implicit none
      integer nterms,nterms2,ier,l,m,jnew,knew
      integer  ldc,imarray,lmarray,imarray1,lmarray1,iephi,lephi
      integer  lused
      real *8 x0y0z0(3),xnynzn(3)
      real *8 sc1,sc2,d,theta,phi,ctheta
      complex *16 locold(0:nterms,-nterms:nterms)
      complex *16 local(0:nterms2,-nterms2:nterms2)
      complex *16 imag
c
c     local allocated workspace array
c
      real *8, allocatable :: w(:)
      complex *16, allocatable :: cw(:)
c
      data imag/(0.0d0,1.0d0)/
C
      ldc = max(nterms,nterms2)
      imarray = 1
      lmarray = (ldc+1)*(2*ldc+1) 
      imarray1 = imarray+lmarray
      lmarray1 = (ldc+1)*(2*ldc+1) 
      iephi = imarray1+lmarray1
      lephi = (2*ldc+3)
      lused = iephi+ lephi
      allocate (cw(lused))
      allocate (w(2*ldc+3))
c
      call l3dloclocquad0(sc1,x0y0z0,locold,nterms,sc2,xnynzn,
     1           local,nterms2,cw(imarray),cw(imarray1),ldc,
     2           cw(iephi),w,ier)
      return
      end
c
c
C***********************************************************************
      subroutine l3dloclocquadu_add(sc1,x0y0z0,locold,nterms,
     1           sc2,xnynzn,local,ldc,nterms2,ier)
C***********************************************************************
C
C     memory management wrapper for 
C     subroutine l3dloclocquad0 (below).
C
C     Usage:
C
C     Shift center of a local expansion.
C     This is a reasonably fast "point and shoot" version which
C     first rotates the coordinate system, then shifts along
C     the Z-axis, and then rotates back to the original
C     coordinates.
C
C     INPUT:
C
C     sc1       scaling parameter for locold expansion
C     x0y0z0    center of original expansion
C     locold    coefficients of original expansion
C     nterms    order of original expansion
C     sc2       scaling parameter for local expansion
C     xnynzn    center of shifted expansion
C     ldc       dimension parameter for local expansion
C     nterms2   order of shifted expansion
C
C     OUTPUT:
C
C     local = coefficients of shifted expansion
C     ier   = error return flag
C                   CURRENTLY NOT USED
C***********************************************************************
      implicit none
      integer nterms,nterms2,ldc,ier,l,m,jnew,knew
      real *8 x0y0z0(3),xnynzn(3)
      real *8 sc1,sc2,d,theta,phi,ctheta
      complex *16 locold(0:nterms,-nterms:nterms)
      complex *16 local(0:ldc,-ldc:ldc)
      complex *16 imag
c
c     local allocated workspace array
c
      complex *16, allocatable :: mptemp(:,:)
c
      data imag/(0.0d0,1.0d0)/
C
      allocate( mptemp(0:nterms2,-nterms2:nterms2) )

      call l3dloclocquadu(sc1,x0y0z0,locold,nterms,
     1           sc2,xnynzn,mptemp,nterms2,ier)

      do l = 0,min(ldc,nterms2)
         do m=-l,l
            local(l,m) = local(l,m)+mptemp(l,m)
         enddo
      enddo

      return
      end
c
c
c
c
C***********************************************************************
      subroutine l3dloclocquad0(sc1,x0y0z0,locold,nterms,sc2,
     1           xnynzn,local,nterms2,marray,marray1,
     2           ldc,ephi,fr,ier) 
C***********************************************************************
C
C     Usage:
C
C           Shifts center of a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then doing the shifting
C           along the Z-axis, and then rotating back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     sc1     : scaling parameter for locold expansion
C     x0y0z0  : center of original multiple expansion
C     locold  : coefficients of original multiple expansion
C     nterms  : order of original local expansion
C     sc2     : scaling parameter for local expansion
C     xnynzn  : center of shifted local expansion
C     nterms2 : order of new local expansion
c     marray  : work array
c     marray1 : work array
c     ldc     : dimension parameter for work arrays
c     ephi    : work array 
c     fr      : work array 
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     local   : coefficients of shifted local expansion
c     ier     : error return code 
c               CURRENTLY NOT USED
C
C***********************************************************************
C
      implicit none
      integer nterms,ier,l,m,jnew,knew,nterms2,ldc,mp
      real *8 x0y0z0(3),xnynzn(3),rvec(3)
      real *8 d,theta,ctheta,phi,sc1,sc2
      real *8 fr(0:ldc+1)
      real *8 rshift
      complex *16 locold(0:nterms,-nterms:nterms)
      complex *16 marray1(0:ldc,-ldc:ldc)
      complex *16 local(0:nterms2,-nterms2:nterms2)
      complex *16 marray(0:ldc,-ldc:ldc)
      complex *16 imag,ephi1
      complex *16 ephi(-ldc-1:ldc+1)
      data imag/(0.0d0,1.0d0)/
C
      rvec(1) = xnynzn(1) - x0y0z0(1)
      rvec(2) = xnynzn(2) - x0y0z0(2)
      rvec(3) = xnynzn(3) - x0y0z0(3)
      call cart2polarl(rvec,d,theta,phi)
c
      ephi1 = exp(imag*phi)
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=dconjg(ephi1)
c
c----- create array of powers e^(i*m*phi).
c
      do l = 1,ldc
         ephi(l+1) = ephi(l)*ephi(1)
         ephi(-1-l) = dconjg(ephi(l+1))
      enddo
c
c      a rotation of THETA radians about the Yprime-axis after PHI
c      radians about the z-axis.
c      The PHI rotation is carried out on the fly by multiplying 
c      locold and ephi inside the following loop. 
c
      do l=0,nterms
         do mp=-l,l
            marray1(l,mp) = locold(l,mp)*ephi(mp)
         enddo
      enddo
      do l=0,nterms2
         do mp=-nterms2,nterms2
            marray(l,mp) = 0.0d0
         enddo
      enddo
ccc      t1 = second()
ccc      call rotviarecur3f90(theta,nterms,nterms,nterms,marray1,
      call rotviarecur3f90(theta,nterms,nterms,nterms2,marray1,
     1      ldc,marray,ldc)
ccc      t2 = second()
c
c----- shift the local expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      rshift = d
ccc      t1 = second()
       call l3dlocloczshift(sc1,marray,ldc,nterms,sc2,local,
     1           nterms2,nterms2,rshift,fr,ier) 
ccc      t2 = second()
c
c      reverse THETA rotation.
c      I.e. rotation of -THETA about Yprime axis.
c
ccc      t1 = second()
      call rotviarecur3f90(-theta,nterms2,nterms2,nterms2,local,
     1      nterms2,marray,ldc)
ccc      t2 = second()
ccc      call prin2(' time for second rot is *',t2-t1,1)
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
         do m=-l,l
            local(l,m)=ephi(-m)*marray(l,m)
         enddo
      enddo
      return
      end
c
c
c
c
c***********************************************************************
      subroutine l3dlocloczshift(scale,locold,lmp,nterms,scale2,
     1  local,lmpn,nterms2,zshift,fr,ier) 
c***********************************************************************
c
c     This subroutine converts a multipole expansion centered at the 
c     origin to a local expansion centered at (0,0,zhift).
c     The expansion is rescaled to that of the local expansion.
c
c     INPUT:
c
c     scale    : scaling parameter for locold
c     locold   : coefficients of original multipole exp.
c     lmp      : leading dim of locold (may be a work array)
c     nterms   : number of terms in the orig. expansion
c
c     scale2   : scaling parameter for output expansion (local)
c     lmpn     : leading dim of local (may be a work array)
c     nterms2  : number of terms in output local exp.
c     zshift   : shifting distance along z-axis (assumed positive)
c     fr       : work array
c
c     OUTPUT:
c
c     local    : coefficients of shifted local exp.
c     ier      : error return code
c                 CURRENTLY NOT USED
c-----------------------------------------------------------------------
      implicit none
      integer nmax,nterms,nterms2,nquad,ier,l,lw,m,jnew,knew
      integer lmp,lmpn,ll
      real *8  zshift,d
      real *8  scale,scale2
      real *8 fr(0:nterms+1)
      complex *16 locold(0:lmp,-lmp:lmp)
      complex *16 local(0:lmpn,-lmpn:lmpn)
      real *8, allocatable :: dc(:,:)
      real *8, allocatable :: carray(:,:)
C
C----- shift along z-axis 
C
      nmax = nterms+nterms2
c
      allocate( dc(0:2*nmax,0:2*nmax) )
      allocate( carray(0:2*nmax,0:2*nmax) )
c
      do l = 0,2*nmax
         carray(l,0) = 1.0d0
         dc(l,0) = 1.0d0
      enddo
      do m=1,2*nmax
         carray(m,m) = 1.0d0
         dc(m,m) = 1.0d0
         do l=m+1,2*nmax
	    carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
	    dc(l,m)=sqrt(carray(l,m))
         enddo
      enddo
c
      d = zshift
      fr(0) = 1.0d0
      d = d*scale
      fr(1) = d
      do l=2,nterms+1
         fr(l) = fr(l-1)*d
      enddo
c
      do jnew = 0,nterms2
         do knew = -jnew,jnew
	    local(jnew,knew) = locold(jnew,knew)
	    do l=1,nterms-jnew
	       ll = l+jnew
	       local(jnew,knew)=local(jnew,knew)+locold(ll,knew)*
     1              fr(l)*dc(ll+knew,l)*dc(ll-knew,l)
            enddo
         enddo
      enddo
      do jnew = 0,nterms2
         do knew = -jnew,jnew
	    local(jnew,knew)=local(jnew,knew)*(scale/scale2)**jnew
         enddo
      enddo
      return
      end
C
