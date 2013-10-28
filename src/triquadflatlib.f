cc Copyright (C) 2009: Leslie Greengard and Zydrunas Gimbutas
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
c    $Date: 2010-10-30 22:11:07 -0400 (Sat, 30 Oct 2010) $
c    $Revision: 1346 $
c
c
c     This file contains quadrature subroutines for constant densities
c     on flat triangles. 
c
c
c     Laplace routines (loop over triangles)
c
c     direct3dtrilaps
c     direct3dtritarglaps
c     direct3dtrilapd
c     direct3dtritarglapd
c
c     Helmholtz and difference kernel routines (loop over triangles)
c
c     direct3dtrihelms2
c     direct3dtritarghelms2
c     direct3dtrihelmd2
c     direct3dtritarghelmd2
c
c     Helmholtz: quadratures for single triangles
c
c     triquadhelm
c     triquadhelmd
c     triquadselfhelm
c     triquadselfhelmd
c
c
c***********************************************************************
c
c                    Laplace SLP routines
c
c***********************************************************************
      subroutine direct3dtrilaps(ipatch,ntri,
     1           zparts,charge,triang,rpot,ptfrc)
c***********************************************************************
c
c     computes potential and field at centroid of patch ipatch due to 
c     piecewise-constant charge density on collection of triangles 
c     numbered jpatch = 1,...,ntri. 
c     If ipatch equals jpatch, they are assumed to be the same triangle
c     and a singular quadrature rule is used. Otherwise, they are assumed 
c     to be distinct. In either case, analytic quadratures are used
c     (see triahquad.f)
c     
c     INPUT:
c
c     ipatch            target face (patch) 
c     ntri              number of triangles
c     zparts(3,ntri)    array of triangle centroids
c     charge(ntri)      array of (piecewise constant) SLP strengths
c     triang(3,3,ntri)  array of triangles in standard format
c
c     OUPUT:
c
c     rpot              potential at centroid zparts(*,ipatch)
c     ptfrc(3)          -gradient of potential at centroid
c
c
      implicit none
      integer ifinit, ipatch, jpatch, ntri, nquad, ier
      integer itype, iquad
      real *8   zparts(3,1), triang(3,3,1)
      real *8   point(3)
      real *8   vert1(2), vert2(2), vert3(2),  vertout(3)
      real *8   w(12)
      real *8   x0,y0,z0,val,valx,valy,valz,derx,dery,derz
      complex *16 charge(1), ptfrc(3), rpot
c
      rpot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c  
      do jpatch = 1, ntri
         call tri_ini(triang(1,1,jpatch),triang(1,2,jpatch),
     1                triang(1,3,jpatch),w,vert1,vert2,vert3)
c
         call tri_for(w,zparts(1,ipatch),vertout)
         x0 = vertout(1)
         y0 = vertout(2)
         z0 = vertout(3)
c
         if (ipatch.eq.jpatch) then
            itype = 1
            iquad = 0
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
            itype = 2
            iquad = 0
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
            itype = 3
            iquad = 0
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
            itype = 4
            iquad = 0
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valz)
            valz = -valz
         else
            iquad = 0
            if (z0.gt.0) iquad = 1
            if (z0.lt.0) iquad = -1
            itype = 1
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
            itype = 2
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
            itype = 3
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
            itype = 4
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valz)
            valz = -valz
         endif
c
         rpot = rpot + charge(jpatch)*val
c
         call rotder3d(w,triang(1,1,jpatch),
     $      valx,valy,valz,derx,dery,derz)
         ptfrc(1)=ptfrc(1)-charge(jpatch)*derx
         ptfrc(2)=ptfrc(2)-charge(jpatch)*dery
         ptfrc(3)=ptfrc(3)-charge(jpatch)*derz
c
      enddo
      return
      end
c
c
c
c
c
c***********************************************************************
      subroutine direct3dtrilaps2(ipatch,ntri,
     1           zparts,charge,triang,ifpot,rpot,iffld,ptfrc)
c***********************************************************************
c
c     computes potential and field at centroid of patch ipatch due to 
c     piecewise-constant charge density on collection of triangles 
c     numbered jpatch = 1,...,ntri. 
c     If ipatch equals jpatch, they are assumed to be the same triangle
c     and a singular quadrature rule is used. Otherwise, they are assumed 
c     to be distinct. In either case, analytic quadratures are used
c     (see triahquad.f)
c     
c     INPUT:
c
c     ipatch            target face (patch) 
c     ntri              number of triangles
c     zparts(3,ntri)    array of triangle centroids
c     charge(ntri)      array of (piecewise constant) SLP strengths
c     triang(3,3,ntri)  array of triangles in standard format
c
c     OUPUT:
c
c     rpot              potential at centroid zparts(*,ipatch)
c     ptfrc(3)          -gradient of potential at centroid
c
c
      implicit none
      integer ifpot,iffld
      integer ifinit, ipatch, jpatch, ntri, nquad, ier
      integer itype, iquad
      real *8   zparts(3,1), triang(3,3,1)
      real *8   point(3)
      real *8   vert1(2), vert2(2), vert3(2),  vertout(3)
      real *8   w(12)
      real *8   x0,y0,z0,val,valx,valy,valz,derx,dery,derz
      complex *16 charge(1), ptfrc(3), rpot
c
      rpot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c  
      do jpatch = 1, ntri
         call tri_ini(triang(1,1,jpatch),triang(1,2,jpatch),
     1                triang(1,3,jpatch),w,vert1,vert2,vert3)
c
         call tri_for(w,zparts(1,ipatch),vertout)
         x0 = vertout(1)
         y0 = vertout(2)
         z0 = vertout(3)
c
         if (ipatch.eq.jpatch) then
            iquad = 0
            if( ifpot .eq. 1 ) then
            itype = 1
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
            endif
            if( iffld .eq. 1 ) then
            itype = 2
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
            itype = 3
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
            itype = 4
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valz)
            valz = -valz
            endif
         else
            iquad = 0
            if (z0.gt.0) iquad = 1
            if (z0.lt.0) iquad = -1
            if( ifpot .eq. 1 ) then
            itype = 1
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
            endif
            if( iffld .eq. 1 ) then
            itype = 2
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
            itype = 3
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
            itype = 4
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valz)
            valz = -valz
            endif
         endif
c
         if( ifpot .eq. 1 ) then
         rpot = rpot + charge(jpatch)*val
         endif
c
         if( iffld .eq. 1 ) then
         call rotder3d(w,triang(1,1,jpatch),
     $      valx,valy,valz,derx,dery,derz)
         ptfrc(1)=ptfrc(1)-charge(jpatch)*derx
         ptfrc(2)=ptfrc(2)-charge(jpatch)*dery
         ptfrc(3)=ptfrc(3)-charge(jpatch)*derz
         endif
c
      enddo
      return
      end
c
c
c
c
c
c***********************************************************************
      subroutine direct3dtritarglaps(ntri,targ,charge,triang,rpot,
     1           ptfrc)
c***********************************************************************
c
c     computes potential and field at arbitrary point TARG not
c     lying on the surface due to piecewise-constant charge density
c     on collection of triangles.
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     ntri              number of triangles
c     targ(3)           target location
c     charge(ntri)      array of SLP strengths (constant)
c     triang(3,3,ntri)  array of triangles in standard format
c
c     OUPUT:
c
c     rpot              potential at targ
c     ptfrc(3)          -gradient of potential at targ
c
c
c
      implicit none
      integer ifinit, jpatch, ntri, nquad, ier,j
      integer iquad, itype
      real *8   targ(3), triang(3,3,1)
      real *8   point(3)
      real *8   vert1(2), vert2(2), vert3(2),  vertout(3)
      real *8   w(12)
      real *8   x0,y0,z0, val, valx,valy,valz,derx,dery,derz
      complex *16 charge(1), ptfrc(3), rpot
c
      rpot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c  
      do jpatch = 1, ntri
         call tri_ini(triang(1,1,jpatch),triang(1,2,jpatch),
     1                triang(1,3,jpatch),w,vert1,vert2,vert3)
         call tri_for(w,targ,vertout)
         x0 = vertout(1)
         y0 = vertout(2)
         z0 = vertout(3)
c
         iquad = 0
         if (z0.gt.0) iquad = 1
         if (z0.lt.0) iquad = -1
         itype = 1
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      val)
         itype = 2
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valx)
         itype = 3
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valy)
         itype = 4
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valz)
         valz = -valz
c
         rpot = rpot + charge(jpatch)*val
c
         call rotder3d(w,triang(1,1,jpatch),
     $      valx,valy,valz,derx,dery,derz)
         ptfrc(1)=ptfrc(1)-charge(jpatch)*derx
         ptfrc(2)=ptfrc(2)-charge(jpatch)*dery
         ptfrc(3)=ptfrc(3)-charge(jpatch)*derz
c
      enddo
      return
      end
c
c
c
c
c***********************************************************************
      subroutine direct3dtritarglaps2(ntri,targ,charge,triang,
     $     ifpot,rpot,iffld,ptfrc)
c***********************************************************************
c
c     computes potential and field at arbitrary point TARG not
c     lying on the surface due to piecewise-constant charge density
c     on collection of triangles.
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     ntri              number of triangles
c     targ(3)           target location
c     charge(ntri)      array of SLP strengths (constant)
c     triang(3,3,ntri)  array of triangles in standard format
c
c     OUPUT:
c
c     rpot              potential at targ
c     ptfrc(3)          -gradient of potential at targ
c
c
c
      implicit none
      integer ifpot,iffld
      integer ifinit, jpatch, ntri, nquad, ier,j
      integer iquad, itype
      real *8   targ(3), triang(3,3,1)
      real *8   point(3)
      real *8   vert1(2), vert2(2), vert3(2),  vertout(3)
      real *8   w(12)
      real *8   x0,y0,z0, val, valx,valy,valz,derx,dery,derz
      complex *16 charge(1), ptfrc(3), rpot
c
      rpot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c  
      do jpatch = 1, ntri
         call tri_ini(triang(1,1,jpatch),triang(1,2,jpatch),
     1                triang(1,3,jpatch),w,vert1,vert2,vert3)
         call tri_for(w,targ,vertout)
         x0 = vertout(1)
         y0 = vertout(2)
         z0 = vertout(3)
c
         iquad = 0
         if (z0.gt.0) iquad = 1
         if (z0.lt.0) iquad = -1
         if( ifpot .eq. 1) then
         itype = 1
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      val)
         endif
         if( iffld .eq. 1) then
         itype = 2
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valx)
         itype = 3
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valy)
         itype = 4
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valz)
         valz = -valz
         endif
c
         if( ifpot .eq. 1) then
         rpot = rpot + charge(jpatch)*val
         endif
c
         if( iffld .eq. 1) then
         call rotder3d(w,triang(1,1,jpatch),
     $      valx,valy,valz,derx,dery,derz)
         ptfrc(1)=ptfrc(1)-charge(jpatch)*derx
         ptfrc(2)=ptfrc(2)-charge(jpatch)*dery
         ptfrc(3)=ptfrc(3)-charge(jpatch)*derz
         endif
c
      enddo
      return
      end
c
c
c***********************************************************************
c
c                 Laplace DLP routines
c
c***********************************************************************
      subroutine direct3dtrilapd(ipatch,ntri,
     1           zparts,dipstr,triang,trinorm,rpot,ptfrc)
c***********************************************************************
c
c     computes potential and field at centroid of patch ipatch due to 
c     piecewise constant dipole layer on collection of triangles.
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     ipatch            target face (patch) 
c     ntri              number of triangles
c     zparts(3,ntri)    array of triangle centroids
c     dipstr(ntri)      array of (piecewise constant) DLP strengths
c     triang(3,3,ntri)  array of triangles in standard format
c     trinorm(3,ntri)   array of triangle normals
c
c     OUPUT:
c
c     rpot              potential at centroid zparts(*,ipatch)
c     ptfrc(3)          -gradient of potential at centroid
c
c
      implicit none
      integer ifinit, ipatch, jpatch, ntri, nquad, ier
      integer itype, iquad
      real *8   zparts(3,1), triang(3,3,1)
      real *8   trinorm(3,1)
      real *8   point(3)
      real *8   vert1(2), vert2(2), vert3(2),  vertout(3)
      real *8   w(12)
      real *8   x0,y0,z0,val,valx,valy,valz,derx,dery,derz
      complex *16 dipstr(1), ptfrc(3), rpot
c
      rpot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c  
      do jpatch = 1, ntri
         call tri_ini(triang(1,1,jpatch),triang(1,2,jpatch),
     1                triang(1,3,jpatch),w,vert1,vert2,vert3)
c
         call tri_for(w,zparts(1,ipatch),vertout)
         x0 = vertout(1)
         y0 = vertout(2)
         z0 = vertout(3)
         x0 = x0
         y0 = y0
         z0 = -z0
c
         if (ipatch.eq.jpatch) then
            itype = 4
            iquad = 0
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
            val = -val
            itype = 5
            iquad = 0
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
            itype = 6
            iquad = 0
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
            itype = 7
            iquad = 0
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valz)
         else
            iquad = 0
            if (z0.gt.0) iquad = +1
            if (z0.lt.0) iquad = -1
            itype = 4
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
            val = -val
            itype = 5
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
            itype = 6
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
            itype = 7
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valz)
         endif

         rpot = rpot + dipstr(jpatch)*val
         call rotder3d(w,triang(1,1,jpatch),
     $      valx,valy,valz,derx,dery,derz)
         ptfrc(1)=ptfrc(1)+dipstr(jpatch)*derx
         ptfrc(2)=ptfrc(2)+dipstr(jpatch)*dery
         ptfrc(3)=ptfrc(3)+dipstr(jpatch)*derz
c
      enddo 
      return
      end
c
c
c
c
c
      subroutine direct3dtrilapd2(ipatch,ntri,
     1     zparts,dipstr,triang,trinorm,
     $     ifpot,rpot,iffld,ptfrc)
c***********************************************************************
c
c     computes potential and field at centroid of patch ipatch due to 
c     piecewise constant dipole layer on collection of triangles.
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     ipatch            target face (patch) 
c     ntri              number of triangles
c     zparts(3,ntri)    array of triangle centroids
c     dipstr(ntri)      array of (piecewise constant) DLP strengths
c     triang(3,3,ntri)  array of triangles in standard format
c     trinorm(3,ntri)   array of triangle normals
c
c     OUPUT:
c
c     rpot              potential at centroid zparts(*,ipatch)
c     ptfrc(3)          -gradient of potential at centroid
c
c
      implicit none
      integer ifpot,iffld
      integer ifinit, ipatch, jpatch, ntri, nquad, ier
      integer itype, iquad
      real *8   zparts(3,1), triang(3,3,1)
      real *8   trinorm(3,1)
      real *8   point(3)
      real *8   vert1(2), vert2(2), vert3(2),  vertout(3)
      real *8   w(12)
      real *8   x0,y0,z0,val,valx,valy,valz,derx,dery,derz
      complex *16 dipstr(1), ptfrc(3), rpot
c
      rpot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c  
      do jpatch = 1, ntri
         call tri_ini(triang(1,1,jpatch),triang(1,2,jpatch),
     1                triang(1,3,jpatch),w,vert1,vert2,vert3)
c
         call tri_for(w,zparts(1,ipatch),vertout)
         x0 = vertout(1)
         y0 = vertout(2)
         z0 = vertout(3)
         x0 = x0
         y0 = y0
         z0 = -z0
c
         if (ipatch.eq.jpatch) then
            iquad = 0
            if( ifpot .eq. 1 ) then
            itype = 4
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
            val = -val
            endif
            if( iffld .eq. 1 ) then
            itype = 5
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
            itype = 6
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
            itype = 7
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valz)
            endif
         else
            iquad = 0
            if (z0.gt.0) iquad = +1
            if (z0.lt.0) iquad = -1
            if( ifpot .eq. 1 ) then
            itype = 4
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
            val = -val
            endif
            if( iffld .eq. 1 ) then
            itype = 5
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
            itype = 6
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
            itype = 7
            call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valz)
            endif
         endif

         if( ifpot .eq. 1 ) then
         rpot = rpot + dipstr(jpatch)*val
         endif

         if( iffld .eq. 1 ) then
         call rotder3d(w,triang(1,1,jpatch),
     $      valx,valy,valz,derx,dery,derz)
         ptfrc(1)=ptfrc(1)+dipstr(jpatch)*derx
         ptfrc(2)=ptfrc(2)+dipstr(jpatch)*dery
         ptfrc(3)=ptfrc(3)+dipstr(jpatch)*derz
         endif
c
      enddo 
      return
      end
c
c
c
c
c
c***********************************************************************
      subroutine direct3dtritarglapd(ntri,targ,dipstr,triang,trinorm,
     1           rpot,ptfrc)
c***********************************************************************
c
c     computes potential and field at arbitrary point TARG due to 
c     piecewise constant dipole layer on collection of triangles.
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     ntri              number of triangles
c     targ(3)           target location
c     dipstr(ntri)      array of SLP strengths (constant)
c     triang(3,3,ntri)  array of triangles in standard format
c     trinorm(3,ntri)   array of triangle normals
c
c     OUPUT:
c
c     rpot              potential at targ
c     ptfrc(3)          -gradient of potential at targ
c
c
      implicit none
      integer ifinit, jpatch, ntri, nquad, ier,j
      integer iquad, itype
      real *8   targ(3), triang(3,3,1), trinorm(3,1)
      real *8   point(3)
      real *8   vert1(2), vert2(2), vert3(2),  vertout(3)
      real *8   w(12)
      real *8   x0,y0,z0, val, valx,valy,valz,derx,dery,derz
        complex *16 dipstr(1), ptfrc(3), rpot
c
      rpot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c  
      do jpatch = 1, ntri
         call tri_ini(triang(1,1,jpatch),triang(1,2,jpatch),
     1                triang(1,3,jpatch),w,vert1,vert2,vert3)
         call tri_for(w,targ,vertout)
         x0 = vertout(1)
         y0 = vertout(2)
         z0 = vertout(3)
         x0 = x0
         y0 = y0
         z0 = -z0
c
         iquad = 0
         if (z0.gt.0) iquad = +1
         if (z0.lt.0) iquad = -1
         itype = 4
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      val)
         val = -val
         itype = 5
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valx)
         itype = 6
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valy)
         itype = 7
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valz)

         rpot = rpot + dipstr(jpatch)*val
         call rotder3d(w,triang(1,1,jpatch),
     $      valx,valy,valz,derx,dery,derz)
         ptfrc(1)=ptfrc(1)+dipstr(jpatch)*derx
         ptfrc(2)=ptfrc(2)+dipstr(jpatch)*dery
         ptfrc(3)=ptfrc(3)+dipstr(jpatch)*derz
      enddo
      return
      end
c
c
c***********************************************************************
      subroutine direct3dtritarglapd2(ntri,targ,dipstr,triang,trinorm,
     1           ifpot,rpot,iffld,ptfrc)
c***********************************************************************
c
c     computes potential and field at arbitrary point TARG due to 
c     piecewise constant dipole layer on collection of triangles.
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     ntri              number of triangles
c     targ(3)           target location
c     dipstr(ntri)      array of SLP strengths (constant)
c     triang(3,3,ntri)  array of triangles in standard format
c     trinorm(3,ntri)   array of triangle normals
c
c     OUPUT:
c
c     rpot              potential at targ
c     ptfrc(3)          -gradient of potential at targ
c
c
      implicit none
      integer ifpot,iffld
      integer ifinit, jpatch, ntri, nquad, ier,j
      integer iquad, itype
      real *8   targ(3), triang(3,3,1), trinorm(3,1)
      real *8   point(3)
      real *8   vert1(2), vert2(2), vert3(2),  vertout(3)
      real *8   w(12)
      real *8   x0,y0,z0, val, valx,valy,valz,derx,dery,derz
        complex *16 dipstr(1), ptfrc(3), rpot
c
      rpot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c  
      do jpatch = 1, ntri
         call tri_ini(triang(1,1,jpatch),triang(1,2,jpatch),
     1                triang(1,3,jpatch),w,vert1,vert2,vert3)
         call tri_for(w,targ,vertout)
         x0 = vertout(1)
         y0 = vertout(2)
         z0 = vertout(3)
         x0 = x0
         y0 = y0
         z0 = -z0
c
         iquad = 0
         if (z0.gt.0) iquad = +1
         if (z0.lt.0) iquad = -1
c
         if( ifpot .eq. 1 ) then
         itype = 4
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      val)
         val = -val
         endif
         if( iffld .eq. 1 ) then
         itype = 5
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valx)
         itype = 6
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valy)
         itype = 7
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      valz)
         endif

         if( ifpot .eq. 1 ) then
         rpot = rpot + dipstr(jpatch)*val
         endif

         if( iffld .eq. 1 ) then
         call rotder3d(w,triang(1,1,jpatch),
     $      valx,valy,valz,derx,dery,derz)
         ptfrc(1)=ptfrc(1)+dipstr(jpatch)*derx
         ptfrc(2)=ptfrc(2)+dipstr(jpatch)*dery
         ptfrc(3)=ptfrc(3)+dipstr(jpatch)*derz
         endif

      enddo
      return
      end
c
c
c***********************************************************************
c
c              Helmholtz/difference SLP routines
c
c***********************************************************************
      subroutine direct3dtrihelms2(ifl,ipatch,ntri,zk,nquad,
     1           zparts,charge,triang,pot,ptfrc)
c***********************************************************************
c
c     computes Helmholtz (or difference) potential and field at 
c     centroid of patch ipatch due to piecewise-constant charge density 
c     on collection of triangles numbered jpatch = 1,...,ntri. 
c     If ipatch equals jpatch, they are assumed to be the same triangle
c     and a singular quadrature rule is used. Otherwise, they are assumed 
c     to be distinct. 
c     
c
c     INPUT:
c
c     ifl               kernel flag.
c                       ifl = 0  means Helmholtz kernel
c                       ifl = 1  means difference kernel
c                                (Helmholtz-Laplace)
c     ipatch            target face (patch) 
c     ntri              number of triangles
c     zk                Helmholtz parameter
c     zparts(3,ntri)    array of triangle centroids
c     charge(ntri)      array of (piecewise constant) SLP strengths
c     triang(3,3,ntri)  array of triangles in standard format
c
c     OUPUT:
c
c     pot              potential at centroid zparts(*,ipatch)
c     ptfrc(3)         -gradient of potential at centroid
c
c
c     Number of quadrature nodes used to integrate diff of Greens
c     functions on triangles is set here as nquad, followed
c     by a call to legewhts. wts, xnodes, bb are dimensioned as
c     100 - typical values are nquad = 10 for about 4 digits,
c     nquad = 20 for about 8 digits, nquad = 30 for about 12 digits.
c
c
      implicit none
      integer  ifl, ifinit, ipatch, jpatch, ntri, nquad, ier
      integer  iquad, itype, nnodes
      real *8  zparts(3,1), triang(3,3,1)
      real *8  point(3)
      real *8  vert1(2), vert2(2), vert3(2),  vertout(3)
      real *8  w(12), val
      real *8  x0,y0,z0
      real *8  wts(100),xnodes(100),bb(100)
      real *8  weights(500),rnodes(2,500)
      complex *16  zk,pot,ptfrc(3),charge(1)
      complex *16   zval,zvalx,zvaly,zvalz,zderx,zdery,zderz
c
ccc      nquad = 20
ccc      ifinit = 1
ccc      call legewhts(nquad,xnodes,wts,ifinit)
      pot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c  
      do jpatch = 1, ntri
         call tri_ini(triang(1,1,jpatch),triang(1,2,jpatch),
     1                triang(1,3,jpatch),w,vert1,vert2,vert3)
c
         call triasymq(nquad,vert1,vert2,vert3,rnodes,weights,
     1                nnodes)
c
         call tri_for(w,zparts(1,ipatch),vertout)
         x0 = vertout(1)
         y0 = vertout(2)
         z0 = vertout(3)
c
         if (ipatch.eq.jpatch) then
            call triquadselfhelm(ifl,vert1,vert2,vert3,x0,y0,zk,
     1   	   zval,zvalx,zvaly,zvalz,nnodes,weights,rnodes,ier)
         else
            call triquadhelm(ifl,vert1,vert2,vert3,x0,y0,z0,zk,
     1   	   zval,zvalx,zvaly,zvalz,nnodes,weights,rnodes,ier)
         endif
         pot = pot + charge(jpatch)*zval
c
         call rotder3dz(w,
     $      zvalx,zvaly,zvalz,zderx,zdery,zderz)
         ptfrc(1)=ptfrc(1)-charge(jpatch)*zderx
         ptfrc(2)=ptfrc(2)-charge(jpatch)*zdery
         ptfrc(3)=ptfrc(3)-charge(jpatch)*zderz
c
      enddo
      return
      end
c
c
c
c
c
c***********************************************************************
      subroutine direct3dtritarghelms2(ifl,ntri,targ,zk,nquad,
     1           charge,triang,pot,ptfrc)
c***********************************************************************
c
c     computes Helmholtz (or difference) potential and field at 
c     arbitrary point TARG due to piecewise-constant charge density 
c     on collection of triangles numbered jpatch = 1,...,ntri. 
c     
c
c     INPUT:
c
c     ifl               kernel flag.
c                       ifl = 0  means Helmholtz kernel
c                       ifl = 1  means difference kernel
c                                (Helmholtz-Laplace)
c     ntri              number of triangles
c     targ(3)           target location
c     zk                Helmholtz parameter
c     charge(ntri)      array of SLP strengths (constant)
c     triang(3,3,ntri)  array of triangles in standard format
c
c     OUPUT:
c
c     pot              potential at targ
c     ptfrc(3)         -gradient of potential at targ
c
c
c     Number of quadrature nodes used to integrate Greens
c     function on triangles is set here as nquad, followed
c     by a call to gaussq. wts, xnodes, bb are dimensioned as
c     100 - typical values are nquad = 10 for 4 digits,
c     nquad = 20 for 8 digits, nquad = 30 for 12 digits.
c
c
      implicit none
      integer ifl, ifinit, jpatch, ntri, nquad, ier,j, nnodes
      real *8   targ(3), triang(3,3,1)
      real *8   point(3)
      real *8   vert1(2), vert2(2), vert3(2),  vertout(3)
      real *8   w(12)
      real *8   x0,y0,z0
      real *8  wts(100),xnodes(100),bb(100)
      real *8  weights(500),rnodes(2,500)
      complex *16  zk, pot,ptfrc(3),charge(1)
      complex *16  zval, zvalx,zvaly,zvalz,zderx,zdery,zderz
c
ccc      nquad = 20
ccc      ifinit = 1
ccc      call legewhts(nquad,xnodes,wts,ifinit)
      pot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c  
      do jpatch = 1, ntri
         call tri_ini(triang(1,1,jpatch),triang(1,2,jpatch),
     1                triang(1,3,jpatch),w,vert1,vert2,vert3)
         call triasymq(nquad,vert1,vert2,vert3,rnodes,weights,
     1                nnodes)
c
         call tri_for(w,targ,vertout)
         x0 = vertout(1)
         y0 = vertout(2)
         z0 = vertout(3)
         call triquadhelm(ifl,vert1,vert2,vert3,x0,y0,z0,zk,
     1   	   zval,zvalx,zvaly,zvalz,nnodes,weights,rnodes,ier)
         pot = pot + charge(jpatch)*zval
         call rotder3dz(w,
     $      zvalx,zvaly,zvalz,zderx,zdery,zderz)
         ptfrc(1)=ptfrc(1)-charge(jpatch)*zderx
         ptfrc(2)=ptfrc(2)-charge(jpatch)*zdery
         ptfrc(3)=ptfrc(3)-charge(jpatch)*zderz
      enddo
      return
      end
c
c***********************************************************************
      subroutine direct3dtritarghelms3(ifl,ntri,targ,zk,nquad,
     1           charge,triang,ifpot,pot,iffld,ptfrc)
c***********************************************************************
c
c     computes Helmholtz (or difference) potential and field at 
c     arbitrary point TARG due to piecewise-constant charge density 
c     on collection of triangles numbered jpatch = 1,...,ntri. 
c     
c
c     INPUT:
c
c     ifl               kernel flag.
c                       ifl = 0  means Helmholtz kernel
c                       ifl = 1  means difference kernel
c                                (Helmholtz-Laplace)
c     ntri              number of triangles
c     targ(3)           target location
c     zk                Helmholtz parameter
c     charge(ntri)      array of SLP strengths (constant)
c     triang(3,3,ntri)  array of triangles in standard format
c
c     OUPUT:
c
c     pot              potential at targ
c     ptfrc(3)         -gradient of potential at targ
c
c
c     Number of quadrature nodes used to integrate Greens
c     function on triangles is set here as nquad, followed
c     by a call to gaussq. wts, xnodes, bb are dimensioned as
c     100 - typical values are nquad = 10 for 4 digits,
c     nquad = 20 for 8 digits, nquad = 30 for 12 digits.
c
c
      implicit none
      integer ifl, ifinit, jpatch, ntri, nquad, ier,j, nnodes
      real *8   targ(3), triang(3,3,1)
      real *8   point(3)
      real *8   vert1(2), vert2(2), vert3(2),  vertout(3)
      real *8   w(12)
      real *8   x0,y0,z0
      real *8  wts(100),xnodes(100),bb(100)
      real *8  weights(500),rnodes(2,500)
      complex *16  zk, pot,ptfrc(3),charge(1)
      complex *16  zval, zvalx,zvaly,zvalz,zderx,zdery,zderz
      integer ifpot,iffld
c
ccc      nquad = 20
ccc      ifinit = 1
ccc      call legewhts(nquad,xnodes,wts,ifinit)
      pot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c  
      do jpatch = 1, ntri
         call tri_ini(triang(1,1,jpatch),triang(1,2,jpatch),
     1                triang(1,3,jpatch),w,vert1,vert2,vert3)
         call triasymq(nquad,vert1,vert2,vert3,rnodes,weights,
     1                nnodes)
c
         call tri_for(w,targ,vertout)
         x0 = vertout(1)
         y0 = vertout(2)
         z0 = vertout(3)
         call triquadhelm2(ifl,vert1,vert2,vert3,x0,y0,z0,zk,
     1      ifpot,iffld,zval,zvalx,zvaly,zvalz,
     $      nnodes,weights,rnodes,ier)
         if( ifpot .eq. 1 ) then
         pot = pot + charge(jpatch)*zval
         endif
         if( iffld .eq. 1 ) then
         call rotder3dz(w,
     $      zvalx,zvaly,zvalz,zderx,zdery,zderz)
         ptfrc(1)=ptfrc(1)-charge(jpatch)*zderx
         ptfrc(2)=ptfrc(2)-charge(jpatch)*zdery
         ptfrc(3)=ptfrc(3)-charge(jpatch)*zderz
         endif
      enddo
      return
      end
c
c***********************************************************************
c
c                Helmholtz/difference DLP routines.
c
c***********************************************************************
      subroutine direct3dtrihelmd2(ifl,ipatch,ntri,zk,nquad,
     1           zparts,dipstr,triang,trinorm,pot,ptfrc)
c***********************************************************************
c
c     computes Helmholtz (or difference) potential and field at centroid 
c     of patch ipatch due to piecewise-constant dipole density on 
c     collection of triangles numbered jpatch = 1,...,ntri. 
c     If ipatch equals jpatch, they are assumed to be the same triangle
c     and a singular quadrature rule is used. Otherwise, they are assumed 
c     to be distinct. 
c     
c     INPUT:
c
c     ifl               kernel flag.
c                       ifl = 0  means Helmholtz kernel
c                       ifl = 1  means difference kernel
c                                (Helmholtz-Laplace)
c     ipatch            target face (patch) 
c     ntri              number of triangles
c     zk                Helmholtz parameter
c     zparts(3,ntri)    array of triangle centroids
c     dipstr(ntri)      array of (piecewise constant) DLP strengths
c     triang(3,3,ntri)  array of triangles in standard format
c     trinorm(3,ntri)   array of triangle normals
c
c     OUPUT:
c
c     pot              potential at centroid zparts(*,ipatch)
c     ptfrc(3)         -gradient of potential at centroid
c
c
c     Number of quadrature nodes used to integrate diff of Greens
c     functions on triangles is set here as nquad, followed
c     by a call to legewhts. wts, xnodes, bb are dimensioned as
c     100 - typical values are nquad = 10 for about 4 digits,
c     nquad = 20 for about 8 digits, nquad = 30 for about 12 digits.
c
c
      implicit none
      integer ifl, ifinit, ipatch, jpatch, ntri, nquad, ier, nnodes
      real *8   zparts(3,1), triang(3,3,1), trinorm(3,1)
      real *8   point(3)
      real *8   vert1(2), vert2(2), vert3(2),  vertout(3)
      real *8   w(12)
      real *8   x0,y0,z0
      real *8  wts(100),xnodes(100),bb(100)
      real *8  weights(500),rnodes(2,500)
      complex *16  zk,pot,ptfrc(3),dipstr(1)
      complex *16   zval,zvalx,zvaly,zvalz,zderx,zdery,zderz
c
ccc      nquad = 20
ccc      ifinit = 1
ccc      call legewhts(nquad,xnodes,wts,ifinit)
      pot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c  
      do jpatch = 1, ntri
         call tri_ini(triang(1,1,jpatch),triang(1,2,jpatch),
     1                triang(1,3,jpatch),w,vert1,vert2,vert3)
c
         call triasymq(nquad,vert1,vert2,vert3,rnodes,weights,
     1                nnodes)
         call tri_for(w,zparts(1,ipatch),vertout)
         x0 = vertout(1)
         y0 = vertout(2)
         z0 = vertout(3)
c
         if (ipatch.eq.jpatch) then
            call triquadselfhelmd(ifl,vert1,vert2,vert3,x0,y0,zk,
     1   	   zval,zvalx,zvaly,zvalz,nnodes,weights,rnodes,ier)
         else
            call triquadhelmd(ifl,vert1,vert2,vert3,x0,y0,z0,zk,
     1   	   zval,zvalx,zvaly,zvalz,nnodes,weights,rnodes,ier)
         endif
c
c       rotate the derivative vector, normal derivative is (0,0,1)
c       in the rotated coordinate system
c
         pot = pot - dipstr(jpatch)*zval
c
         call rotder3dz(w,
     $      zvalx,zvaly,zvalz,zderx,zdery,zderz)
c
         ptfrc(1) = ptfrc(1) + dipstr(jpatch)*zderx
         ptfrc(2) = ptfrc(2) + dipstr(jpatch)*zdery
         ptfrc(3) = ptfrc(3) + dipstr(jpatch)*zderz
c
      enddo
      return
      end
c
c
c
c***********************************************************************
      subroutine direct3dtritarghelmd2(ifl,ntri,targ,zk,nquad,
     1           dipstr,triang,trinorm,pot,ptfrc)
c***********************************************************************
c
c     computes Helmholtz (or difference) potential and field at arbitrary 
c     point TARG due to piecewise-constant dipole density on collection 
c     of triangles numbered jpatch = 1,...,ntri. 
c     
c     INPUT:
c
c     ifl               kernel flag.
c                       ifl = 0  means Helmholtz kernel
c                       ifl = 1  means difference kernel
c                                (Helmholtz-Laplace)
c     ntri              number of triangles
c     targ(3)           target location
c     zk                Helmholtz parameter
c     dipstr(ntri)      array of SLP strengths (constant)
c     triang(3,3,ntri)  array of triangles in standard format
c     trinorm(3,ntri)   array of triangle normals
c
c     OUPUT:
c
c     pot              potential at targ
c     ptfrc(3)         -gradient of potential at targ
c
c
c     Number of quadrature nodes used to integrate diff of Greens
c     functions on triangles is set here as nquad, followed
c     by a call to legewhts. wts, xnodes, bb are dimensioned as
c     100 - typical values are nquad = 10 for 4 digits,
c     nquad = 20 for 8 digits, nquad = 30 for 12 digits.
c
c
      implicit none
      integer ifl, ifinit, jpatch, ntri, nquad, ier,j, nnodes
      real *8   targ(3), triang(3,3,1), trinorm(3,1)
      real *8   point(3)
      real *8   vert1(2), vert2(2), vert3(2),  vertout(3)
      real *8   w(12)
      real *8   x0,y0,z0
      real *8  wts(100),xnodes(100),bb(100)
      real *8  weights(500),rnodes(2,500)
      complex *16  zk, pot,ptfrc(3),dipstr(1)
      complex *16  zval, zvalx,zvaly,zvalz,zderx,zdery,zderz
c
ccc      nquad = 20
ccc      ifinit = 1
ccc      call legewhts(nquad,xnodes,wts,ifinit)
      pot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c  
      do jpatch = 1, ntri
         call tri_ini(triang(1,1,jpatch),triang(1,2,jpatch),
     1                triang(1,3,jpatch),w,vert1,vert2,vert3)
         call triasymq(nquad,vert1,vert2,vert3,rnodes,weights,
     1                nnodes)
c
         call tri_for(w,targ,vertout)
         x0 = vertout(1)
         y0 = vertout(2)
         z0 = vertout(3)
         call triquadhelmd(ifl,vert1,vert2,vert3,x0,y0,z0,zk,
     1   	   zval,zvalx,zvaly,zvalz,nnodes,weights,rnodes,ier)
c
c       rotate the derivative vector, normal derivative is (0,0,1)
c       in the rotated coordinate system
c
         pot = pot - dipstr(jpatch)*zval
c
         call rotder3dz(w,
     $      zvalx,zvaly,zvalz,zderx,zdery,zderz)
c
         ptfrc(1) = ptfrc(1) + dipstr(jpatch)*zderx
         ptfrc(2) = ptfrc(2) + dipstr(jpatch)*zdery
         ptfrc(3) = ptfrc(3) + dipstr(jpatch)*zderz
c
      enddo
      return
      end
c
c
c
c
c***********************************************************************
      subroutine direct3dtritarghelmd3(ifl,ntri,targ,zk,nquad,
     1           dipstr,triang,trinorm,ifpot,pot,iffld,ptfrc)
c***********************************************************************
c
c     computes Helmholtz (or difference) potential and field at arbitrary 
c     point TARG due to piecewise-constant dipole density on collection 
c     of triangles numbered jpatch = 1,...,ntri. 
c     
c     INPUT:
c
c     ifl               kernel flag.
c                       ifl = 0  means Helmholtz kernel
c                       ifl = 1  means difference kernel
c                                (Helmholtz-Laplace)
c     ntri              number of triangles
c     targ(3)           target location
c     zk                Helmholtz parameter
c     dipstr(ntri)      array of SLP strengths (constant)
c     triang(3,3,ntri)  array of triangles in standard format
c     trinorm(3,ntri)   array of triangle normals
c
c     OUPUT:
c
c     pot              potential at targ
c     ptfrc(3)         -gradient of potential at targ
c
c
c     Number of quadrature nodes used to integrate diff of Greens
c     functions on triangles is set here as nquad, followed
c     by a call to legewhts. wts, xnodes, bb are dimensioned as
c     100 - typical values are nquad = 10 for 4 digits,
c     nquad = 20 for 8 digits, nquad = 30 for 12 digits.
c
c
      implicit none
      integer ifl, ifinit, jpatch, ntri, nquad, ier,j, nnodes
      real *8   targ(3), triang(3,3,1), trinorm(3,1)
      real *8   point(3)
      real *8   vert1(2), vert2(2), vert3(2),  vertout(3)
      real *8   w(12)
      real *8   x0,y0,z0
      real *8  wts(100),xnodes(100),bb(100)
      real *8  weights(500),rnodes(2,500)
      complex *16  zk, pot,ptfrc(3),dipstr(1)
      complex *16  zval, zvalx,zvaly,zvalz,zderx,zdery,zderz
      integer ifpot,iffld
c
ccc      nquad = 20
ccc      ifinit = 1
ccc      call legewhts(nquad,xnodes,wts,ifinit)
      pot = 0.0d0
      ptfrc(1) = 0.0d0
      ptfrc(2) = 0.0d0
      ptfrc(3) = 0.0d0
c  
      do jpatch = 1, ntri
         call tri_ini(triang(1,1,jpatch),triang(1,2,jpatch),
     1                triang(1,3,jpatch),w,vert1,vert2,vert3)
         call triasymq(nquad,vert1,vert2,vert3,rnodes,weights,
     1                nnodes)
c
         call tri_for(w,targ,vertout)
         x0 = vertout(1)
         y0 = vertout(2)
         z0 = vertout(3)
         call triquadhelmd2(ifl,vert1,vert2,vert3,x0,y0,z0,zk,
     1      ifpot,iffld,zval,zvalx,zvaly,zvalz,
     $      nnodes,weights,rnodes,ier)
c
c       rotate the derivative vector, normal derivative is (0,0,1)
c       in the rotated coordinate system
c
         if( ifpot .eq. 1 ) then
         pot = pot - dipstr(jpatch)*zval
         endif
c
         if( iffld .eq. 1 ) then
         call rotder3dz(w,
     $      zvalx,zvaly,zvalz,zderx,zdery,zderz)
c
         ptfrc(1) = ptfrc(1) + dipstr(jpatch)*zderx
         ptfrc(2) = ptfrc(2) + dipstr(jpatch)*zdery
         ptfrc(3) = ptfrc(3) + dipstr(jpatch)*zderz
         endif
c
      enddo
      return
      end
c
c
c
c
c
c***********************************************************************
c
c                  Symmetric Quadrature routines
c
c***********************************************************************
      subroutine triquadhelm(ifl,vert1,vert2,vert3,x0,y0,z0,zk,
     1           zval,zvalx,zvaly,zvalz,nquad,wts,rnodes,ier)
c***********************************************************************
        implicit real *8 (a-h,o-z)
        dimension vert1(2),vert2(2),vert3(2)
        dimension wts(1),rnodes(2,1)
        complex *16 zk, zval,zvalx,zvaly,zvalz
        complex *16 zsumacross,zsumacrossx,zsumacrossy,zsumacrossz
        complex *16 gdiff,gdiffx,gdiffy,gdiffz
        complex *16 zsum,zsumx,zsumy,zsumz,zarg,zz,zexp,eye
c
c       If ifl = 0 then
c
c       this subroutine evaluates the integral
c       
c       \int_T exp(i zk*r)/r                        (A)
c
c       and its derivatives
c       
c       where r = sqrt((x-x0)^2+(y-y0)^2+z0^2)
c
c       over the triangle T with the vertices (vert1,vert2,vert3)
c
c       If ifl = 1 then it evaluates the integral of the difference
c       between the Laplace and the Helmholtz kernel
c
c       \int_T  [exp(i zk*r)/r - 1/r]               (B)
c
c
c       Algorithm:  
c       Analytic integration for Laplace kernel (if ifl=0),
c       followed by Gaussian quadrature on (B).
c
c     INPUT:
c
c     vert1, vert2, vert3 are vertices which have been mapped
c         to the x,y plane with vert1 at the origin, vert2 on the 
c         positive x axis and vert3 in the upper half plane.
c
c     x0,y0,z0  are coordinates of the mapped target point
c     nquad     order of Gaussian rule 
c     wts       standard Gaussian weights using nquad pts.
c     xnodes    standard Gaussian nodes   using nquad pts.
c
c     OUTPUT:
c
c     zval is value of potential at x0,y0,z0
c     zvalx is value of x-derivative of potential at x0,y0,z0
c     zvaly is value of y-derivative of potential at x0,y0,z0
c     zvalz is value of z-derivative of potential at x0,y0,z0
c     ier  is error return code:
c          0 = normal execution
c          1 = vert3 is in lower half plane => not positively oriented.
c
c------------------------------------------------------------------
c
	eye = dcmplx(0.0d0,1.0d0)
        a=x0
        b=y0
        c=z0
c
c
      ier = 0
      if (vert3(2).le.0) then
c
c     triangle is not positively oriented - return with error code 1.
c
         ier = 1
	 return
      endif
c
c    Get Laplace quadrature if ifl is zero.
c
      if (ifl.eq.0) then
         iquad = 0
         if (z0.gt.0) iquad = 1
         if (z0.lt.0) iquad = -1
         itype = 1
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      val)
         zsum = val
         itype = 2
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      val)
         zsumx = val
         itype = 3
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      val)
         zsumy = val
         itype = 4
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      val)
         zsumz = -val
      else
         zsum = 0.0d0
         zsumx = 0.0d0
         zsumy = 0.0d0
         zsumz = 0.0d0
      endif
ccc      call prin2(' zsum is *',zsum,1)
c
c     now compute integrals over the triangle of 
c     [ exp(i*zk*r)/r - 1/r] (and derivatives) and add to 
c     the preceding results.
c
c     In the inner loop, when abs(zk*rr) is small, the exponential
c     is approximated by a four-term Taylor series and handled
c     analytically from there. This avoids catastrophic cancellation.
c
      do i = 1,nquad
	 x = rnodes(1,i)
	 y = rnodes(2,i)
	 dx = (x0-x)
	 dy = (y0-y)
	 dz = z0
	 rr2 = dx**2 + dy**2 + dz**2
	 rr = dsqrt(rr2)
	 rr3 = rr*rr2
	 zarg = eye*zk*rr
	 zz = eye*zk
	 if (abs(zarg) .lt. 0.01d0) then
            gdiff = zz*(1.0d0+zarg/2.0d0+zarg**2/6.0d0+
     1                 zarg**3/24.0d0 + zarg**4/120.0d0)
            gdiffx = (dx/rr)*( zz**2/2.0d0 + rr*zz**3/3.0d0 + 
     1           	rr2*zz**4/8.0d0 + rr3*zz**5/30.0d0)
            gdiffy = (dy/rr)*( zz**2/2.0d0 + rr*zz**3/3.0d0 + 
     1           	rr2*zz**4/8.0d0 + rr3*zz**5/30.0d0)
            gdiffz = (dz/rr)*( zz**2/2.0d0 + rr*zz**3/3.0d0 + 
     1           	rr2*zz**4/8.0d0 + rr3*zz**5/30.0d0)
	 else
	    zexp = exp(zarg)
	    gdiff = (zexp - 1.0d0)/rr
	    gdiffx = dx*(zexp*(zarg- 1.0d0) + 1.0d0)/rr3
	    gdiffy = dy*(zexp*(zarg- 1.0d0) + 1.0d0)/rr3
	    gdiffz = dz*(zexp*(zarg- 1.0d0) + 1.0d0)/rr3
         endif
	 ww = wts(i)
	 zsum = zsum + gdiff*ww
	 zsumx = zsumx + gdiffx*ww
	 zsumy = zsumy + gdiffy*ww
	 zsumz = zsumz + gdiffz*ww
      enddo
ccc      call prin2(' zsum is *',zsum,1)
ccc      call prin2(' zsumx is *',zsumx,1)
ccc      call prin2(' zsumy is *',zsumy,1)
ccc      call prin2(' zsumz is *',zsumz,1)
      zval = zsum
      zvalx = zsumx
      zvaly = zsumy
      zvalz = zsumz
      return
      end
c
c
c
c
c
      subroutine triquadhelm2(ifl,vert1,vert2,vert3,x0,y0,z0,zk,
     1     ifpot,iffld,zval,zvalx,zvaly,zvalz,
     $     nquad,wts,rnodes,ier)
c***********************************************************************
        implicit real *8 (a-h,o-z)
        dimension vert1(2),vert2(2),vert3(2)
        dimension wts(1),rnodes(2,1)
        complex *16 zk, zval,zvalx,zvaly,zvalz
        complex *16 zsumacross,zsumacrossx,zsumacrossy,zsumacrossz
        complex *16 gdiff,gdiffx,gdiffy,gdiffz
        complex *16 zsum,zsumx,zsumy,zsumz,zarg,zz,zexp,eye
c
c       If ifl = 0 then
c
c       this subroutine evaluates the integral
c       
c       \int_T exp(i zk*r)/r                        (A)
c
c       and its derivatives
c       
c       where r = sqrt((x-x0)^2+(y-y0)^2+z0^2)
c
c       over the triangle T with the vertices (vert1,vert2,vert3)
c
c       If ifl = 1 then it evaluates the integral of the difference
c       between the Laplace and the Helmholtz kernel
c
c       \int_T  [exp(i zk*r)/r - 1/r]               (B)
c
c
c       Algorithm:  
c       Analytic integration for Laplace kernel (if ifl=0),
c       followed by Gaussian quadrature on (B).
c
c     INPUT:
c
c     vert1, vert2, vert3 are vertices which have been mapped
c         to the x,y plane with vert1 at the origin, vert2 on the 
c         positive x axis and vert3 in the upper half plane.
c
c     x0,y0,z0  are coordinates of the mapped target point
c     nquad     order of Gaussian rule 
c     wts       standard Gaussian weights using nquad pts.
c     xnodes    standard Gaussian nodes   using nquad pts.
c
c     OUTPUT:
c
c     zval is value of potential at x0,y0,z0
c     zvalx is value of x-derivative of potential at x0,y0,z0
c     zvaly is value of y-derivative of potential at x0,y0,z0
c     zvalz is value of z-derivative of potential at x0,y0,z0
c     ier  is error return code:
c          0 = normal execution
c          1 = vert3 is in lower half plane => not positively oriented.
c
c------------------------------------------------------------------
c
	eye = dcmplx(0.0d0,1.0d0)
        a=x0
        b=y0
        c=z0
c
c
      ier = 0
      if (vert3(2).le.0) then
c
c     triangle is not positively oriented - return with error code 1.
c
         ier = 1
	 return
      endif
c
         zsum = 0.0d0
         zsumx = 0.0d0
         zsumy = 0.0d0
         zsumz = 0.0d0
c
c    Get Laplace quadrature if ifl is zero.
c
      if (ifl.eq.0) then
         iquad = 0
         if (z0.gt.0) iquad = 1
         if (z0.lt.0) iquad = -1
         if( ifpot .eq. 1 ) then
         itype = 1
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      val)
         zsum = val
         endif
         if( iffld .eq. 1 ) then
         itype = 2
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      val)
         zsumx = val
         itype = 3
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      val)
         zsumy = val
         itype = 4
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1      val)
         zsumz = -val
         endif
      else
         zsum = 0.0d0
         zsumx = 0.0d0
         zsumy = 0.0d0
         zsumz = 0.0d0
      endif
ccc      call prin2(' zsum is *',zsum,1)
c
c     now compute integrals over the triangle of 
c     [ exp(i*zk*r)/r - 1/r] (and derivatives) and add to 
c     the preceding results.
c
c     In the inner loop, when abs(zk*rr) is small, the exponential
c     is approximated by a four-term Taylor series and handled
c     analytically from there. This avoids catastrophic cancellation.
c
      do i = 1,nquad
	 x = rnodes(1,i)
	 y = rnodes(2,i)
	 dx = (x0-x)
	 dy = (y0-y)
	 dz = z0
	 rr2 = dx**2 + dy**2 + dz**2
	 rr = dsqrt(rr2)
	 rr3 = rr*rr2
	 zarg = eye*zk*rr
	 zz = eye*zk
	 if (abs(zarg) .lt. 0.01d0) then
            if( ifpot .eq. 1 ) then
            gdiff = zz*(1.0d0+zarg/2.0d0+zarg**2/6.0d0+
     1                 zarg**3/24.0d0 + zarg**4/120.0d0)
            endif
            if( iffld .eq. 1 ) then
            gdiffx = (dx/rr)*( zz**2/2.0d0 + rr*zz**3/3.0d0 + 
     1           	rr2*zz**4/8.0d0 + rr3*zz**5/30.0d0)
            gdiffy = (dy/rr)*( zz**2/2.0d0 + rr*zz**3/3.0d0 + 
     1           	rr2*zz**4/8.0d0 + rr3*zz**5/30.0d0)
            gdiffz = (dz/rr)*( zz**2/2.0d0 + rr*zz**3/3.0d0 + 
     1           	rr2*zz**4/8.0d0 + rr3*zz**5/30.0d0)
            endif
	 else
	    zexp = exp(zarg)
            if( ifpot .eq. 1 ) then
	    gdiff = (zexp - 1.0d0)/rr
            endif
            if( iffld .eq. 1 ) then
	    gdiffx = dx*(zexp*(zarg- 1.0d0) + 1.0d0)/rr3
	    gdiffy = dy*(zexp*(zarg- 1.0d0) + 1.0d0)/rr3
	    gdiffz = dz*(zexp*(zarg- 1.0d0) + 1.0d0)/rr3
            endif
         endif
	 ww = wts(i)
            if( ifpot .eq. 1 ) then
            zsum = zsum + gdiff*ww
            endif
            if( iffld .eq. 1 ) then
            zsumx = zsumx + gdiffx*ww
            zsumy = zsumy + gdiffy*ww
            zsumz = zsumz + gdiffz*ww
            endif
      enddo
ccc      call prin2(' zsum is *',zsum,1)
ccc      call prin2(' zsumx is *',zsumx,1)
ccc      call prin2(' zsumy is *',zsumy,1)
ccc      call prin2(' zsumz is *',zsumz,1)
c
        if( ifpot .eq. 1 ) then
        zval = zsum
        else
        zval = 0
        endif
        if( iffld .eq. 1 ) then
        zvalx = zsumx
        zvaly = zsumy
        zvalz = zsumz
        else
        zvalx = 0
        zvaly = 0
        zvalz = 0
        endif
c
      return
      end
c
c
c
c
c
c***********************************************************************
      subroutine triquadhelmd(ifl,vert1,vert2,vert3,x0,y0,z0,zk,
     1                   zval,zvalx,zvaly,zvalz,nquad,
     2                   wts,rnodes,ier)
c***********************************************************************
      implicit real *8 (a-h,o-z)
      dimension vert1(2),vert2(2),vert3(2)
      dimension wts(1),rnodes(2,1)
      complex *16 zk, zval,zvalx,zvaly,zvalz
      complex *16 zsumacross,zsumacrossx,zsumacrossy,zsumacrossz
      complex *16 zsacross,zsacrossx,zsacrossy,zsacrossz
      complex *16 gdiff,gdiffx,gdiffy,gdiffz
      complex *16 gd,gdx,gdy,gdz
      complex *16 zsum,zsumx,zsumy,zsumz,zarg,zz,zexp,eye
      complex *16 zs,zsx,zsy,zsz
      complex *16 ztaylor1
      complex *16 ztaylor2
c
c
c     If ifl = 0 then
c
c       this subroutine evaluates the double layer potential D_k on
c       a flat triangle as well as its gradient.
c       
c       \int_T d/dz  [ exp(i zk*r)/r ]                  (A)
c
c     where r = sqrt((x-x0)^2+(y-y0)^2+z0^2)
c     over the triangle T with the vertices (vert1,vert2,vert3)
c
c     If ifl = 1 then it evaluates the integral of the difference
c       between the Laplace and the Helmholtz kernel
c
c       \int_T d/dz  [exp(i zk*r)/r - 1/r]               (B)
c
c
c     Method: Analytic integration for Laplace kernel (if ifl=0),
c     followed by Gaussian quadrature on (B).
c
c     INPUT:
c
c     vert1, vert2, vert3 are vertices which have been mapped
c         to the x,y plane with vert1 at the origin, vert2 on the 
c         positive x axis and vert3 in the upper half plane.
c
c     x0,y0,z0  are coordinates of the mapped target point
c     nquad     order of Gaussian rule to be used in y-direction.  
c     wts       standard Gaussian weights using nquad pts.
c     xnodes    standard Gaussian nodes   using nquad pts.
c
c     OUTPUT:
c
c     zval is value of potential at x0,y0,z0
c     zvalx is value of x-derivative of potential at x0,y0,z0
c     zvaly is value of y-derivative of potential at x0,y0,z0
c     zvalz is value of z-derivative of potential at x0,y0,z0
c     ier  is error return code:
c          0 = normal execution
c          1 = vert3 is in lower half plane => not positively oriented.
c
c------------------------------------------------------------------
c
c   see quadnotes.txt for derivative/Taylor series calculations
c
	eye = dcmplx(0.0d0,1.0d0)
        a=x0
        b=y0
        c=z0
c
c
      ier = 0
      if (vert3(2).le.0) then
c
c     triangle is not positively oriented - return with error code 1.
c
         ier = 1
	 return
      endif
c
c    Carry out anaylitc quadrature for Laplace case if ifl equals 0.
c
      zs = 0.0d0
      zsx = 0.0d0
      zsy = 0.0d0
      zsz = 0.0d0
      if (ifl.eq.0) then
         iquad = 0
         if (z0.gt.0) iquad = +1
         if (z0.lt.0) iquad = -1
         itype = 4
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
         val = -val
         itype = 5
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
         itype = 6
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
         itype = 7
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valz)
      else
         val = 0.0d0
         valx = 0.0d0
         valy = 0.0d0
         valz = 0.0d0
      endif
c
c     now compute integrals over the triangle of 
c     d/dz [ exp(i*zk*r)/r - 1/r] (and derivatives) and add to 
c     the preceding results.
c
c     In the inner loop, when abs(zk*rr) is small, the exponential
c     is approximated by a four-term Taylor series and handled
c     analytically from there. This avoids catastrophic cancellation.
c
      do i = 1,nquad
	 x = rnodes(1,i)
	 y = rnodes(2,i)
	 dx = (x0-x)
	 dy = (y0-y)
	 dz = z0
	 rr2 = dx**2 + dy**2 + dz**2
	 rr = dsqrt(rr2)
	 rr3 = rr*rr2
	 zarg = eye*zk*rr
	 zz = eye*zk
	 zexp = exp(zarg)
	 if (abs(zarg) .lt. 0.01d0) then
            ztaylor1 = ( zz**2/2.0d0 + rr*zz**3/3.0d0 + 
     1           	rr2*zz**4/8.0d0 + rr3*zz**5/30.0d0)
            gd = (dz/rr)*ztaylor1
            gdx = (dx*dz/rr3)*(zz*zz*zexp - 3*ztaylor1)
            gdy = (dy*dz/rr3)*(zz*zz*zexp - 3*ztaylor1)
            gdz = ztaylor1/rr-(dz*dz/rr3)*(zexp*zk*zk+3*ztaylor1)
	 else
	    gd = (dz/rr)*(zexp*(zarg- 1.0d0) + 1.0d0)/rr2
            gdx = (dx*dz/rr3)*(zz*zz*zexp - 
     1               3*(zexp*(zarg-1.0d0)+1.0d0)/rr2)
            gdy = (dy*dz/rr3)*(zz*zz*zexp - 
     1               3*(zexp*(zarg-1.0d0)+1.0d0)/rr2)
            gdz = (zexp*(zarg-1.0d0)+1.0d0)/rr3
            gdz = gdz-zexp*(zk*zk*dz*dz)/rr3
            gdz = gdz-3*(dz*dz/rr3)*(zexp*(zarg-1.0d0)+1.0d0)/rr2
         endif
	 ww = wts(i)
	 zs = zs + gd*ww
	 zsx = zsx + gdx*ww
	 zsy = zsy + gdy*ww
	 zsz = zsz + gdz*ww
      enddo
      zval = val+zs
      zvalx = -valx+zsx
      zvaly = -valy+zsy
      zvalz = valz+zsz
      return
      end
c
c
c
c
c
c***********************************************************************
      subroutine triquadhelmd2(ifl,vert1,vert2,vert3,x0,y0,z0,zk,
     1                   ifpot,iffld,zval,zvalx,zvaly,zvalz,nquad,
     2                   wts,rnodes,ier)
c***********************************************************************
      implicit real *8 (a-h,o-z)
      dimension vert1(2),vert2(2),vert3(2)
      dimension wts(1),rnodes(2,1)
      complex *16 zk, zval,zvalx,zvaly,zvalz
      complex *16 zsumacross,zsumacrossx,zsumacrossy,zsumacrossz
      complex *16 zsacross,zsacrossx,zsacrossy,zsacrossz
      complex *16 gdiff,gdiffx,gdiffy,gdiffz
      complex *16 gd,gdx,gdy,gdz
      complex *16 zsum,zsumx,zsumy,zsumz,zarg,zz,zexp,eye
      complex *16 zs,zsx,zsy,zsz
      complex *16 ztaylor1
      complex *16 ztaylor2
c
c
c     If ifl = 0 then
c
c       this subroutine evaluates the double layer potential D_k on
c       a flat triangle as well as its gradient.
c       
c       \int_T d/dz  [ exp(i zk*r)/r ]                  (A)
c
c     where r = sqrt((x-x0)^2+(y-y0)^2+z0^2)
c     over the triangle T with the vertices (vert1,vert2,vert3)
c
c     If ifl = 1 then it evaluates the integral of the difference
c       between the Laplace and the Helmholtz kernel
c
c       \int_T d/dz  [exp(i zk*r)/r - 1/r]               (B)
c
c
c     Method: Analytic integration for Laplace kernel (if ifl=0),
c     followed by Gaussian quadrature on (B).
c
c     INPUT:
c
c     vert1, vert2, vert3 are vertices which have been mapped
c         to the x,y plane with vert1 at the origin, vert2 on the 
c         positive x axis and vert3 in the upper half plane.
c
c     x0,y0,z0  are coordinates of the mapped target point
c     nquad     order of Gaussian rule to be used in y-direction.  
c     wts       standard Gaussian weights using nquad pts.
c     xnodes    standard Gaussian nodes   using nquad pts.
c
c     OUTPUT:
c
c     zval is value of potential at x0,y0,z0
c     zvalx is value of x-derivative of potential at x0,y0,z0
c     zvaly is value of y-derivative of potential at x0,y0,z0
c     zvalz is value of z-derivative of potential at x0,y0,z0
c     ier  is error return code:
c          0 = normal execution
c          1 = vert3 is in lower half plane => not positively oriented.
c
c------------------------------------------------------------------
c
c   see quadnotes.txt for derivative/Taylor series calculations
c
	eye = dcmplx(0.0d0,1.0d0)
        a=x0
        b=y0
        c=z0
c
c
      ier = 0
      if (vert3(2).le.0) then
c
c     triangle is not positively oriented - return with error code 1.
c
         ier = 1
	 return
      endif
c
c    Carry out anaylitc quadrature for Laplace case if ifl equals 0.
c
      zs = 0.0d0
      zsx = 0.0d0
      zsy = 0.0d0
      zsz = 0.0d0
      if (ifl.eq.0) then
         iquad = 0
         if (z0.gt.0) iquad = +1
         if (z0.lt.0) iquad = -1
         if( ifpot .eq. 1 ) then
         itype = 4
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
         val = -val
         endif
         if( iffld .eq. 1 ) then
         itype = 5
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
         itype = 6
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
         itype = 7
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valz)
         endif
      else
         val = 0.0d0
         valx = 0.0d0
         valy = 0.0d0
         valz = 0.0d0
      endif
c
c     now compute integrals over the triangle of 
c     d/dz [ exp(i*zk*r)/r - 1/r] (and derivatives) and add to 
c     the preceding results.
c
c     In the inner loop, when abs(zk*rr) is small, the exponential
c     is approximated by a four-term Taylor series and handled
c     analytically from there. This avoids catastrophic cancellation.
c
      do i = 1,nquad
	 x = rnodes(1,i)
	 y = rnodes(2,i)
	 dx = (x0-x)
	 dy = (y0-y)
	 dz = z0
	 rr2 = dx**2 + dy**2 + dz**2
	 rr = dsqrt(rr2)
	 rr3 = rr*rr2
	 zarg = eye*zk*rr
	 zz = eye*zk
	 zexp = exp(zarg)
	 if (abs(zarg) .lt. 0.01d0) then
            ztaylor1 = ( zz**2/2.0d0 + rr*zz**3/3.0d0 + 
     1           	rr2*zz**4/8.0d0 + rr3*zz**5/30.0d0)
            if( ifpot .eq. 1 ) then
            gd = (dz/rr)*ztaylor1
            endif
            if( iffld .eq. 1 ) then
            gdx = (dx*dz/rr3)*(zz*zz*zexp - 3*ztaylor1)
            gdy = (dy*dz/rr3)*(zz*zz*zexp - 3*ztaylor1)
            gdz = ztaylor1/rr-(dz*dz/rr3)*(zexp*zk*zk+3*ztaylor1)
            endif
	 else
            if( ifpot .eq. 1 ) then
	    gd = (dz/rr)*(zexp*(zarg- 1.0d0) + 1.0d0)/rr2
            endif
            if( iffld .eq. 1 ) then
            gdx = (dx*dz/rr3)*(zz*zz*zexp - 
     1               3*(zexp*(zarg-1.0d0)+1.0d0)/rr2)
            gdy = (dy*dz/rr3)*(zz*zz*zexp - 
     1               3*(zexp*(zarg-1.0d0)+1.0d0)/rr2)
            gdz = (zexp*(zarg-1.0d0)+1.0d0)/rr3
            gdz = gdz-zexp*(zk*zk*dz*dz)/rr3            
            gdz = gdz-3*(dz*dz/rr3)*(zexp*(zarg-1.0d0)+1.0d0)/rr2
            endif
         endif
	 ww = wts(i)
         if( ifpot .eq. 1 ) then
	 zs = zs + gd*ww
         endif
         if( iffld .eq. 1 ) then
	 zsx = zsx + gdx*ww
	 zsy = zsy + gdy*ww
	 zsz = zsz + gdz*ww
         endif
      enddo
c
        if( ifpot .eq. 1 ) then
        zval = val+zs
        else
        zval = 0
        endif
        if( iffld .eq. 1 ) then
        zvalx = -valx+zsx
        zvaly = -valy+zsy
        zvalz = valz+zsz
        else
        zvalx = 0
        zvaly = 0
        zvalz = 0
        endif
c
      return
      end
c
c
c
c
c
c***********************************************************************
      subroutine triquadselfhelm(ifl,vert1,vert2,vert3,x0,y0,zk,
     1                   zval,zvalx,zvaly,zvalz,nquad,
     2                   wts,rnodes,ier)
c***********************************************************************
        implicit real *8 (a-h,o-z)
        dimension vert1(2),vert2(2),vert3(2)
        dimension wts(1),rnodes(2,1)
        complex *16 zk, zval,zvalx,zvaly,zvalz
        complex *16 zsumacross,zsumacrossx,zsumacrossy,zsumacrossz
        complex *16 gdiff,gdiffx,gdiffy,gdiffz
        complex *16 zsum,zsumx,zsumy,zsumz,zarg,zz,zexp,eye
c
c
c       If (ifl .eq. 0) then this subroutine evaluates the integral
c       
c       \int_T exp(i zk*r)/r 
c
c       where r = sqrt((x-x0)^2+(y-y0)^2+z0^2)
c
c       over the triangle T with the vertices (vert1,vert2,vert3)
c
c       If (ifl .neq. 0) then this subroutine evaluates the integral
c       
c       \int_T  [ exp(i zk*r)/r - 1/r ]                  (B)
c
c       Algorithm:  Analytic integration for Laplace kernel 
c       (if ifl.eq.0), followed by Gaussian quadrature for (B).
c
c     INPUT:
c
c     vert1, vert2, vert3 are vertices which have been mapped
c         to the x,y plane with vert1 at the origin, vert2 on the 
c         positive x axis and vert3 in the upper half plane.
c
c     x0,y0,z0  are coordinates of the mapped target point
c     nquad     order of Gaussian rule to be used in y-direction.  
c     wts       standard Gaussian weights using nquad pts.
c     xnodes    standard Gaussian nodes   using nquad pts.
c
c     OUTPUT:
c
c     zval is value of potential at x0,y0,z0
c     zvalx is value of x-derivative of potential at x0,y0,z0
c     zvaly is value of y-derivative of potential at x0,y0,z0
c     zvalz is value of z-derivative of potential at x0,y0,z0
c     ier  is error return code:
c          0 = normal execution
c          1 = vert3 is in lower half plane => not positively oriented.
c
c------------------------------------------------------------------
c
      eye = dcmplx(0.0d0,1.0d0)
      a=x0
      b=y0
      z0=0
      c=z0
c
ccc      call prinf(' nquad is *',nquad,1)
c
      ier = 0
      if (vert3(2).le.0) then
c
c     triangle is not positively oriented - return with error code 1.
c
         ier = 1
	 return
      endif
c
c    Laplace quadrature if fil equals zero.
c
      if (ifl.eq.0) then
         itype = 1
         iquad = 0
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
         zsum = val
         itype = 2
         iquad = 0
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
         zsumx = valx
         itype = 3
         iquad = 0
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
         zsumy = valy
ccc         itype = 4
ccc         iquad = 0
ccc         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
ccc     1   	   valz)
ccc         zsumz = -valz
         zsumz = 0.0d0
      else
         zsum = 0.0d0
         zsumx = 0.0d0
         zsumy = 0.0d0
         zsumz = 0.0d0
      endif
c
c     now compute integrals over the triangle of 
c     [ exp(i*zk*r)/r - 1/r] (and derivatives) and add to 
c     the preceding results.
c
      itype = 8
      iquad = 0
      call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valx)
      valx = -valx
      zsumx = zsumx - (zk**2)*valx/2.0d0
      itype = 9
      iquad = 0
      call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valy)
      valy = -valy
      zsumy = zsumy - (zk**2)*valy/2.0d0
c
c     In the inner loop, when abs(zk*rr) is small, the exponential
c     is approximated by a four-term Taylor series and handled
c     analytically from there. This avoids catastrophic cancellation.
c
      do i = 1,nquad
	 x = rnodes(1,i)
	 y = rnodes(2,i)
	 dx = (x0-x)
	 dy = (y0-y)
	 dz = z0
	 rr2 = dx**2 + dy**2 + dz**2
	 rr = dsqrt(rr2)
	 rr3 = rr*rr2
	 zarg = eye*zk*rr
	 zz = eye*zk
	 if (abs(zarg) .lt. 0.01d0) then
            gdiff = zz*(1.0d0+zarg/2.0d0+zarg**2/6.0d0+
     1                 zarg**3/24.0d0 + zarg**4/120.0d0)
ccc            gdiffx = (dx/rr)*( zz**2/2.0d0 + rr*zz**3/3.0d0 + 
ccc     1           	rr2*zz**4/8.0d0 + rr3*zz**5/30.0d0)
            gdiffx = dx*(zz**3/3.0d0 + 
     1           	rr*zz**4/8.0d0 + rr2*zz**5/30.0d0)
ccc            gdiffy = (dy/rr)*( zz**2/2.0d0 + rr*zz**3/3.0d0 + 
ccc     1           	rr2*zz**4/8.0d0 + rr3*zz**5/30.0d0)
            gdiffy = dy*(zz**3/3.0d0 + 
     1           	rr*zz**4/8.0d0 + rr2*zz**5/30.0d0)
ccc            gdiffz = (dz/rr)*( zz**2/2.0d0 + rr*zz**3/3.0d0 + 
ccc     1           	rr2*zz**4/8.0d0 + rr3*zz**5/30.0d0)
	 else
	    zexp = exp(zarg)
	    gdiff = (zexp - 1.0d0)/rr
	    gdiffx = dx*(zexp*(zarg- 1.0d0) + 1.0d0)/rr3
	    gdiffx = gdiffx - (dx/rr)*(zz**2)/2.0d0
	    gdiffy = dy*(zexp*(zarg- 1.0d0) + 1.0d0)/rr3
	    gdiffy = gdiffy - (dy/rr)*(zz**2)/2.0d0
ccc	    gdiffz = dz*(zexp*(zarg- 1.0d0) + 1.0d0)/rr3
         endif
	 ww = wts(i)
	 zsum = zsum + gdiff*ww
	 zsumx = zsumx + gdiffx*ww
	 zsumy = zsumy + gdiffy*ww
ccc	 zsumz = zsumz + gdiffz*ww
ccc         write(11,*) x,y,dble(gdiffx)
      enddo
ccc      stop
      zval = zsum
      zvalx = zsumx
      zvaly = zsumy
ccc      zvalz = zsumz
      zvalz = 0.0d0
      return
      end
c
c
c
c
c
c
c
c***********************************************************************
      subroutine triquadselfhelmd(ifl,vert1,vert2,vert3,x0,y0,zk,
     1                   zval,zvalx,zvaly,zvalz,nquad,
     2                   wts,rnodes,ier)
c***********************************************************************
        implicit real *8 (a-h,o-z)
        dimension vert1(2),vert2(2),vert3(2)
        dimension wts(1),rnodes(2,1)
        complex *16 zk, zval,zvalx,zvaly,zvalz
        complex *16 zsumacross,zsumacrossx,zsumacrossy,zsumacrossz
        complex *16 gdiff,gdiffx,gdiffy,gdiffz
        complex *16 zsum,zsumx,zsumy,zsumz,zarg,zz,zexp,eye
        complex *16 zt1,zt2
c
c
c     If (ifl.eq.0) then this subroutine evaluates the double layer 
c     potential D_k on a flat triangle as well as its gradient.
c       
c       \int_T d/dz  [ exp(i zk*r)/r ]
c
c     where r = sqrt((x-x0)^2+(y-y0)^2+z0^2)
c
c    over the triangle T with the vertices (vert1,vert2,vert3)
c
c     If (ifl.ne.0) this subroutine evaluates the difference of
c      double layer potentials (Laplace - Helmholtz).
c       
c    In fact the potential is zero as are the 
c    x and y derivatives, but the z derivative
c    is not. If (ifl.eq.0), we carry out analytic integration
c    for the Laplace kernal, followed by Gaussian quadrature on 
c    difference kernel (after removing leading order term
c    which is ( (1/2)zk^2/r). )
c
c     INPUT:
c
c     vert1, vert2, vert3 are vertices which have been mapped
c         to the x,y plane with vert1 at the origin, vert2 on the 
c         positive x axis and vert3 in the upper half plane.
c
c     x0,y0,z0  are coordinates of the mapped target point
c     nquad     order of Gaussian rule to be used in y-direction.  
c     wts       standard Gaussian weights using nquad pts.
c     xnodes    standard Gaussian nodes   using nquad pts.
c
c     OUTPUT:
c
c     zval is value of potential at x0,y0,z0
c     zvalx is value of x-derivative of potential at x0,y0,z0
c     zvaly is value of y-derivative of potential at x0,y0,z0
c     zvalz is value of z-derivative of potential at x0,y0,z0
c     ier  is error return code:
c          0 = normal execution
c          1 = vert3 is in lower half plane => not positively oriented.
c
c------------------------------------------------------------------
c
c   see quadnotes.txt for derivative/Taylor series calculations
c
      eye = dcmplx(0.0d0,1.0d0)
      a=x0
      b=y0
      z0 = 0
      c=z0
c
      if (ifl.eq.0) then
         itype = 7
         iquad = 0
         call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   valz)
         zsumz = valz
      else
         zsumz = 0.0d0
      endif
      itype = 1
      iquad = 0
      call triahquad(itype,iquad,vert1,vert2,vert3,x0,y0,z0,
     1   	   val)
      zsumz = zsumz - (zk**2)*val/2.0d0
c
c     now compute integral over the triangle of 
c     d/dz d/dz [ exp(i*zk*r)/r - 1/r] + (1/2)*zk**2/r and add to 
c     the preceding results.
c
c     In the inner loop, when abs(zk*rr) is small, the exponential
c     is approximated by a four-term Taylor series and handled
c     analytically from there. This avoids catastrophic cancellation.
c
      do i = 1,nquad
	 x = rnodes(1,i)
	 y = rnodes(2,i)
	 dx = (x0-x)
	 dy = (y0-y)
	 dz = z0
	 rr2 = dx**2 + dy**2 + dz**2
	 rr = dsqrt(rr2)
	 rr3 = rr*rr2
	 zarg = eye*zk*rr
	 zz = eye*zk
	 zexp = exp(zarg)
         if (abs(zarg).lt.0.01d0) then
            gdiffz = zz**3/3.0d0 + rr*zz**4/8.0d0 + 
     1           	rr2*zz**5/30.0d0 + rr3*zz**6/144.0d0 
         else
            gdiffz = (zexp*(zarg-1.0d0)+1.0d0)/rr3
            gdiffz = gdiffz + (zk**2)/(2.0d0*rr)
         endif
	 ww = wts(i)
	 zsumz = zsumz + gdiffz*ww
      enddo
      zval = 0.0d0
      zvalx = 0.0d0
      zvaly = 0.0d0
      zvalz = zsumz
      return
      end
c
c
c
c
c
c***********************************************************************
      subroutine triquadselfhelmold(v1,v2,v3,x0,y0,zk,
     1                   zval,zvalx,zvaly,zvalz,nquad,wts,xnodes,ier)
c***********************************************************************
      implicit real *8 (a-h,o-z)
      dimension v1(2),v2(2),v3(2)
      dimension wts(1),xnodes(1)
      complex *16 zval,zvalx,zvaly,zvalz,zk,zarg,zexp,eye
      complex *16 sum,zz
c
c
c       this subroutine evaluates the integral
c       
c       \int_T exp(-rlambda*r)/r
c
c       where r = sqrt((x-x0)^2+(y-y0)^2)
c
c       over the triangle T with the vertices (v1,v2,v3)
c
c       at the point P(x0,y0) lying in T (not necessarily at the
c       centroid).
c       The integral is computed in three pieces 
c
c
c                            v3
c                           / |\
c                          /  | \
c                         /   |  \
c                        /    |   \ 
c                       / III | I  \ 
c                      /      P     \ 
c                     /  _ /    \ _   \ 
c                    / /    II      \  \
c                  v1------------------v2
c
c       with the coordinate system rotated differently each time 
c       with P at the "origin" so that the 
c       calculation is easy in polar coordinates 
c
c       \int_0^\theta  \int_0^r(phi)  1/r r dr dphi
c       = \int_0^\theta  r(phi) dphi
c
c     When handling triangle I, for example, P-v2 is considered
c     phi = 0 and and phi sweeps counterclockwise until P-v3 is 
c     reached at the maximal angular sweep (denoted theta below).
c     The analogous thing is done for II and III.
c
c
c     INPUT:
c
c     v1, v2, v3 are vertices which have been mapped
c         to the x,y plane 
c     x0,y0,z0  are coordinates of the mapped interior target point
c     nquad     order of Gaussian rule to be used in theta-direction
c     wts       standard Gaussian weights using nquad pts.
c     xnodes    standard Gaussian nodes   using nquad pts.
c
c     OUTPUT:
c
c     val is value of potential at x0,y0
c     valx is value of x-derivative of potential at x0,y0
c          NOT YET IMPLEMENTED - NOT NEEDED FOR MOST BOUNDARY
c          VALUE PROBLEMS - code returns 0.
c     valy is value of y-derivative of potential at x0,y0
c          NOT YET IMPLEMENTED - NOT NEEDED FOR MOST BOUNDARY
c          VALUE PROBLEMS - code returns 0.
c     valz is value of z-derivative of potential at x0,y0 which 
c             is zero.
c     ier  no error return code at present
c          0 = normal execution
c
c------------------------------------------------------------------
c
c
      eye = dcmplx(0.0d0,1.0d0)
      ier = 0
c
c     process subtriangle I: (P-V2-V3)
c
      a11 = v3(1)-v2(1)
      a21 = v3(2)-v2(2)
      e11 = v2(1)-x0
      e12 = v2(2)-y0
      dpv2 = dsqrt(e11**2 +e12**2)
      dpv1 = dsqrt( (v1(1)-x0)**2 + (v1(2)-y0)**2)
      dpv3 = dsqrt( (v3(1)-x0)**2 + (v3(2)-y0)**2)
      e11 = e11/dpv2
      e12 = e12/dpv2
      e21 = -e12
      e22 =  e11
      rhs1 = x0-v2(1)
      rhs2 = y0-v2(2)
      ctheta = ((v2(1)-x0)*(v3(1)-x0)+(v2(2)-y0)*(v3(2)-y0))
      ctheta = ctheta/dpv2
      ctheta = ctheta/dpv3
      theta = acos(ctheta)
      D = theta
      u = D/2.0d0
      v = D/2.0d0
      sum = 0.0d0
      do i = 1,nquad
	 phi = u*xnodes(i)+v
         a12 = -cos(phi)*e11-sin(phi)*e21
         a22 = -cos(phi)*e12-sin(phi)*e22
	 det = a11*a22- a12*a21 
	 rho = (-a21*rhs1+a11*rhs2)/det
	 zarg = eye*zk*rho
	 if ( abs(zarg) .lt. 1.0d-2) then
	    zz = rho*(1.0d0+zarg/2.0d0+zarg**2/6.0d0+
     1                 zarg**3/24.0d0+zarg**4/120.0d0)
         else
	    zz = (exp(zarg)-1.0d0)/(eye*zk)
         endif
         sum = sum + D*wts(i)*zz/2.0d0
      enddo
c
c     process subtriangle II: (P-V1-V2)
c
      a11 = v2(1)-v1(1)
      a21 = v2(2)-v1(2)
      e11 = v1(1)-x0
      e12 = v1(2)-y0
      e11 = e11/dpv1
      e12 = e12/dpv1
      e21 = -e12
      e22 =  e11
      rhs1 = x0-v1(1)
      rhs2 = y0-v1(2)
      ctheta = ((v1(1)-x0)*(v2(1)-x0)+(v1(2)-y0)*(v2(2)-y0))
      ctheta = ctheta/dpv1
      ctheta = ctheta/dpv2
      theta = acos(ctheta)
      D = theta
      u = D/2.0d0
      v = D/2.0d0
      do i = 1,nquad
	 phi = u*xnodes(i)+v
         a12 = -cos(phi)*e11-sin(phi)*e21
         a22 = -cos(phi)*e12-sin(phi)*e22
	 det = a11*a22- a12*a21 
	 rho = (-a21*rhs1+a11*rhs2)/det
	 zarg = eye*zk*rho
	 if ( abs(zarg) .lt. 1.0d-2) then
	    zz = rho*(1.0d0+zarg/2.0d0+zarg**2/6.0d0+
     1                 zarg**3/24.0d0+zarg**4/120.0d0)
         else
	    zz = (exp(zarg)-1.0d0)/(eye*zk)
         endif
         sum = sum + D*wts(i)*zz/2.0d0
      enddo
c
c     process subtriangle III: (P-V3-V1)
c
      a11 = v1(1)-v3(1)
      a21 = v1(2)-v3(2)
      e11 = v3(1)-x0
      e12 = v3(2)-y0
      e11 = e11/dpv3
      e12 = e12/dpv3
      e21 = -e12
      e22 =  e11
      rhs1 = x0-v3(1)
      rhs2 = y0-v3(2)
      ctheta = ((v3(1)-x0)*(v1(1)-x0)+(v3(2)-y0)*(v1(2)-y0))
      ctheta = ctheta/dpv3
      ctheta = ctheta/dpv1
      theta = acos(ctheta)
      D = theta
      u = D/2.0d0
      v = D/2.0d0
      do i = 1,nquad
	 phi = u*xnodes(i)+v
         a12 = -cos(phi)*e11-sin(phi)*e21
         a22 = -cos(phi)*e12-sin(phi)*e22
	 det = a11*a22- a12*a21 
	 rho = (-a21*rhs1+a11*rhs2)/det
	 zarg = eye*zk*rho
	 if ( abs(zarg) .lt. 1.0d-2) then
	    zz = rho*(1.0d0+zarg/2.0d0+zarg**2/6.0d0+
     1                 zarg**3/24.0d0+zarg**4/120.0d0)
         else
	    zz = (exp(zarg)-1.0d0)/(eye*zk)
         endif
         sum = sum + D*wts(i)*zz/2.0d0
      enddo
      zval = sum
      zvalx = 0.0d0
      zvaly = 0.0d0
      zvalz = 0.0d0
      return
      end
c
c
c
c
c
c
c
