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
c     This is the first release of the FMM3D library, together with
c     associated subroutines, which computes N-body interactions
c     and layer potentials governed by the Laplace or Helmholtz equations.
c  
c
c
        program test_htria
        implicit real *8 (a-h,o-z)
c
c     Testing code for FMM: tests evaluation of single and double
c     layer potentials against O(N^2) direct method using quadratures 
c     on all triangles
c
c     We also test Green's identity on a surface:
c
c     u = S (du/dn) - D (u)  or, taking  the limit to the boundary,
c
c     1/2 u = S (du/dn) - D^* (u)  on surface, where D^* denotes the 
c     principal value integral.
c
c     Convention: the normal vectors are outward from the sphere - hence INTO
c     the domain of interest (the exterior).
c
c     u, du/dn are computed from 
c
c     call lpotfld3d(iff,source,ch,centroids(1,i),potd,fldd)
c
c     using a single interior source. 
c     Since fldd = - grad(phi) and trinorm is -(exterior normal) with
c     respect to the domain, we have (fldd \cdot trinorm = du/dn)
c     as conventionally understood.
c
c     Since triangles are oriented so that normals point into the exterior
c     domain, signs must be flipped in Greens identity: i.e.
c
c     1/2 u =  -S (du/dn) + D^* (u)  on surface.
c
c     taking the limit from the exterior.
c
      parameter (maxvert = 150000) 
      parameter (ntrimax =  50000) 
      integer itrivert(3,maxvert)
      real *8 verts(3,maxvert)
      real *8 centroids(3,ntrimax)
      real *8 triangles(3,3,ntrimax)
      real *8 trinorm(3,ntrimax)
      real *8 triarea(ntrimax)
      real *8 source(3,maxvert)
      real *8 dipvec(3,maxvert)
c
      complex *16 ch, charge(ntrimax)
      complex *16 dipstr(ntrimax)
      complex *16 pot(ntrimax)
      complex *16 potn(ntrimax)
      complex *16 fld(3,ntrimax)
      complex *16 pothalf(ntrimax)
      complex *16 pot2(ntrimax)
      complex *16 pot3(ntrimax)
      complex *16 potn2(ntrimax)
      complex *16 fld2(3,ntrimax)
      complex *16 fld3(3,ntrimax)
      complex *16 pottestfmm(1000)
      complex *16 fldtestfmm(3,1000)
      complex *16 pottestdir(1000)
      complex *16 fldtestdir(3,1000)
      complex *16 ptemp,ftemp(3)
      complex *16 zk,ima,potd,fldd(3)
c
      real *8 target(3,2 000 000)
      complex *16 pottarg(2 000 000)
      complex *16 fldtarg(3,2 000 000)
c
      real *8 center(3),corners(3,8)
      real *8 iseptarg(ntrimax)
c
      real *8 target1(3,2 000 000)
      complex *16 pottarg1(2 000 000)
      complex *16 fldtarg1(3,2 000 000)
c
      data ima/(0.0d0,1.0d0)/
c
c
      done=1
      pi=4*atan(done)
      zk= 1.0d0 + ima*0.1d0
c       
c     Initialize simple printing routines. The parameters to prini
c     define output file numbers using standard Fortran conventions.
c
c     Calling prini(6,13) causes printing to the screen and to 
c     file fort.13.     
c
      call prini(6,13)
c
c
c     Get geometry in .tri format.
c        
      ir = 17
      open (unit=ir,file='sphere2880.a.tri')
c
      call atriread(ir,verts,nverts,itrivert,ntri)
c
c     create triangles, centroids, and normals from tri format data
c
      do itri = 1,ntri
         v11 = verts(1,itrivert(1,itri))
         v21 = verts(2,itrivert(1,itri))
         v31 = verts(3,itrivert(1,itri))
         v12 = verts(1,itrivert(2,itri))
         v22 = verts(2,itrivert(2,itri))
         v32 = verts(3,itrivert(2,itri))
         v13 = verts(1,itrivert(3,itri))
         v23 = verts(2,itrivert(3,itri))
         v33 = verts(3,itrivert(3,itri))
         
         triangles(1,1,itri) = v11
         triangles(2,1,itri) = v21
         triangles(3,1,itri) = v31
         triangles(1,2,itri) = v12
         triangles(2,2,itri) = v22
         triangles(3,2,itri) = v32
         triangles(1,3,itri) = v13
         triangles(2,3,itri) = v23
         triangles(3,3,itri) = v33
         centroids(1,itri) = (v11+v12+v13)/3
         centroids(2,itri) = (v21+v22+v23)/3
         centroids(3,itri) = (v31+v32+v33)/3
         call triangle_norm(triangles(1,1,itri),trinorm(1,itri))
         call triangle_area(triangles(1,1,itri),triarea(itri))
      enddo
c
      call prinf('ntri=*',ntri,1)
c
      iprec=1
      call prinf('iprec=*',iprec,1)
c
      ifpot=1
      iffld=1
      ifslp=1
      ifdlp=1
      ifpottarg=1
      iffldtarg=1

      call prinf('ifpot=*',ifpot,1)
      call prinf('iffld=*',iffld,1)
      call prinf('ifslp=*',ifslp,1)
      call prinf('ifdlp=*',ifdlp,1)
c
      if( ifslp .eq. 1 ) then
         do i=1,ntri
            charge(i)=1
         enddo
      endif
c
      if( ifdlp .eq. 1 ) then
         do i=1,ntri
            dipstr(i)=1
         enddo
      endif
c
      do i=1,ntri
         target(1,i)=centroids(1,i) - 2
         target(2,i)=centroids(2,i)
         target(3,i)=centroids(3,i)
      enddo
      ntarget = ntri
c
c     set source in interior of sphere (external to domain of interest)
c
      source(1,1) = 0.0d0
      source(2,1) = 0.2d0
      source(3,1) = 0.0d0
c       
      do i=1,ntri
         iff = 1
         ch = 1.0d0
         call hpotfld3d(iff,source(1,1),ch,centroids(1,i),zk,potd,fldd)
c       
c     define u as dipole strength and fld \cdot normal = -du/dn as charge strength.
c
         dipstr(i) = potd/(4*pi)
         pothalf(i) = potd/2
         charge(i) = fldd(1)*trinorm(1,i) +
     1              fldd(2)*trinorm(2,i) +
     1              fldd(3)*trinorm(3,i)
         charge(i) = charge(i)/(4*pi)
      enddo
c
ccc      call prin2('half of pot is =*',pothalf,2*ntri)
ccc      call prin2('charge is =*',charge,2*ntri)
ccc      call prin2('dipstr is =*',dipstr,2*ntri)
        call prinf('ntarget=*',ntarget,1)
c
      t1=second()
C$        t1=omp_get_wtime()
c
      call hfmm3dtriatarg(ier,iprec,zk,ntri,
     $   triangles,trinorm,centroids,
     $   ifslp,charge,ifdlp,dipstr,trinorm,
     $   ifpot,pot,iffld,fld,
     $   ntarget,target,ifpottarg,pottarg,iffldtarg,fldtarg)
c
      t2=second()
C$        t2=omp_get_wtime()        
c
c
      if( iffld .eq. 1 ) then
        do i=1,ntri
           potn(i) = fld(1,i)*trinorm(1,i) +
     1        fld(2,i)*trinorm(2,i) +
     1        fld(3,i)*trinorm(3,i)
        enddo
      endif
c
      call prin2('after fmm, time (sec)=*',t2-t1,1)
      call prin2('after fmm, speed (triangles+targets/sec)=*',
     $     (ntri+ntarget)/(t2-t1),1)
c
c
      err = 0.0d0
      err7 = 0.0d0
      denom = 0.0d0
      denom7 = 0.0d0
      do i=1,ntri
         err = err + abs(pothalf(i)-pot(i))**2
         denom = denom + abs(pothalf(i))**2
         pot2(i) = fld(1,i)*trinorm(1,i) +
     1      fld(2,i)*trinorm(2,i) +
     1      fld(3,i)*trinorm(3,i)
         pot2(i) = potn(i)/ (2*pi)
         err7 = err7 + abs(pot2(i)-charge(i))**2
         denom7 = denom7 + abs(charge(i))**2
      enddo
      err = sqrt(err/denom)
      err7 = sqrt(err7/denom7)
      call prin2(' Greens identity error in ext is *',err,1)
      call prin2(' Greens identity error in dudn is *',err7,1)
c        
c
      do i=1,ntri
         if (ifpot.eq.1) pot2(i)=0
         if (iffld .eq. 1) then
            fld2(1,i)=0
            fld2(2,i)=0
            fld2(3,i)=0
         endif
      enddo
c        
      m=min(ntri,10)
c
      t1=second()
C$        t1=omp_get_wtime()
c
      nqtri = 6
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,ftemp,ifl) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
      do j=1,m
         if (ifslp .eq. 1) ifl = 0
         if (ifslp .eq. 2) ifl = 1
         if (ifslp .ne. 0) then
            call direct3dtrihelms2(ifl,j,ntri,zk,nqtri,centroids,charge,
     $          triangles,ptemp,ftemp)
            if( ifpot .eq. 1 ) pot2(j) = pot2(j) + ptemp
            if( iffld .eq. 1 ) then
               fld2(1,j) = fld2(1,j) + ftemp(1)
               fld2(2,j) = fld2(2,j) + ftemp(2)
               fld2(3,j) = fld2(3,j) + ftemp(3)
            endif
         endif
         if (ifdlp .ne. 0) then
            call direct3dtrihelmd2(ifl,j,ntri,zk,nqti,centroids,dipstr,
     $        triangles,trinorm,ptemp,ftemp)
            if( ifpot .eq. 1) pot2(j) = pot2(j) + ptemp
            if( iffld .eq. 1) then
               fld2(1,j) = fld2(1,j) + ftemp(1)
               fld2(2,j) = fld2(2,j) + ftemp(2)
               fld2(3,j) = fld2(3,j) + ftemp(3)
            endif
         endif
      enddo
C$OMP END PARALLEL DO
c
      if (iffld.eq.1) then
         do i=1,m
            potn2(i) = fld2(1,i)*trinorm(1,i) +
     1        fld2(2,i)*trinorm(2,i) +
     1        fld2(3,i)*trinorm(3,i)
         enddo
      endif
c       
      t2=second()
C$        t2=omp_get_wtime()
c
      ifprint = 1
      call prinf(' ifpot=*',ifpot,1)
c
      if (ifprint .eq. 1) then
         if( ifpot .eq. 1 ) call prin2('from fmm, pot=*',pot,2*m)
         if( ifpot .eq. 1 ) call prin2('directly, pot2=*',pot2,2*m)
         if( iffld .eq. 1 ) call prin2('from fmm, fld=*',fld,6*m)
         if( iffld .eq. 1 ) call prin2('directly, fld2=*',fld2,6*m)
      endif
c
      call prin2('directly, time (sec)=*',
     $     (t2-t1)*dble(ntri)/dble(m),1)
      call prin2('directly, speed (triangles/sec)=*',
     $     m/(t2-t1),1)
c       
      if (ifpot .eq. 1) then
         call h3derror(pot,pot2,m,aerr,rerr)
         call prin2('abs error in potential=*',aerr,1)
         call prin2('relative error in potential=*',rerr,1)
      endif
c       
      if (iffld .eq. 1) then
         call h3derror(fld,fld2,3*m,aerr,rerr)
         call prin2('abs error in field=*',aerr,1)
         call prin2('relative error in field=*',rerr,1)
      endif
c
      do i=1,ntarget
         if (ifpottarg.ne.0) pot2(i)=0
         if (iffldtarg .eq. 1) then
            fld2(1,i)=0
            fld2(2,i)=0
            fld2(3,i)=0
         endif
      enddo
c        
      m=min(ntarget,10)
      t1=second()
C$        t1=omp_get_wtime()
c
      nqtri = 6
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,ftemp,ifl) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
      do j=1,m
         if (ifslp .eq. 1) ifl = 0
         if (ifslp .eq. 2) ifl = 1
         if (ifslp .ne. 0) then
            call direct3dtritarghelms2(ifl,ntri,target(1,j),zk,nqtri,
     $           charge,triangles,ptemp,ftemp)
            if( ifpottarg .eq. 1 ) pot2(j) = pot2(j) + ptemp
            if( iffldtarg .eq. 1 ) then
               fld2(1,j) = fld2(1,j) + ftemp(1)
               fld2(2,j) = fld2(2,j) + ftemp(2)
               fld2(3,j) = fld2(3,j) + ftemp(3)
            endif
         endif
         if (ifdlp .ne. 0) then
            call direct3dtritarghelmd2(ifl,ntri,target(1,j),zk,nqtri,
     $         dipstr,triangles,trinorm,ptemp,ftemp)
            if( ifpottarg .eq. 1) pot2(j) = pot2(j) + ptemp
            if( iffldtarg .eq. 1) then
              fld2(1,j) = fld2(1,j) + ftemp(1)
              fld2(2,j) = fld2(2,j) + ftemp(2)
              fld2(3,j) = fld2(3,j) + ftemp(3)
            endif
        endif
      enddo
C$OMP END PARALLEL DO
c
      t2=second()
C$        t2=omp_get_wtime()
c
      ifprint = 1
c
      if (ifprint .eq. 1) then
         if( ifpottarg .eq. 1 ) call prin2('fmm, pottarg=*',pottarg,2*m)
         if( ifpottarg .eq. 1 ) call prin2('directly, pot2=*',pot2,2*m)
         if( iffldtarg .eq. 1 ) call prin2('fmm, fldtarg=*',fldtarg,6*m)
         if( iffldtarg .eq. 1 ) call prin2('directly, fld2=*',fld2,6*m)
      endif
c
      call prin2('directly, time (sec)=*',
     $     (t2-t1)*dble(ntri)/dble(m),1)
      call prin2('directly, speed (targets/sec)=*',
     $     m/(t2-t1),1)
c       
      if (ifpottarg .eq. 1) then
         call h3derror(pottarg,pot2,m,aerr,rerr)
         call prin2('abs error in target potential=*',aerr,1)
         call prin2('relative error in target potential=*',rerr,1)
      endif
c       
      if (iffldtarg .eq. 1) then
         call h3derror(fldtarg,fld2,3*m,aerr,rerr)
         call prin2('abs error in target field=*',aerr,1)
         call prin2('relative error in target field=*',rerr,1)
      endif
c       
      call prinf('=== test multipole expansion in far field*',i,0)
c
c     finally, test the far field evaluation via multipole expansion 
c     by congtructing a distant sphere of targets.
c
      t1=second()
C$        t1=omp_get_wtime()
c       
c
      do i=1,ntarget
         target(1,i)=centroids(1,i) - 30
         target(2,i)=centroids(2,i)
         target(3,i)=centroids(3,i)
      enddo
c      
      call hfmm3dtriampftarg(ier,iprec,zk,ntri,
     $     triangles,trinorm,centroids,
     $     ifslp,charge,ifdlp,dipstr,trinorm,
     $     ntarget,target,ifpottarg,pot3,iffldtarg,fld3)
c       
c
      t2=second()
C$        t2=omp_get_wtime()
c
      call prin2('via mp, time (sec)=*',
     $     (t2-t1),1)
      call prin2('via mp, speed (targets/sec)=*',
     $     ntarget/(t2-t1),1)
c       
      if (ifprint .eq. 1) then
c
         if( ifpottarg .eq. 1 ) 
     $     call prin2('via mp, pottarg=*',pot3,2*m)
         if( iffldtarg .eq. 1 ) 
     $     call prin2('via mp, fldtarg=*',fld3,6*m)
      endif
c
      call prinf(' ifpot=*',ifpot,1)
      do i=1,ntarget
         if (ifpottarg.eq.1) pot2(i)=0
         if (iffldtarg .eq. 1) then
            fld2(1,i)=0
            fld2(2,i)=0
            fld2(3,i)=0
         endif
      enddo
      t1=second()
C$        t1=omp_get_wtime()
c
      nqtri = 12
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,ftemp,ifl) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
      do j=1,m
         if (ifslp .eq. 1) ifl = 0
         if (ifslp .eq. 2) ifl = 1
         if (ifslp .ne. 0) then
            call direct3dtritarghelms2(ifl,ntri,target(1,j),zk,nqtri,
     $           charge,triangles,ptemp,ftemp)
            if( ifpottarg .eq. 1 ) pot2(j) = pot2(j) + ptemp
            if( iffldtarg .eq. 1 ) then
               fld2(1,j) = fld2(1,j) + ftemp(1)
               fld2(2,j) = fld2(2,j) + ftemp(2)
               fld2(3,j) = fld2(3,j) + ftemp(3)
            endif
         endif
         if (ifdlp .ne. 0) then
            call direct3dtritarghelmd2(ifl,ntri,target(1,j),zk,nqtri,
     $         dipstr,triangles,trinorm,ptemp,ftemp)
            if( ifpottarg .eq. 1) pot2(j) = pot2(j) + ptemp
            if( iffldtarg .eq. 1) then
              fld2(1,j) = fld2(1,j) + ftemp(1)
              fld2(2,j) = fld2(2,j) + ftemp(2)
              fld2(3,j) = fld2(3,j) + ftemp(3)
            endif
        endif
      enddo
C$OMP END PARALLEL DO
c
      t2=second()
C$        t2=omp_get_wtime()
c

      if (ifpottarg .eq. 1) then
        call h3derror(pot3,pot2,m,aerr,rerr)
        call prin2('abs error in target potential=*',aerr,1)
        call prin2('relative error in target potential=*',rerr,1)
      endif
c       
      if (iffldtarg .eq. 1) then
        call h3derror(fld3,fld2,3*m,aerr,rerr)
        call prin2('abs error in target field=*',aerr,1)
        call prin2('relative error in target field=*',rerr,1)
      endif
c
      stop
      end
c
c
c
c
c
c
      subroutine h3derror(pot1,pot2,n,ae,re)
      implicit real *8 (a-h,o-z)
c
c     evaluate absolute and relative errors
c
      complex *16 pot1(n),pot2(n)
c
      d=0
      a=0
c       
      do i=1,n
         d=d+abs(pot1(i)-pot2(i))**2
         a=a+abs(pot1(i))**2
      enddo
c       
      d=d/n
      d=sqrt(d)
      a=a/n
      a=sqrt(a)
c       
      ae=d
      re=d/a
c       
      return
      end
c
c
