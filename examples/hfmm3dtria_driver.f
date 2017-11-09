cc Copyright (C) 2009-2012: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This software is being released under a modified FreeBSD license
cc (see COPYING in home directory). 
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
      ms = m
      mt = 0
      nqtri = 6

      call h3dtriadirect_test(ms,mt,nqtri,zk,ntri,
     $   triangles,trinorm,centroids,
     $   ifslp,charge,ifdlp,dipstr,trinorm,
     $   ifpot,pot2,iffld,fld2,
     $   ntarget,target,ifpottarg,pot3,iffldtarg,fld3)

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



        subroutine h3dtriadirect_test(ms,mt,nqtri,zk,nsource,
     $     triaflat,trianorm,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ntarget,
     $     target,ifpottarg,pottarg,iffldtarg,fldtarg)
        implicit real *8 (a-h,o-z)
c
c       Helmholtz interactions in R^3: evaluate all pairwise triangle
c       interactions (including self-interaction) + interactions with targets
c       via direct O(N^2) algorithm.
c
c       This is the principal subroutine for evaluating 
c       Helmholtz layer potentials on (flat) triangulated surfaces.
c       It permits the evaluation of a single layer potential
c       with piecewise constant density defined by <<charge>>
c       and a dipole layer potential with piecewise constant density 
c       and dipole orientation defined by <<dipstr,dipvec>>.
c
c       We use (exp(ikr)/r) for the Green's function,
c       without the (1/4 pi) scaling.  Self-interactions are included.
c   
c       It is capable of evaluating the layer potentials either on 
c       or off the surface (or both).            
c
c       The quadratures are approximate: we subtract the 1/r singular
c       part and proceed with a quadrature for smooth functions of each
c       triangle. The accuracy of this scheme is controlled by the
c       user-defined parameter nqtri.
c
c       INPUT PARAMETERS:
c       
c       ms: integer: numer of sources to be tested
c       mt: integer: numer of targets to be tested
c
c       nqtri: integer: number of quadrature nodes.
c          Suggested values for nqtri, assuming resonably regular triangulation.
c          iprec:  FMM precision flag
c                 -2 => tolerance =.5d0
c                 -1 => tolerance =.5d-1
c                  0 => tolerance =.5d-2
c                  1 => tolerance =.5d-3
c                  2 => tolerance =.5d-6
c                  3 => tolerance =.5d-9
c                  4 => tolerance =.5d-12
c                  5 => tolerance =.5d-15
c          if( iprec .eq. -2 ) nqtri=1
c          if( iprec .eq. -1 ) nqtri=2
c          if( iprec .eq.  0 ) nqtri=3
c          if( iprec .ge.  1 ) nqtri=6
c       zk: complex *16: Helmholtz parameter
c       nsource: integer:  number of triangles
c       triaflat: real *8 (3,3,nsource): triangle coordinate array
c       trianorm: real *8 (3,nsource): triangle normals
c       source: real *8 (3,nsource):  triangle centroids
c       ifcharge:  single layer potential (SLP) flag
c                  ifcharge = 1   =>  include SLP contribution
c                                     otherwise do not
c                  ifcharge = 2   =>  include SLP contribution and
c                                     subtract Laplace SLP
c       charge: complex *16 (nsource): piecewise constant SLP strength
c       ifdipole:  dipole layer potential (DLP) flag
c                  ifdipole = 1   =>  include DLP contribution
c                                     otherwise do not
c                  ifdipole = 2   =>  include DLP contribution and
c                                     subtract Laplace SLP
c       dipstr: complex *16 (nsource): piecewise constant DLP strengths
c       dipvec: real *8 (3,nsource): piecewise constant dipole orientation 
c                                    vectors. 
c
c           NOTE: In the present version, dipvec MUST BE SET EQUAL
c                 to the triangle normal. It is here as an additional
c                 parameter for future use, where an arbitrarily 
c                 oriented dipole vector is permitted.
c
c       ifpot:  potential flag (1=compute potential, otherwise no)
c       iffld:  field flag (1=compute field, otherwise no)
c       ntarget: integer:  number of targets
c       target: real *8 (3,ntarget):  target locations
c       ifpottarg:  target potential flag 
c                   (1=compute potential, otherwise no)
c       iffldtarg:  target field flag 
c                   (1=compute field, otherwise no)
c
c       OUTPUT PARAMETERS:
c
c       pot: complex *16 (nsource): potential at triangle centroids
c       fld: complex *16 (3,nsource): field (-gradient) at triangle centroids 
c       pottarg: complex *16 (ntarget): potential at target locations 
c       fldtarg: complex *16 (3,ntarget): field (-gradient) at target locations 
c
c
        real *8 triaflat(3,3,1),trianorm(3,1)
c
        real *8 source(3,1),dipvec(3,1)
        complex *16 charge(1),dipstr(1),zk
        real *8 target(3,1)
c
        complex *16 pot(1),fld(3,1),pottarg(1),fldtarg(3,1)
        complex *16 ptemp,ftemp(3)
c
c
        do i=1,nsource
        if( ifpot .eq. 1) pot(i)=0
        if( iffld .eq. 1) then
           fld(1,i)=0
           fld(2,i)=0
           fld(3,i)=0
        endif
        enddo
c       
        do i=1,ntarget
        if( ifpottarg .eq. 1) pottarg(i)=0
        if( iffldtarg .eq. 1) then
           fldtarg(1,i)=0
           fldtarg(2,i)=0
           fldtarg(3,i)=0
        endif
        enddo
c
        ione = 1
        if( ifpot .eq. 1 .or. iffld .eq. 1 ) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,ftemp,ifl) 
        do j=1,ms
        do i=1,nsource
c
        if (ifcharge.eq.1) ifl=0
        if (ifcharge.eq.2) ifl=1
        if (ifcharge.ne.0) then
        if( i .eq. j ) then
        call direct3dtrihelms2(ifl,ione,ione,zk,nqtri,
     1     source(1,j),charge(j),triaflat(1,1,j),ptemp,ftemp)
        else
        call direct3dtritarghelms3(ifl,ione,source(1,j),zk,
     $     nqtri,charge(i),triaflat(1,1,i),ifpot,ptemp,iffld,ftemp)
        endif
        if (ifpot.eq.1) pot(j)=pot(j)+ptemp
        if (iffld.eq.1) then
        fld(1,j)=fld(1,j)+ftemp(1)
        fld(2,j)=fld(2,j)+ftemp(2)
        fld(3,j)=fld(3,j)+ftemp(3)
        endif
        endif
c
        if (ifdipole.eq.1) ifl=0
        if (ifdipole.eq.2) ifl=1
        if (ifdipole.ne.0) then
        if( i .eq. j ) then
        call direct3dtrihelmd2(ifl,ione,ione,zk,nqtri,
     $     source(1,j),dipstr(j),triaflat(1,1,j),
     $     trianorm(1,j),ptemp,ftemp)
        else
        call direct3dtritarghelmd3(ifl,ione,source(1,j),zk,
     $     nqtri,dipstr(i),triaflat(1,1,i),
     $     trianorm(1,i),ifpot,ptemp,iffld,ftemp)
        endif
        if (ifpot.eq.1) pot(j)=pot(j)+ptemp
        if (iffld.eq.1) then
        fld(1,j)=fld(1,j)+ftemp(1)
        fld(2,j)=fld(2,j)+ftemp(2)
        fld(3,j)=fld(3,j)+ftemp(3)
        endif
        endif
c
        enddo
        enddo
C$OMP END PARALLEL DO
        endif

        if( ifpottarg .eq. 1 .or. iffldtarg .eq. 1 ) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,ftemp,ifl) 
        do j=1,mt
        do i=1,nsource
c
        if (ifcharge.eq.1) ifl=0
        if (ifcharge.eq.2) ifl=1
        if (ifcharge.ne.0) then
        call direct3dtritarghelms3(ifl,ione,target(1,j),zk,
     $     nqtri,charge(i),triaflat(1,1,i),
     $     ifpottarg,ptemp,iffldtarg,ftemp)
        if (ifpottarg.eq.1) pottarg(j)=pottarg(j)+ptemp
        if (iffldtarg.eq.1) then
        fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
        fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
        fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
        endif
        endif
c
        if (ifdipole.eq.1) ifl=0
        if (ifdipole.eq.2) ifl=1
        if (ifdipole.ne.0) then
        call direct3dtritarghelmd3(ifl,ione,target(1,j),zk,
     $     nqtri,dipstr(i),triaflat(1,1,i),
     $     trianorm(1,i),ifpottarg,ptemp,iffldtarg,ftemp)
        if (ifpottarg.eq.1) pottarg(j)=pottarg(j)+ptemp
        if (iffldtarg.eq.1) then
        fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
        fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
        fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
        endif
        endif
c
        enddo
        enddo
C$OMP END PARALLEL DO
        endif
c
        return
        end
c
c
c
c
c
