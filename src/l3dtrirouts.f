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
c    Low level FMM routines for triangles in Laplace regime.
c    (forming multipole expansions and forming taylor expansions).
c
c       SLP - single layer potential
c       DLP - double layer potential
c
c      l3dformmptrisone: form a multipole expansion due to SLP 
c                      on a single triangle.
c
c      l3dformmptris_add: *increment* a multipole expansion due to SLP 
c                      on a collection of triangles.
c
c      l3dformmptris: form a multipole expansion due to SLP 
c                      on a collection of triangles.
c
c      l3dformmptridone: form a multipole expansion due to DLP 
c                      on a single triangle.
c
c      l3dformmptrid_add: *increment* a multipole expansion due to DLP 
c                      on a collection of triangles.
c
c      l3dformmptrid: form a multipole expansion due to DLP 
c                      on a collection of triangles.
c
c      l3dformtatrisone: form a local expansion due to SLP 
c                      on a single triangle.
c
c      l3dformtatris_add: *increment* a local expansion due to SLP 
c                      on a collection of triangles.
c
c      l3dformtatris: form a local expansion due to SLP 
c                      on a collection of triangles.
c
c      l3dformtatridone: form a local expansion due to DLP 
c                      on a single triangle.
c
c      l3dformtatrid_add: *increment* a local expansion due to DLP 
c                      on a collection of triangles.
c
c      l3dformtatrid: form a local expansion due to DLP 
c                      on a collection of triangles.
c
c
c      l3dformmptrisone2: form a multipole expansion due to SLP 
c                      on a single triangle (OPTIMIZED VERSION).
c
c      l3dformmptris2_add: *increment* a multipole expansion due to SLP 
c                      on a collection of triangles (OPTIMIZED VERSION).
c
c      l3dformmptridone2: form a multipole expansion due to DLP 
c                      on a single triangle (OPTIMIZED VERSION).
c
c      l3dformmptrid2_add: *increment* a multipole expansion due to DLP 
c                      on a collection of triangles (OPTIMIZED VERSION).
c
c      l3dformtatrisone2: form a local expansion due to SLP 
c                      on a single triangle (OPTIMIZED VERSION).
c
c      l3dformtatris2_add: *increment* a local expansion due to SLP 
c                      on a collection of triangles (OPTIMIZED VERSION).
c
c      l3dformtatridone2: form a local expansion due to DLP 
c                      on a single triangle (OPTIMIZED VERSION).
c
c      l3dformtatrid2_add: *increment* a local expansion due to DLP 
c                      on a collection of triangles (OPTIMIZED VERSION).
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine l3dformmptrisone(ier, scale, triang,
     1          charge, x0y0z0, norder, nterms, mpole)
      implicit none
c
c     Form  multipole expansion due to SLP on  single triangle.
c
c     INPUT:
c
c     scale        = scaling parameter
c     triang(i,j)  = ith coord of jth vertex
c     charge       = SLP strength
c     x0y0z0       = center of the expansion
c     norder       = order of Gaussian rule on triangle
c     nterms       = order of multipole expansion
c
c     OUTPUT:
c
c     mpole        =  induced multipole expansion 
c     lused        =  amount of workspace w used
c     ier          =  error return code
c
c---------------------------------------------------------------------------
c
      integer  ier, ifinit, norder, nnodes, nterms, lused
      integer  i, j, l, m
c
      real *8  triang(3,3), x0y0z0(3)
      real *8  scale
      real *8  rnodes(2,500), weights(500), zparts(3,500)
      real *8  vert1(2),vert2(2),vert3(2),vertout(3)
      real *8  w(20)
c
      complex *16 eye, charge, ctemp(500)
      complex *16  mpole(0:nterms,-nterms:nterms)
c
      data eye/(0.0d0,1.0d0)/
c
      call tri_ini(triang(1,1),triang(1,2),triang(1,3),w,
     1             vert1,vert2,vert3)
c
      call triasymq(norder,vert1,vert2,vert3,rnodes,weights,nnodes)
ccc      call prin2(' weights are *',weights,nnodes)
c
      do i = 1,nnodes
         vertout(1) = rnodes(1,i)
         vertout(2) = rnodes(2,i)
         vertout(3) = 0
         call tri_bak(w,vertout,zparts(1,i))
         ctemp(i) = charge*weights(i)
      enddo
ccc      call prin2(' ctemp are *',ctemp,2*nnodes)
ccc      call prin2(' zparts are *',zparts,3*nnodes)
ccc      call prin2(' zk is *',zk,2)
ccc      call prin2(' scale is *',scale,1)
ccc      call prin2(' x0y0z0 is *',x0y0z0,3)
ccc      call prinf(' nterms is *',nterms,1)
      call l3dformmp(ier,scale,zparts,ctemp,nnodes,x0y0z0,
     1               nterms,mpole)
      return
      end
c
c
c
c
c
      subroutine l3dformmptris(ier, scale, triang, charge, 
     1          ntri, x0y0z0, norder, nterms, mpole, mptemp)
      implicit none
c
c
c     This subroutine INCREMENTS the multipole expansion to include
c     contributions from SLP on collection of triangles.
c
c     INPUT:
c
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     charge             = array of SLP densities
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of multipole expansion
c     mptemp             = work array to hold temp multipole expansion
c
c     OUTPUT:
c
c     mpole              =  coefficients of multipole expansion 
c                           are INCREMENTED by contributions from each
c                           triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c---------------------------------------------------------------------------
c
      integer  nterms, ntri, ier, lused
      integer  i, norder
      real *8  triang(3,3,1), x0y0z0(3), scale
      complex *16  zk, charge(1)
      complex *16  mpole(0:nterms,-nterms:nterms), mptemp(1)
c
      do i = 1,ntri
         call l3dformmptrisone(ier,scale,triang(1,1,i),charge(i), 
     1          x0y0z0, norder, nterms, mptemp)
         call l3dadd(mptemp,mpole,nterms)
      enddo
      return
      end
c
c
c
      subroutine l3dformmptris_add(ier, scale, triang, charge, 
     1          ntri, x0y0z0, norder, nterms, mpole)
      implicit none
c
c
c     This subroutine INCREMENTS the multipole expansion to include
c     contributions from SLP on collection of triangles.
c
c     INPUT:
c
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     charge             = array of SLP densities
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of multipole expansion
c     mptemp             = work array to hold temp multipole expansion
c
c     OUTPUT:
c
c     mpole              =  coefficients of multipole expansion 
c                           are INCREMENTED by contributions from each
c                           triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c---------------------------------------------------------------------------
c
      integer  nterms, ntri, ier, lused
      integer  i, norder
      real *8  triang(3,3,1), x0y0z0(3), scale
      complex *16  zk, charge(1)
      complex *16  mpole(0:nterms,-nterms:nterms)
      complex *16, allocatable :: mptemp(:,:)
c
      allocate( mptemp(0:nterms,-nterms:nterms) )
c
      do i = 1,ntri
         call l3dformmptrisone(ier,scale,triang(1,1,i),charge(i), 
     1          x0y0z0, norder, nterms, mptemp)
         call l3dadd(mptemp,mpole,nterms)
      enddo
      return
      end
c
c
c
c
c
c
      subroutine l3dformmptridone(ier, scale, triang, trinorm,
     1          dipstr, x0y0z0, norder, nterms, mpole)
      implicit none
c
c     Form multipole expansion due to DLP on single triangle.
c
c     INPUT:
c
c     scale        = scaling parameter
c     triang(i,j)  = ith coord of jth vertex
c     trinorm(i)   = ith coord of normal
c     dipstr       = DLP strength
c     x0y0z0       = center of the expansion
c     norder       = order of Gaussian rule on triangle
c     nterms       = order of multipole expansion
c
c     OUTPUT:
c
c     mpole        =  induced multipole expansion 
c     lused        =  amount of workspace w used
c     ier          =  error return code
c
c---------------------------------------------------------------------------
c
      integer  ier, ifinit, norder, nnodes, nterms, lused
      integer  i, j, l, m
c
      real *8  triang(3,3), trinorm(3), x0y0z0(3)
      real *8  scale
      real *8  rnodes(2,500), weights(500), zparts(3,500)
      real *8  dipvec(3,500)
      real *8  vert1(2),vert2(2),vert3(2),vertout(3)
      real *8  w(20)
c
      complex *16 eye, dipstr, zk, ctemp(500)
      complex *16  mpole(0:nterms,-nterms:nterms)
c
      data eye/(0.0d0,1.0d0)/
c
      call tri_ini(triang(1,1),triang(1,2),triang(1,3),w,
     1             vert1,vert2,vert3)
c
      call triasymq(norder,vert1,vert2,vert3,rnodes,weights,nnodes)
ccc        call prinf('nnodes=*',nnodes,1)
c
        dipvec(1,1)=trinorm(1)
        dipvec(2,1)=trinorm(2)
        dipvec(3,1)=trinorm(3)
c
      do i = 1,nnodes
         vertout(1) = rnodes(1,i)
         vertout(2) = rnodes(2,i)
         vertout(3) = 0
         ctemp(i) = dipstr*weights(i)
         call tri_bak(w,vertout,zparts(1,i))
         dipvec(1,i) = dipvec(1,1)
         dipvec(2,i) = dipvec(2,1)
         dipvec(3,i) = dipvec(3,1)
      enddo
      call l3dformmp_dp(ier,scale,zparts,ctemp,dipvec,nnodes,x0y0z0,
     1               nterms,mpole)
      return
      end
c
c
c
c
c
      subroutine l3dformmptrid(ier,scale,triang,trinorm,dipstr,
     1          ntri,x0y0z0,norder,nterms,mpole,mptemp)
      implicit none
c
c
c     This subroutine INCREMENTS the multipole expansion about x0y0z0 
c     to include contributions from DLP on collection of triangles.
c
c     INPUT:
c
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     dipstr             = array of DLP densities
c     trinorm            = normal to triangle
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of multipole expansion
c     mptemp             = work array to hold temp multipole expansion
c
c     OUTPUT:
c
c     mpole              =  coefficients of multipole expansion 
c                           are INCREMENTED by contributions from each
c                           DLP triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c---------------------------------------------------------------------------
c
      integer  nterms, ntri, ier, lused
      integer  i, norder
      real *8  triang(3,3,1), trinorm(3,1), x0y0z0(3), scale
      complex *16  eye, zk, dipstr(1)
      complex *16  mpole(0:nterms,-nterms:nterms), mptemp(1)
c
      do i = 1,ntri
         call l3dformmptridone(ier,scale,triang(1,1,i),trinorm(1,i),
     1          dipstr(i),x0y0z0, norder,nterms,mptemp)
         call l3dadd(mptemp,mpole,nterms)
      enddo
      return
      end
c
c
      subroutine l3dformmptrid_add(ier,scale,triang,trinorm,dipstr,
     1          ntri,x0y0z0,norder,nterms,mpole)
      implicit none
c
c
c     This subroutine INCREMENTS the multipole expansion about x0y0z0 
c     to include contributions from DLP on collection of triangles.
c
c     INPUT:
c
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     dipstr             = array of DLP densities
c     trinorm            = normal to triangle
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of multipole expansion
c     mptemp             = work array to hold temp multipole expansion
c
c     OUTPUT:
c
c     mpole              =  coefficients of multipole expansion 
c                           are INCREMENTED by contributions from each
c                           DLP triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c---------------------------------------------------------------------------
c
      integer  nterms, ntri, ier, lused
      integer  i, norder
      real *8  triang(3,3,1), trinorm(3,1), x0y0z0(3), scale
      complex *16  eye, zk, dipstr(1)
      complex *16  mpole(0:nterms,-nterms:nterms)
      complex *16, allocatable :: mptemp(:,:)
c
      allocate( mptemp(0:nterms,-nterms:nterms) )
c
      do i = 1,ntri
         call l3dformmptridone(ier,scale,triang(1,1,i),trinorm(1,i),
     1          dipstr(i),x0y0z0, norder,nterms,mptemp)
         call l3dadd(mptemp,mpole,nterms)
      enddo
c
      return
      end
c
c
c
c     FORMTA routines: Individual triangle to local expansion
c     for list 3/4 processing.
c
c
c
      subroutine l3dformtatrisone(ier, scale, triang,
     1          charge, x0y0z0, norder, nterms, local)
      implicit none
c
c     Form local expansion about x0y0z0 due to SLP on single triangle.
c
c     INPUT:
c
c     scale        = scaling parameter
c     triang(i,j)  = ith coord of jth vertex
c     charge       = SLP strength
c     x0y0z0       = center of the expansion
c     norder       = order of Gaussian rule on triangle
c     nterms       = order of multipole expansion
c
c     OUTPUT:
c
c     local        =  induced local expansion 
c     lused        =  amount of workspace w used
c     ier          =  error return code
c
c---------------------------------------------------------------------------
c
      integer  ier, ifinit, norder, nnodes, nterms, lused
      integer  i, j, l, m
c
      real *8  triang(3,3), x0y0z0(3)
      real *8  scale
      real *8  rnodes(2,500), weights(500), zparts(3,500)
      real *8  vert1(2),vert2(2),vert3(2),vertout(3)
      real *8  w(20)
c
      complex *16 eye, charge, zk, ctemp(500)
      complex *16  local(0:nterms,-nterms:nterms)
c
      data eye/(0.0d0,1.0d0)/
c
      call tri_ini(triang(1,1),triang(1,2),triang(1,3),w,
     1             vert1,vert2,vert3)
c
      call triasymq(norder,vert1,vert2,vert3,rnodes,weights,nnodes)
ccc      call prin2(' weights are *',weights,nnodes)
c
      do i = 1,nnodes
         vertout(1) = rnodes(1,i)
         vertout(2) = rnodes(2,i)
         vertout(3) = 0
         call tri_bak(w,vertout,zparts(1,i))
         ctemp(i) = charge*weights(i)
      enddo
ccc      call prin2(' ctemp are *',ctemp,2*nnodes)
ccc      call prin2(' zparts are *',zparts,3*nnodes)
ccc      call prin2(' zk is *',zk,2)
ccc      call prin2(' scale is *',scale,1)
ccc      call prin2(' x0y0z0 is *',x0y0z0,3)
ccc      call prinf(' nterms is *',nterms,1)
      call l3dformta(ier,scale,zparts,ctemp,nnodes,x0y0z0,
     1               nterms,local)
      return
      end
c
c
c
c
c
      subroutine l3dformtatris(ier, scale, triang, charge, 
     1          ntri, x0y0z0, norder, nterms, local, mptemp)
      implicit none
c
c
c     This subroutine INCREMENTS the local expansion to include
c     contributions from SLP on collection of triangles.
c
c     INPUT:
c
c     zk                 = Helmholtz parameter
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     charge             = array of SLP densities
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of multipole expansion
c     mptemp             = work array to hold temp multipole expansion
c
c     OUTPUT:
c
c     local              =  coefficients of local expansion 
c                           are INCREMENTED by contributions from each
c                           SLP triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c---------------------------------------------------------------------------
c
      integer  nterms, ntri, ier, lused
      integer  i, norder
      real *8  triang(3,3,1), x0y0z0(3), scale
      complex *16  zk, charge(1)
      complex *16  local(0:nterms,-nterms:nterms), mptemp(1)
c
      do i = 1,ntri
         call l3dformtatrisone(ier,scale,triang(1,1,i),charge(i), 
     1          x0y0z0, norder, nterms, mptemp)
         call l3dadd(mptemp,local,nterms)
      enddo
      return
      end
c
c
      subroutine l3dformtatris_add(ier, scale, triang, charge, 
     1          ntri, x0y0z0, norder, nterms, local)
      implicit none
c
c
c     This subroutine INCREMENTS the local expansion to include
c     contributions from SLP on collection of triangles.
c
c     INPUT:
c
c     zk                 = Helmholtz parameter
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     charge             = array of SLP densities
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of multipole expansion
c     mptemp             = work array to hold temp multipole expansion
c
c     OUTPUT:
c
c     local              =  coefficients of local expansion 
c                           are INCREMENTED by contributions from each
c                           SLP triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c---------------------------------------------------------------------------
c
      integer  nterms, ntri, ier, lused
      integer  i, norder
      real *8  triang(3,3,1), x0y0z0(3), scale
      complex *16  zk, charge(1)
      complex *16  local(0:nterms,-nterms:nterms)
      complex *16, allocatable :: mptemp(:,:)
c
      allocate( mptemp(0:nterms,-nterms:nterms) )
c
      do i = 1,ntri
         call l3dformtatrisone(ier,scale,triang(1,1,i),charge(i), 
     1          x0y0z0, norder, nterms, mptemp)
         call l3dadd(mptemp,local,nterms)
      enddo
      return
      end
c
c
c
c
c
c
      subroutine l3dformtatridone(ier, scale, triang, trinorm,
     1          dipstr, x0y0z0, norder, nterms, local)
      implicit none
c
c     Form local expansion about x0y0z0 due to DLP on single triangle.
c
c     INPUT:
c
c     scale        = scaling parameter
c     triang(i,j)  = ith coord of jth vertex
c     trinorm(i)   = ith coord of normal
c     dipstr       = DLP strength
c     x0y0z0       = center of the expansion
c     norder       = order of Gaussian rule on triangle
c     nterms       = order of multipole expansion
c
c     OUTPUT:
c
c     local        =  induced local expansion 
c     lused        =  amount of workspace w used
c     ier          =  error return code
c
c---------------------------------------------------------------------------
c
      integer  ier, ifinit, norder, nnodes, nterms, lused
      integer  i, j, l, m
c
      real *8  triang(3,3), trinorm(3), x0y0z0(3)
      real *8  scale
      real *8  rnodes(2,500), weights(500), zparts(3,500)
      real *8  dipvec(3,500)
      real *8  vert1(2),vert2(2),vert3(2),vertout(3)
      real *8  w(20)
c
      complex *16 eye, dipstr, zk, ctemp(500)
      complex *16  local(0:nterms,-nterms:nterms)
c
      data eye/(0.0d0,1.0d0)/
c
      call tri_ini(triang(1,1),triang(1,2),triang(1,3),w,
     1             vert1,vert2,vert3)
c
      call triasymq(norder,vert1,vert2,vert3,rnodes,weights,nnodes)
c
        dipvec(1,1)=trinorm(1)
        dipvec(2,1)=trinorm(2)
        dipvec(3,1)=trinorm(3)
c
      do i = 1,nnodes
         vertout(1) = rnodes(1,i)
         vertout(2) = rnodes(2,i)
         vertout(3) = 0
         ctemp(i) = dipstr*weights(i)
         call tri_bak(w,vertout,zparts(1,i))
         dipvec(1,i) = dipvec(1,1)
         dipvec(2,i) = dipvec(2,1)
         dipvec(3,i) = dipvec(3,1)
      enddo
      call l3dformta_dp(ier,scale,zparts,ctemp,dipvec,nnodes,x0y0z0,
     1               nterms,local)
      return
      end
c
c
c
c
c
      subroutine l3dformtatrid(ier,scale,triang,trinorm,dipstr,
     1          ntri,x0y0z0,norder,nterms,local,mptemp)
      implicit none
c
c
c     This subroutine INCREMENTS the local expansion about x0y0z0 
c     to include contributions from DLP on collection of triangles.
c
c     INPUT:
c
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     dipstr             = array of DLP densities
c     trinorm            = normal to triangle
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of multipole expansion
c     mptemp             = work array to hold temp multipole expansion
c
c     OUTPUT:
c
c     local              =  coefficients of local expansion 
c                           are INCREMENTED by contributions from each
c                           DLP triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c-----------------------------------------------------------------------
c
      integer  nterms, ntri, ier, lused
      integer  i, norder
      real *8  triang(3,3,1), trinorm(3,1), x0y0z0(3), scale
      complex *16  eye, zk, dipstr(1)
      complex *16  local(0:nterms,-nterms:nterms), mptemp(1)
c
      do i = 1,ntri
         call l3dformtatridone(ier,scale,triang(1,1,i),trinorm(1,i),
     1          dipstr(i),x0y0z0, norder,nterms,mptemp)
         call l3dadd(mptemp,local,nterms)
      enddo
      return
      end
c
c
      subroutine l3dformtatrid_add(ier,scale,triang,trinorm,dipstr,
     1          ntri,x0y0z0,norder,nterms,local)
      implicit none
c
c
c     This subroutine INCREMENTS the local expansion about x0y0z0 
c     to include contributions from DLP on collection of triangles.
c
c     INPUT:
c
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     dipstr             = array of DLP densities
c     trinorm            = normal to triangle
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of multipole expansion
c     mptemp             = work array to hold temp multipole expansion
c
c     OUTPUT:
c
c     local              =  coefficients of local expansion 
c                           are INCREMENTED by contributions from each
c                           DLP triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c-----------------------------------------------------------------------
c
      integer  nterms, ntri, ier, lused
      integer  i, norder
      real *8  triang(3,3,1), trinorm(3,1), x0y0z0(3), scale
      complex *16  eye, zk, dipstr(1)
      complex *16  local(0:nterms,-nterms:nterms)
      complex *16, allocatable :: mptemp(:,:)
c
      allocate( mptemp(0:nterms,-nterms:nterms) )
c
      do i = 1,ntri
         call l3dformtatridone(ier,scale,triang(1,1,i),trinorm(1,i),
     1          dipstr(i),x0y0z0, norder,nterms,mptemp)
         call l3dadd(mptemp,local,nterms)
      enddo
      return
      end
c
c
c***********************************************************************
c
c    Low level FMM routines for triangles in Laplace regime.
c    (forming multipole expansions and forming taylor expansions).
c       
c    Optimized routines, use precomputed recursion coefficients for 
c    Legendre functions...
c
c
      subroutine l3dformmptrisone2(ier, scale, triang,
     1          charge, x0y0z0, norder, nterms, mpole, wlege,nlege)
      implicit real *8 (a-h,o-z)
c
c     Form  multipole expansion due to SLP on  single triangle.
c
c     INPUT:
c
c     scale        = scaling parameter
c     triang(i,j)  = ith coord of jth vertex
c     charge       = SLP strength
c     x0y0z0       = center of the expansion
c     norder       = order of Gaussian rule on triangle
c     nterms       = order of multipole expansion
c
c     OUTPUT:
c
c     mpole        =  induced multipole expansion 
c     lused        =  amount of workspace w used
c     ier          =  error return code
c
c---------------------------------------------------------------------------
c
      integer  ier, ifinit, norder, nnodes, nterms, lused
      integer  i, j, l, m
c
      real *8  triang(3,3), x0y0z0(3)
      real *8  scale
      real *8  rnodes(2,500), weights(500), zparts(3,500)
      real *8  vert1(2),vert2(2),vert3(2),vertout(3)
      real *8  w(20), wlege(1)
c
      complex *16 eye, charge, ctemp(500)
      complex *16  mpole(0:nterms,-nterms:nterms)
c
      data eye/(0.0d0,1.0d0)/
c
      call tri_ini(triang(1,1),triang(1,2),triang(1,3),w,
     1             vert1,vert2,vert3)
c
      call triasymq(norder,vert1,vert2,vert3,rnodes,weights,nnodes)
ccc      call prin2(' weights are *',weights,nnodes)
c
      do i = 1,nnodes
         vertout(1) = rnodes(1,i)
         vertout(2) = rnodes(2,i)
         vertout(3) = 0
         call tri_bak(w,vertout,zparts(1,i))
         ctemp(i) = charge*weights(i)
      enddo
ccc      call prin2(' ctemp are *',ctemp,2*nnodes)
ccc      call prin2(' zparts are *',zparts,3*nnodes)
ccc      call prin2(' zk is *',zk,2)
ccc      call prin2(' scale is *',scale,1)
ccc      call prin2(' x0y0z0 is *',x0y0z0,3)
ccc      call prinf(' nterms is *',nterms,1)
      call l3dformmp_trunc(ier,scale,zparts,ctemp,nnodes,x0y0z0,
     1               nterms,nterms,mpole,wlege,nlege)
      return
      end
c
c
c
c
c
      subroutine l3dformmptris2_add(ier, scale, triang, charge, 
     1          ntri, x0y0z0, norder, nterms, mpole, wlege,nlege)
      implicit real *8 (a-h,o-z)
c
c
c     This subroutine INCREMENTS the multipole expansion to include
c     contributions from SLP on collection of triangles.
c
c     INPUT:
c
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     charge             = array of SLP densities
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of multipole expansion
c     mptemp             = work array to hold temp multipole expansion
c
c     OUTPUT:
c
c     mpole              =  coefficients of multipole expansion 
c                           are INCREMENTED by contributions from each
c                           triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c---------------------------------------------------------------------------
c
      integer  nterms, ntri, ier, lused
      integer  i, norder
      real *8  triang(3,3,1), x0y0z0(3), scale, wlege(1)
      complex *16  zk, charge(1)
      complex *16  mpole(0:nterms,-nterms:nterms)
      complex *16, allocatable :: mptemp(:,:)
c
      allocate( mptemp(0:nterms,-nterms:nterms) )
c
      do i = 1,ntri
         call l3dformmptrisone2(ier,scale,triang(1,1,i),charge(i), 
     1          x0y0z0, norder, nterms, mptemp, wlege,nlege)
         call l3dadd(mptemp,mpole,nterms)
      enddo
      return
      end
c
c
c
c
c
c
      subroutine l3dformmptridone2(ier, scale, triang, trinorm,
     1          dipstr, x0y0z0, norder, nterms, mpole, wlege, nlege)
      implicit real *8 (a-h,o-z)
c
c     Form multipole expansion due to DLP on single triangle.
c
c     INPUT:
c
c     scale        = scaling parameter
c     triang(i,j)  = ith coord of jth vertex
c     trinorm(i)   = ith coord of normal
c     dipstr       = DLP strength
c     x0y0z0       = center of the expansion
c     norder       = order of Gaussian rule on triangle
c     nterms       = order of multipole expansion
c
c     OUTPUT:
c
c     mpole        =  induced multipole expansion 
c     lused        =  amount of workspace w used
c     ier          =  error return code
c
c---------------------------------------------------------------------------
c
      integer  ier, ifinit, norder, nnodes, nterms, lused
      integer  i, j, l, m
c
      real *8  triang(3,3), trinorm(3), x0y0z0(3)
      real *8  scale
      real *8  rnodes(2,500), weights(500), zparts(3,500)
      real *8  dipvec(3,500)
      real *8  vert1(2),vert2(2),vert3(2),vertout(3)
      real *8  w(20), wlege(1)
c
      complex *16 eye, dipstr, zk, ctemp(500)
      complex *16  mpole(0:nterms,-nterms:nterms)
c
      data eye/(0.0d0,1.0d0)/
c
      call tri_ini(triang(1,1),triang(1,2),triang(1,3),w,
     1             vert1,vert2,vert3)
c
      call triasymq(norder,vert1,vert2,vert3,rnodes,weights,nnodes)
ccc        call prinf('nnodes=*',nnodes,1)
c
        dipvec(1,1)=trinorm(1)
        dipvec(2,1)=trinorm(2)
        dipvec(3,1)=trinorm(3)
c
      do i = 1,nnodes
         vertout(1) = rnodes(1,i)
         vertout(2) = rnodes(2,i)
         vertout(3) = 0
         ctemp(i) = dipstr*weights(i)
         call tri_bak(w,vertout,zparts(1,i))
         dipvec(1,i) = dipvec(1,1)
         dipvec(2,i) = dipvec(2,1)
         dipvec(3,i) = dipvec(3,1)
      enddo
      call l3dformmp_dp_trunc
     $   (ier,scale,zparts,ctemp,dipvec,nnodes,x0y0z0,
     1               nterms,nterms,mpole,wlege,nlege)
      return
      end
c
c
c
c
c
      subroutine l3dformmptrid2_add(ier,scale,triang,trinorm,dipstr,
     1          ntri,x0y0z0,norder,nterms,mpole,wlege,nlege)
      implicit real *8 (a-h,o-z)
c
c
c     This subroutine INCREMENTS the multipole expansion about x0y0z0 
c     to include contributions from DLP on collection of triangles.
c
c     INPUT:
c
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     dipstr             = array of DLP densities
c     trinorm            = normal to triangle
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of multipole expansion
c     mptemp             = work array to hold temp multipole expansion
c
c     OUTPUT:
c
c     mpole              =  coefficients of multipole expansion 
c                           are INCREMENTED by contributions from each
c                           DLP triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c---------------------------------------------------------------------------
c
      integer  nterms, ntri, ier, lused
      integer  i, norder
      real *8  triang(3,3,1), trinorm(3,1), x0y0z0(3), scale, wlege(1)
      complex *16  eye, zk, dipstr(1)
      complex *16  mpole(0:nterms,-nterms:nterms)
      complex *16, allocatable :: mptemp(:,:)
c
      allocate( mptemp(0:nterms,-nterms:nterms) )
c
      do i = 1,ntri
         call l3dformmptridone2(ier,scale,triang(1,1,i),trinorm(1,i),
     1          dipstr(i),x0y0z0, norder,nterms,mptemp,wlege,nlege)
         call l3dadd(mptemp,mpole,nterms)
      enddo
c
      return
      end
c
c
c
c     FORMTA routines: Individual triangle to local expansion
c     for list 3/4 processing.
c
c
c
      subroutine l3dformtatrisone2(ier, scale, triang,
     1          charge, x0y0z0, norder, nterms, local, wlege, nlege)
      implicit real *8 (a-h,o-z)
c
c     Form local expansion about x0y0z0 due to SLP on single triangle.
c
c     INPUT:
c
c     scale        = scaling parameter
c     triang(i,j)  = ith coord of jth vertex
c     charge       = SLP strength
c     x0y0z0       = center of the expansion
c     norder       = order of Gaussian rule on triangle
c     nterms       = order of multipole expansion
c
c     OUTPUT:
c
c     local        =  induced local expansion 
c     lused        =  amount of workspace w used
c     ier          =  error return code
c
c---------------------------------------------------------------------------
c
      integer  ier, ifinit, norder, nnodes, nterms, lused
      integer  i, j, l, m
c
      real *8  triang(3,3), x0y0z0(3)
      real *8  scale
      real *8  rnodes(2,500), weights(500), zparts(3,500)
      real *8  vert1(2),vert2(2),vert3(2),vertout(3)
      real *8  w(20), wlege(1)
c
      complex *16 eye, charge, zk, ctemp(500)
      complex *16  local(0:nterms,-nterms:nterms)
c
      data eye/(0.0d0,1.0d0)/
c
      call tri_ini(triang(1,1),triang(1,2),triang(1,3),w,
     1             vert1,vert2,vert3)
c
      call triasymq(norder,vert1,vert2,vert3,rnodes,weights,nnodes)
ccc      call prin2(' weights are *',weights,nnodes)
c
      do i = 1,nnodes
         vertout(1) = rnodes(1,i)
         vertout(2) = rnodes(2,i)
         vertout(3) = 0
         call tri_bak(w,vertout,zparts(1,i))
         ctemp(i) = charge*weights(i)
      enddo
ccc      call prin2(' ctemp are *',ctemp,2*nnodes)
ccc      call prin2(' zparts are *',zparts,3*nnodes)
ccc      call prin2(' zk is *',zk,2)
ccc      call prin2(' scale is *',scale,1)
ccc      call prin2(' x0y0z0 is *',x0y0z0,3)
ccc      call prinf(' nterms is *',nterms,1)
      call l3dformta_trunc(ier,scale,zparts,ctemp,nnodes,x0y0z0,
     1               nterms,nterms,local,wlege,nlege)
      return
      end
c
c
c
c
c
      subroutine l3dformtatris2_add(ier, scale, triang, charge, 
     1          ntri, x0y0z0, norder, nterms, local, wlege, nlege)
      implicit real *8 (a-h,o-z)
c
c
c     This subroutine INCREMENTS the local expansion to include
c     contributions from SLP on collection of triangles.
c
c     INPUT:
c
c     zk                 = Helmholtz parameter
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     charge             = array of SLP densities
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of multipole expansion
c     mptemp             = work array to hold temp multipole expansion
c
c     OUTPUT:
c
c     local              =  coefficients of local expansion 
c                           are INCREMENTED by contributions from each
c                           SLP triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c---------------------------------------------------------------------------
c
      integer  nterms, ntri, ier, lused
      integer  i, norder
      real *8  triang(3,3,1), x0y0z0(3), scale, wlege(1)
      complex *16  zk, charge(1)
      complex *16  local(0:nterms,-nterms:nterms)
      complex *16, allocatable :: mptemp(:,:)
c
      allocate( mptemp(0:nterms,-nterms:nterms) )
c
      do i = 1,ntri
         call l3dformtatrisone2(ier,scale,triang(1,1,i),charge(i), 
     1          x0y0z0, norder, nterms, mptemp, wlege, nlege)
         call l3dadd(mptemp,local,nterms)
      enddo
      return
      end
c
c
c
c
c
c
      subroutine l3dformtatridone2(ier, scale, triang, trinorm,
     1          dipstr, x0y0z0, norder, nterms, local,wlege,nlege)
      implicit real *8 (a-h,o-z)
c
c     Form local expansion about x0y0z0 due to DLP on single triangle.
c
c     INPUT:
c
c     scale        = scaling parameter
c     triang(i,j)  = ith coord of jth vertex
c     trinorm(i)   = ith coord of normal
c     dipstr       = DLP strength
c     x0y0z0       = center of the expansion
c     norder       = order of Gaussian rule on triangle
c     nterms       = order of multipole expansion
c
c     OUTPUT:
c
c     local        =  induced local expansion 
c     lused        =  amount of workspace w used
c     ier          =  error return code
c
c---------------------------------------------------------------------------
c
      integer  ier, ifinit, norder, nnodes, nterms, lused
      integer  i, j, l, m
c
      real *8  triang(3,3), trinorm(3), x0y0z0(3)
      real *8  scale
      real *8  rnodes(2,500), weights(500), zparts(3,500)
      real *8  dipvec(3,500)
      real *8  vert1(2),vert2(2),vert3(2),vertout(3)
      real *8  w(20), wlege(1)
c
      complex *16 eye, dipstr, zk, ctemp(500)
      complex *16  local(0:nterms,-nterms:nterms)
c
      data eye/(0.0d0,1.0d0)/
c
      call tri_ini(triang(1,1),triang(1,2),triang(1,3),w,
     1             vert1,vert2,vert3)
c
      call triasymq(norder,vert1,vert2,vert3,rnodes,weights,nnodes)
c
        dipvec(1,1)=trinorm(1)
        dipvec(2,1)=trinorm(2)
        dipvec(3,1)=trinorm(3)
c
      do i = 1,nnodes
         vertout(1) = rnodes(1,i)
         vertout(2) = rnodes(2,i)
         vertout(3) = 0
         ctemp(i) = dipstr*weights(i)
         call tri_bak(w,vertout,zparts(1,i))
         dipvec(1,i) = dipvec(1,1)
         dipvec(2,i) = dipvec(2,1)
         dipvec(3,i) = dipvec(3,1)
      enddo
      call l3dformta_dp_trunc
     $   (ier,scale,zparts,ctemp,dipvec,nnodes,x0y0z0,
     1               nterms,nterms,local,wlege,nlege)
      return
      end
c
c
c
c
c
      subroutine l3dformtatrid2_add(ier,scale,triang,trinorm,dipstr,
     1          ntri,x0y0z0,norder,nterms,local,wlege,nlege)
      implicit real *8 (a-h,o-z)
c
c
c     This subroutine INCREMENTS the local expansion about x0y0z0 
c     to include contributions from DLP on collection of triangles.
c
c     INPUT:
c
c     scale              = scaling parameter
c     triang(i,j,k)      = ith coord of jth vertex of kth triangle
c     dipstr             = array of DLP densities
c     trinorm            = normal to triangle
c     ntri               = number of triangles
c     x0y0z0             = center of the expansion
c     norder             = order of Gaussian rule used 
c     nterms             = order of multipole expansion
c     mptemp             = work array to hold temp multipole expansion
c
c     OUTPUT:
c
c     local              =  coefficients of local expansion 
c                           are INCREMENTED by contributions from each
c                           DLP triangle.
c     lused              =  amount of workspace w used
c     ier                =  error return code
c
c-----------------------------------------------------------------------
c
      integer  nterms, ntri, ier, lused
      integer  i, norder
      real *8  triang(3,3,1), trinorm(3,1), x0y0z0(3), scal, wlege(1)
      complex *16  eye, zk, dipstr(1)
      complex *16  local(0:nterms,-nterms:nterms)
      complex *16, allocatable :: mptemp(:,:)
c
      allocate( mptemp(0:nterms,-nterms:nterms) )
c
      do i = 1,ntri
         call l3dformtatridone2(ier,scale,triang(1,1,i),trinorm(1,i),
     1          dipstr(i),x0y0z0, norder,nterms,mptemp,wlege,nlege)
         call l3dadd(mptemp,local,nterms)
      enddo
      return
      end
c
c

