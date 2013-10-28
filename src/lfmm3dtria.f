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
c     $Date: 2010-07-10 18:32:21 -0400 (Sat, 10 Jul 2010) $
c     $Revision: 1078 $
c
c       
c     This file contains the main FMM routines and some related
c     subroutines for evaluating Laplace potentials and fields on
c     surfaces defined by a collection of positively oriented flat
c     triangles.  (FORTRAN 90 VERSION)
c
c     lfmm3dtria - Laplace FMM in R^3: evaluate all pairwise triangle
c         interactions (including self-interaction)
c
c     lfmm3dtriaself - Laplace FMM in R^3: evaluate all pairwise triangle
c         interactions (including self-interaction)
c
c     lfmm3dtriatarg - Laplace FMM in R^3: evaluate all pairwise
c         triangle interactions (including self-interaction) +
c         interactions with targets
c
c     lfmm3dtriampftarg - post-processing routine: evaluate 
c         interactions with targets (off-surface)
c
c     l3dtriadirect - Laplace interactions in R^3: evaluate all pairwise
c         triangle interactions (including self-interaction) +
c         interactions with targets via direct O(N^2) algorithm
c
c
c     Note: dipole vectors MUST BE SET EQUAL to the triangle normals.
c
c
c***********************************************************************
        subroutine lfmm3dtria(ier,iprec,nsource,
     $     triaflat,trianorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld)
c***********************************************************************
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the harmonic potential and 
c       field due to a collection of flat triangles with constant
c       single and/or double layer densities. We use (1/r) for the Green's 
c       function, without the (1/4 pi) scaling.
c       Self-interactions are included.
c   
c       The main FMM routine permits both evaluation on surface
c       and at a collection of off-surface targets. 
c       This subroutine is used to simplify the user interface 
c       (by setting the number of targets to zero) and calling the more 
c       general FMM.
c
c       See below for explanation of calling sequence arguments.
c  
        lused7=0
c
        ntarget=0
        ifpottarg=0
        iffldtarg=0
c
        call lfmm3dtriatarg(ier,iprec,nsource,
     $     triaflat,trianorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,
     $     ntarget,target,ifpottarg,pottarg,iffldtarg,fldtarg)
c
        return
        end
c
c
c
c
c
c***********************************************************************
        subroutine lfmm3dtriaself(ier,iprec,nsource,
     $     triaflat,trianorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld)
c***********************************************************************
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the harmonic potential and 
c       field due to a collection of flat triangles with constant
c       single and/or double layer densities. We use (1/r) for the Green's 
c       function, without the (1/4 pi) scaling.
c       Self-interactions are included.
c   
c       The main FMM routine permits both evaluation on surface
c       and at a collection of off-surface targets. 
c       This subroutine is used to simplify the user interface 
c       (by setting the number of targets to zero) and calling the more 
c       general FMM.
c
c       See below for explanation of calling sequence arguments.
c  
        lused7=0
c
        ntarget=0
        ifpottarg=0
        iffldtarg=0
c
        call lfmm3dtriatarg(ier,iprec,nsource,
     $     triaflat,trianorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,
     $     ntarget,target,ifpottarg,pottarg,iffldtarg,fldtarg)
c
        return
        end
c
c
c
c
c
c***********************************************************************
        subroutine lfmm3dtriatarg(ier,iprec,nsource,
     $     triaflat,trianorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,
     $     ntarget,target,ifpottarg,pottarg,iffldtarg,fldtarg)
c***********************************************************************
        implicit real *8 (a-h,o-z)
c       
c       Laplace FMM in R^3: evaluate all pairwise triangle
c       interactions (including self-interaction) + interactions with targets
c
c       This is the principal subroutine for evaluating 
c       harmonic layer potentials on (flat) triangulated surfaces.
c       It permits the evaluation of a single layer potential
c       with piecewise constant density defined by <<charge>>
c       and a dipole layer potential with piecewise constant density 
c       and dipole orientation defined by <<dipstr,dipvec>>.
c
c       We use (1/r) for the Green's function,
c       without the (1/4 pi) scaling.  Self-interactions are included.
c   
c       It is capable of evaluating the layer potentials either on 
c       or off the surface (or both).            
c
c       This is primarily a memory management code. 
c       The actual work is carried out in subroutine lfmm3dtriatargmain.
c
c       INPUT PARAMETERS:
c
c       iprec:  FMM precision flag
c
c                 -2 => tolerance =.5d0
c                 -1 => tolerance =.5d-1
c                  0 => tolerance =.5d-2
c                  1 => tolerance =.5d-3
c                  2 => tolerance =.5d-6
c                  3 => tolerance =.5d-9
c                  4 => tolerance =.5d-12
c                  5 => tolerance =.5d-15
c
c       nsource: integer:  number of triangles
c       triaflat: real *8 (3,3,nsource): triangle coordinate array
c       trianorm: real *8 (3,nsource): triangle normals
c       source: real *8 (3,nsource):  triangle centroids
c       ifcharge:  single layer potential (SLP) flag
c                  ifcharge = 1   =>  include SLP contribution
c                                     otherwise do not
c       charge: complex *16 (nsource): piecewise constant SLP strength
c       ifdipole:  dipole layer potential (DLP) flag
c                  ifdipole = 1   =>  include DLP contribution
c                                     otherwise do not
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
c       ier   =  error return code
c                ier=0     =>  normal execution
c                ier=4     =>  cannot allocate tree workspace
c                ier=8     =>  cannot allocate bulk FMM  workspace
c                ier=16    =>  cannot allocate mpole expansion
c                              workspace in FMM
c
c       pot: complex *16 (nsource): potential at triangle centroids
c       fld: complex *16 (3,nsource): field (-gradient) at triangle centroids 
c       pottarg: complex *16 (ntarget): potential at target locations 
c       fldtarg: complex *16 (3,ntarget): field (-gradient) at target locations 
c
c
cf2py   intent(out) ier
cf2py   intent(in) iprec
cf2py   intent(in) nsource
cf2py   intent(in) triaflat,trianorm,source
cf2py   intent(in) ifcharge,charge
cf2py   check(!ifcharge || (shape(charge,0) == nsource))  charge
cf2py   depend(nsource)  charge
cf2py   intent(in) ifdipole,dipvec,dipstr
cf2py   check(!ifdipole || (shape(dipstr,0) == nsource))  dipstr
cf2py   depend(nsource)  dipstr
cf2py   intent(in) ifpot,iffld
cf2py   intent(out) pot,fld
cf2py   intent(in) ifpottarg
cf2py   intent(in) iffldtarg
cf2py   intent(in) target
cf2py   intent(in) ntarget
cf2py   check((!ifpottarg && !iffldtarg) || (shape(target,0)==3 && shape(target,1) == ntarget))  target
cf2py   check((!ifpottarg) || (shape(pottarg,0)==ntarget))  pottarg
cf2py   check((!iffldtarg) || (shape(fldtarg,0)==3 && shape(fldtarg,1) == ntarget))  fldtarg
c
c       (F2PY workaround: pottarg, fldtarg must be input because f2py
c       refuses to allocate zero-size output arrays.)
c
cf2py   intent(in,out) pottarg,fldtarg

        dimension triaflat(3,3,nsource)
        dimension trianorm(3,nsource)
        dimension source(3,nsource)
        dimension dipvec(3,nsource)
        dimension target(3,ntarget)
        complex *16 charge(nsource)
        complex *16 dipstr(nsource)
        complex *16 ima
        complex *16 pot(nsource)
        complex *16 fld(3,nsource)
        complex *16 pottarg(ntarget)
        complex *16 fldtarg(3,ntarget)
c
        dimension timeinfo(10)
c       
c
c     Note: various arrays dimensioned here to 200.
c     That allows for 200 evels of refinment, which is 
c     more than enough for any non-pathological case.
c
c       
        dimension laddr(2,200)
        dimension nterms(0:200)
        integer box(20)
        integer box1(20)
        dimension scale(0:200)
        dimension bsize(0:200)
        dimension center(3)
        dimension center0(3),corners0(3,8)
        dimension center1(3),corners1(3,8)
        real *8, allocatable :: w(:) 
        real *8, allocatable :: wlists(:) 
        real *8, allocatable :: wrmlexp(:) 
        complex *16 ptemp,ftemp(3)
c       
        data ima/(0.0d0,1.0d0)/
c       
        ier=0
        lused7=0
c       
        done=1
        pi=4*atan(done)
c
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     
      ifprint=0
c
c
c     set fmm tolerance based on iprec flag.
c       
        if( iprec .eq. -2 ) epsfmm=.5d-0 
        if( iprec .eq. -1 ) epsfmm=.5d-1
        if( iprec .eq. 0 ) epsfmm=.5d-2
        if( iprec .eq. 1 ) epsfmm=.5d-3
        if( iprec .eq. 2 ) epsfmm=.5d-6
        if( iprec .eq. 3 ) epsfmm=.5d-9
        if( iprec .eq. 4 ) epsfmm=.5d-12
        if( iprec .eq. 5 ) epsfmm=.5d-15
        if( iprec .eq. 6 ) epsfmm=0
c       
        if (ifprint.eq.1) call prin2('epsfmm=*',epsfmm,1)
c
c     set criterion for box subdivision (number of sources per box)
c
        if( iprec .eq. -2 ) nbox=8/3
        if( iprec .eq. -1 ) nbox=15/3
        if( iprec .eq. 0 ) nbox=30/3
        if( iprec .eq. 1 ) nbox=60/3
        if( iprec .eq. 2 ) nbox=120/3
        if( iprec .eq. 3 ) nbox=240/3
        if( iprec .eq. 4 ) nbox=480/3
        if( iprec .eq. 5 ) nbox=700/3
        if( iprec .eq. 6 ) nbox=nsource+ntarget
c
        if (ifprint.eq.1) call prinf('nbox=*',nbox,1)
c
c     create oct-tree data structure
c
        ntot = 100*(nsource+ntarget)+10000
        if( iprec .eq. -2 ) ntot = ntot * 1.5*1.5*1.5
        if( iprec .eq. -1 ) ntot = ntot * 1.5*1.5
        do ii = 1,10
           allocate (wlists(ntot))
           call lfmm3dparttree(ier,iprec,
     $        nsource,source,ntarget,target,
     $        nbox,epsfmm,iisource,iitarget,iwlists,lwlists,
     $        nboxes,laddr,nlev,center,size,
     $        wlists,ntot,lused7)
           if (ier.ne.0) then
              deallocate(wlists)
              ntot = ntot*1.5
              call prinf(' increasing allocation, ntot is *',ntot,1)
           else
              goto 1200
           endif
        enddo
1200    continue
        if (ier.ne.0) then
           call prinf(' exceeded max allocation, ntot is *',ntot,1)
           ier = 4
           return          
        endif

c
c     lused7 is counter that steps through workspace,
c     keeping track of total memory used.
c
        lused7=1
        do i = 0,nlev
        scale(i) = 1.0d0
        enddo
c       
        if (ifprint.eq.1) call prin2('scale=*',scale,nlev+1)
c
c       carve up workspace further
c
c     itriaflatsort is pointer for sorted triangle coordinates
c     itrianormsort is pointer for sorted triangle normals
c     isourcesort is pointer for sorted triangle centroids
c     itargetsort is pointer for sorted target locations
c     ichargesort is pointer for sorted charge densities
c     idipvecsort is pointer for sorted dipole orientation vectors
c     idipstrsort is pointer for sorted dipole densities
c
c
        itriaflatsort = lused7 + 5
        ltriaflatsort = 3*3*nsource
        itrianormsort = itriaflatsort + ltriaflatsort
        ltrianormsort = 3*nsource
        isourcesort = itrianormsort + ltrianormsort 
        lsourcesort = 3*nsource
        itargetsort = isourcesort+lsourcesort
        ltargetsort = 3*ntarget
        ichargesort = itargetsort+ltargetsort
        lchargesort = 2*nsource
        idipvecsort = ichargesort+lchargesort
        if (ifdipole.eq.1) then
          ldipvec = 3*nsource
          ldipstr = 2*nsource
        else
          ldipvec = 3
          ldipstr = 2
        endif
        idipstrsort = idipvecsort + ldipvec
        lused7 = idipstrsort + ldipstr       
c
c     allocate workspace for potential and field arrays
c
        ipot = lused7
        lpot = 2*nsource
        lused7=lused7+lpot
c       
        ifld = lused7
        if( iffld .eq. 1) then
           lfld = 2*(3*nsource)
        else
           lfld=6
        endif
        lused7=lused7+lfld
c      
        ipottarg = lused7
        lpottarg = 2*ntarget
        lused7=lused7+lpottarg
c       
        ifldtarg = lused7
        if( iffldtarg .eq. 1) then
           lfldtarg = 2*(3*ntarget)
        else
           lfldtarg=6
        endif
        lused7=lused7+lfldtarg
c      
        if (ifprint.eq.1) call prinf(' lused7 is *',lused7,1)
c
c       based on FMM tolerance, compute expansion lengths nterms(i)
c      
        nmax = 0
        do i = 0,nlev
           bsize(i)=size/2.0d0**i
           call l3dterms(epsfmm, nterms(i), ier)
           if (nterms(i).gt. nmax .and. i.ge. 2) nmax = nterms(i)
        enddo
c
        nquad=2*nmax        
c       
c     ixnodes is pointer for quadrature nodes
c     iwhts is pointer for quadrature weights
c
        ixnodes = lused7 
        iwts = ixnodes + nquad
        lused7 = iwts + nquad
c
        if (ifprint.eq.1) call prinf('nterms=*',nterms,nlev+1)
        if (ifprint.eq.1) call prinf('nmax=*',nmax,1)
c
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(2,nboxes).
c
c     iiaddr is pointer to iaddr array, itself contained in workspace.
c     imptemp is pointer for single expansion (dimensioned by nmax)
c   
c       ... allocate iaddr and temporary arrays
c
        iiaddr = lused7 
        imptemp = iiaddr + 2*nboxes
        lmptemp = (nmax+1)*(2*nmax+1)*2
        lused7 = imptemp + lmptemp
        allocate(w(lused7),stat=ier)
        if (ier.ne.0) then
           call prinf(' cannot allocate bulk FMM workspace,
     1                   lused7 is *',lused7,1)
           ier = 8
           return          
        endif
c
c     reorder triangles, centroids etc. so that each box holds
c     contiguous list of source numbers.
c
        call l3dreorder(nsource,source,ifcharge,charge,wlists(iisource),
     $     ifdipole,dipstr,dipvec,
     1     w(isourcesort),w(ichargesort),w(idipvecsort),w(idipstrsort)) 
c       
        call l3dreordertria(nsource,wlists(iisource),
     $     triaflat,w(itriaflatsort),trianorm,w(itrianormsort))
c
        call l3dreordertarg(ntarget,target,wlists(iitarget),
     $     w(itargetsort))
c
      if (ifprint.eq.1) then
        call prinf('finished reordering=*',ier,1)
        call prinf('ier=*',ier,1)
        call prinf('nboxes=*',nboxes,1)
        call prinf('nlev=*',nlev,1)
        call prinf('nboxes=*',nboxes,1)
        call prinf('lused7=*',lused7,1)
      endif
c       
        ifinit=1
        call legewhts(nquad,w(ixnodes),w(iwts),ifinit)
c
ccc        call prin2('xnodes=*',xnodes,nquad)
ccc        call prin2('wts=*',wts,nquad)
c
c     allocate memory need by multipole, local expansions at all
c     levels
c     irmlexp is pointer for workspace need by various fmm routines,
c
        call l3dmpalloc(wlists(iwlists),w(iiaddr),nboxes,lmptot,nterms)
c
        if (ifprint.eq.1) call prinf(' lmptot is *',lmptot,1)
c       
        irmlexp = 1
        lused7 = irmlexp + lmptot 
        if (ifprint.eq.1) call prinf(' lused7 is *',lused7,1)
        allocate(wrmlexp(lused7))
        if (ier.ne.0) then
           call prinf(' cannot allocate mpole expansion workspace,
     1                   lused7 is *',lused7,1)
           ier = 16
           return          
        endif
c
c     Memory allocation is complete. 
c     Call main fmm routine. There are, unfortunately, a lot
c     of parameters here. ifevalfar and ifevalloc determine
c     whether far field and local fields (respectively) are to 
c     be evaluated. Setting both to 1 means that both will be
c     computed (which is the normal scenario).
c
        ifevalfar=1
        ifevalloc=1
c
        call lfmm3dtriatargmain(ier,iprec,
     $     ifevalfar,ifevalloc,
     $     nsource,w(itriaflatsort),w(itrianormsort),
     $     w(isourcesort),wlists(iisource),
     $     ifcharge,w(ichargesort),
     $     ifdipole,w(idipstrsort),w(idipvecsort),
     $     ifpot,w(ipot),iffld,w(ifld),
     $     ntarget,w(itargetsort),wlists(iitarget),
     $     ifpottarg,w(ipottarg),iffldtarg,w(ifldtarg),
     $     epsfmm,w(iiaddr),wrmlexp(irmlexp),w(imptemp),lmptemp,
     $     w(ixnodes),w(iwts),nquad,
     $     nboxes,laddr,nlev,scale,bsize,nterms,
     $     wlists(iwlists),lwlists)
c
c       parameter ier from targmain routine is currently meaningless, reset to 0
        if( ier .ne. 0 ) ier = 0
c
        if (ifprint.eq.1) call prinf('lwlists=*',lwlists,1)
        if (ifprint.eq.1) call prinf('lused total =*',lused7,1)       
c
        if (ifprint.eq.1) 
     $     call prin2('memory / point = *',(lused7)/dble(nsource),1)
c       
ccc        call prin2('after w=*', w(1+lused7-100), 2*100)
c
        if(ifpot .eq. 1) 
     $     call l3dpsort(nsource,wlists(iisource),w(ipot),pot)
        if(iffld .eq. 1) 
     $     call l3dfsort(nsource,wlists(iisource),w(ifld),fld)
c
        if(ifpottarg .eq. 1 )
     $     call l3dpsort(ntarget,wlists(iitarget),w(ipottarg),pottarg)
        if(iffldtarg .eq. 1) 
     $     call l3dfsort(ntarget,wlists(iitarget),w(ifldtarg),fldtarg)
c       
        return
        end
c
c
c
c
c
        subroutine lfmm3dtriatargmain(ier,iprec,
     $     ifevalfar,ifevalloc,
     $     nsource,triaflatsort,trianormsort,sourcesort,isource,
     $     ifcharge,chargesort,
     $     ifdipole,dipstrsort,dipvecsort,
     $     ifpot,pot,iffld,fld,ntarget,
     $     targetsort,itarget,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     epsfmm,iaddr,rmlexp,
     $     mptemp,lmptemp,xnodes,wts,nquad,
     $     nboxes,laddr,nlev,scale,bsize,nterms,
     $     wlists,lwlists)
        implicit real *8 (a-h,o-z)
        dimension triaflatsort(3,3,1),trianormsort(3,1)
        dimension sourcesort(3,1), isource(1)
        complex *16 chargesort(1)
        complex *16 dipstrsort(1)
        dimension dipvecsort(3,1)
        complex *16 ima
        complex *16 pot(1)
        complex *16 fld(3,1)
        dimension targetsort(3,1), itarget(1)
        complex *16 pottarg(1)
        complex *16 fldtarg(3,1)
        dimension wlists(1)
        dimension iaddr(2,nboxes)
        real *8 rmlexp(1)
        complex *16 mptemp(lmptemp)
        dimension xnodes(nquad),wts(nquad)
        dimension timeinfo(10)
        dimension center(3)
        dimension laddr(2,200)
        dimension scale(0:200)
        dimension bsize(0:200)
        dimension nterms(0:200)
        dimension list(10 000)
        complex *16 ptemp,ftemp(3)
        integer box(20)
        dimension center0(3),corners0(3,8)
        integer box1(20)
        dimension center1(3),corners1(3,8)
        dimension itable(-3:3,-3:3,-3:3)
        dimension wlege(40000)
        dimension nterms_eval(4,0:200)
c
        data ima/(0.0d0,1.0d0)/

c
c     INPUT PARAMETERS:
c
c     iprec        precision flag (see above)   ELIMINATE???
c     ifevalfar    far field flag (1 means compute far field, 
c                                  else dont)
c     ifevalloc    local field flag (1 means compute local field, 
c                                    else dont)
c     nsource      number of triangles
c     triaflatsort sorted triangle coordinate array 
c     trianormsort sorted triangle normals array 
c     sourcesort   sorted triangle centroids
c     isource      sorting index for sources
c     ifcharge     flag indicating potential includes contribution
c                  from single layer potential 
c     chargesort   sorted piecewise constant single layer densities
c     ifdipole     flag indicating potential includes contribution
c                  from dipole layer potential 
c     dipstrsort   sorted piecewise constant double layer density
c     dipvecsort   sorted dipole orientation vectors
c     ifpot        potential flag (1 => compute, else do not)
c     iffld        field flag (1 => compute, else do not)
c     ntarget      number of targets
c     targetsort   sorted array of target locations
c     itarget      sorting index for targets
c     ifpottarg    target potential flag (1 => compute, else do not)
c     iffldtarg    target field flag (1 => compute, else do not)
c     epsfmm       FMM tolerance
c     iaddr        iaddr(2,nboxes) array points to mpole/local
c                     expansions for each box
c     rmlexp       workspace to contain mpole/local expansions.
c     mptemp       workspace for single expansion 
c     xnodes       workspace to hold quadrature nodes
c     wts          workspace to hold quadrature weights
c     nquad        number of quadrature nodes
c     nboxes       number of boxes in FMM hierarchy
c     laddr        indexing array for FMM data structure
c     nlev         number of levels in FMM hierarchy
c     scale        array of scaling parameters
c     bsize        box dimension for FMM
c     nterms       array of nterms needed at each level
c     wlists       FMM data structure (real array)
c     lw           length of wlists
c
c
c     OUTPUT PARAMETERS:
c
c     pot          surface potential (if ifpot=1)
c     fld          surface field=-gradient(potential) (if iffld=1)
c     pottarg      target potential (if ifpot=1)
c     fldtarg      target field=-gradient(potential) (if iffld=1)
c     ier          error return code
c
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c       
        ifprint=0
c
c     
c       ... set the potential and field to zero
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
        do i=1,10
        timeinfo(i)=0
        enddo
c
c
        norder=1
        nqtri=1
c
        if( iprec .eq. -2 ) then
        norder=2
        nqtri=2
        endif
c
        if( iprec .eq. -1 ) then
        norder=2
        nqtri=2
        endif
c
        if( iprec .eq. 0 ) then
        norder=4
        nqtri=4
        endif
c
        if( iprec .ge. 1 ) then
        norder=6
        nqtri=6
        endif
c
        if( ifevalfar .eq. 0 ) goto 8000
c       
c       ... initialize Legendre function evaluation routines
c
        nlege=100
        lw7=40000
        call ylgndrfwini(nlege,wlege,lw7,lused7)
ccc        write(*,*)' lused7 from ylgndrfwini is',lused7
c
        do i=0,nlev
        do itype=1,4
        call l3dterms_eval(itype,epsfmm,
     1       nterms_eval(itype,i),ier)
        enddo
        enddo
c
        if (ifprint.ge.2) 
     $     call prinf('nterms_eval=*',nterms_eval,4*(nlev+1))
c
c       ... set all multipole and local expansions to zero
c
        do ibox = 1,nboxes
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        level=box(1)
        call l3dzero(rmlexp(iaddr(1,ibox)),nterms(level))
        call l3dzero(rmlexp(iaddr(2,ibox)),nterms(level))
        enddo
c
c
        if (ifprint.ge.1) 
     $     call prinf('=== STEP 1 (form mp) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 1, locate all charges, assign them to boxes, and
c       form multipole expansions
c
ccc        do 1200 ibox=1,nboxes
        do 1300 ilev=3,nlev+1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,radius)
C$OMP$PRIVATE(ier,i,j,ptemp,ftemp,cd) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
        do 1200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
c
        level=box(1)
c
c
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
c        ipts=box(14)
c        npts=box(15)
c        call prinf('ipts=*',ipts,1)
c        call prinf('npts=*',npts,1)
        npts=box(15)
        if (ifprint .ge. 2) then
           call prinf('npts=*',npts,1)
           call prinf('isource=*',isource(box(14)),box(15))
        endif
        endif
c
c       ... prune all sourceless boxes
c
        if( box(15) .eq. 0 ) goto 1200
c
        if (nkids .eq. 0) then
c
c       ... form multipole expansions
c
	    radius = (corners0(1,1) - center0(1))**2
	    radius = radius + (corners0(2,1) - center0(2))**2
	    radius = radius + (corners0(3,1) - center0(3))**2
	    radius = sqrt(radius)
c
 	    if (ifcharge .eq. 1) then
               call l3dformmptris2_add(ier,scale(level),
     1  	  triaflatsort(1,1,box(14)),chargesort(box(14)),npts,
     2            center0,norder,nterms(level),rmlexp(iaddr(1,ibox)),
     $            wlege,nlege)
            endif
c
            if (ifdipole .eq. 1 ) then
               call l3dformmptrid2_add(ier,scale(level),
     1           triaflatsort(1,1,box(14)),trianormsort(1,box(14)),
     2           dipstrsort(box(14)),npts,center0,norder,
     3           nterms(level),rmlexp(iaddr(1,ibox)),wlege,nlege)
            endif

         endif
c
 1200    continue
C$OMP END PARALLEL DO
 1300    continue
c
         t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
         timeinfo(1)=t2-t1
c       
        if (ifprint.ge.1) 
     $     call prinf('=== STEP 2 (form lo) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 2, adaptive part, form local expansions, 
c           or evaluate the potentials and fields directly
c 
         do 3251 ibox=1,nboxes
c
         call d3tgetb(ier,ibox,box,center0,corners0,wlists)
c
         itype=3
         call d3tgetl(ier,ibox,itype,list,nlist,wlists)
         if (nlist .gt. 0) then 
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('list3=*',list,nlist)
            endif
         endif
c
c       ... prune all sourceless boxes
c
         if( box(15) .eq. 0 ) nlist=0
c
c
c       ... note that lists 3 and 4 are dual
c
c       ... form local expansions for all boxes in list 3
c       ... if target is childless, evaluate directly (if cheaper)
c        
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(level,npts,nkids)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect3,radius)
C$OMP$PRIVATE(ier,i,j,ptemp,ftemp,cd,ilist) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
         do ilist=1,nlist
            jbox=list(ilist)
            call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c        
            level1=box1(1)
c
            ifdirect3 = 0
            if( box1(15) .lt. (nterms(level1)+1)**2/4 .and.
     $          box(15) .lt. (nterms(level1)+1)**2/4 ) ifdirect3 = 1
c
            ifdirect3 = 0
c
            if( ifdirect3 .eq. 0 ) then
c
               npts=box(15)
c
               if( ifcharge .eq. 1 ) then
                  call l3dformtatris2_add(ier,scale(level1),
     1  	    triaflatsort(1,1,box(14)),chargesort(box(14)),
     $              npts,center1,
     2              norder,nterms(level1),
     $              rmlexp(iaddr(2,jbox)),wlege,nlege)
	       endif
               if( ifdipole .eq. 1 ) then
                  call l3dformtatrid2_add(ier,scale(level1),
     1               triaflatsort(1,1,box(14)),trianormsort(1,box(14)),
     2               dipstrsort(box(14)),npts,center1,
     3               norder,nterms(level1),
     $               rmlexp(iaddr(2,jbox)),wlege,nlege)
               endif
c
            else

            call lfmm3dtria_direct(nqtri,box,box1,
     $         triaflatsort,trianormsort,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         ifpot,pot,iffld,fld,
     $         targetsort,ifpottarg,pottarg,iffldtarg,fldtarg)

            endif
         enddo
C$OMP END PARALLEL DO
c
 3251    continue
c
         t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
         timeinfo(2)=t2-t1
c
c
        if( ifprint .ge. 1 )
     $      call prinf('=== STEPS 3,4,5 ====*',i,0)
        ifprune_list2 = 1
        if (ifpot.eq.1) ifprune_list2 = 0
        if (iffld.eq.1) ifprune_list2 = 0
        call lfmm3d_list2
     $     (bsize,nlev,laddr,scale,nterms,rmlexp,iaddr,epsfmm,
     $     timeinfo,wlists,mptemp,lmptemp,
     $     ifprune_list2)

c
        if (ifprint .ge. 1) call prinf('=== STEP 6 (eval mp) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 6, adaptive part, evaluate multipole expansions, 
c           or evaluate the potentials and fields directly
c
         do 3252 ibox=1,nboxes
         call d3tgetb(ier,ibox,box,center0,corners0,wlists)
c
         itype=4
         call d3tgetl(ier,ibox,itype,list,nlist,wlists)
         if (nlist .gt. 0) then 
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('list4=*',list,nlist)
            endif
         endif
c
c       ... prune all sourceless boxes
c
         if( box(15) .eq. 0 ) nlist=0
c
c       ... note that lists 3 and 4 are dual
c
c       ... evaluate multipole expansions for all boxes in list 4 
c       ... if source is childless, evaluate directly (if cheaper)
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1)
C$OMP$PRIVATE(ier,i,j,ptemp,ftemp,cd,ilist) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
         do ilist=1,nlist
            jbox=list(ilist)
            call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c
            level=box(1)
c
            ifdirect4 = 0
c
            if (box1(15) .lt. (nterms(level)+1)**2/4 .and.
     $         box(15) .lt. (nterms(level)+1)**2/4 ) ifdirect4 = 1
c
            ifdirect4 = 0
c
            if (ifdirect4 .eq. 0) then

            call l3dmpevalall_trunc(scale(level),center0,
     $         rmlexp(iaddr(1,ibox)),
     $         nterms(level),nterms_eval(1,level),
     $         sourcesort(1,box1(14)),box1(15),
     $         ifpot,pot(box1(14)),
     $         iffld,fld(1,box1(14)),
     $         wlege,nlege,ier)

            call l3dmpevalall_trunc(scale(level),center0,
     $         rmlexp(iaddr(1,ibox)),
     $         nterms(level),nterms_eval(1,level),
     $         targetsort(1,box1(16)),box1(17),
     $         ifpottarg,pottarg(box1(16)),
     $         iffldtarg,fldtarg(1,box1(16)),
     $         wlege,nlege,ier)

            else
            
            call lfmm3dtria_direct(nqtri,box,box1,
     $         triaflatsort,trianormsort,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         ifpot,pot,iffld,fld,
     $         targetsort,ifpottarg,pottarg,iffldtarg,fldtarg)
c
            endif
        enddo
C$OMP END PARALLEL DO
 3252   continue
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(6)=t2-t1
c

        if (ifprint .ge. 1) call prinf('=== STEP 7 (eval lo) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 7, evaluate local expansions
c       and all fields directly
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,ier)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 6201 ibox=1,nboxes
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
c
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
            npts=box(15)
            if (ifprint .ge. 2) then
               call prinf('npts=*',npts,1)
               call prinf('isource=*',isource(box(14)),box(15))
            endif
        endif
c
        if (nkids .eq. 0) then
c
c       ... evaluate local expansions
c       
        level=box(1)
        npts=box(15)
c       
        if (level .ge. 2) then

        call l3dtaevalall_trunc(scale(level),center0,
     $     rmlexp(iaddr(2,ibox)),
     $     nterms(level),nterms_eval(1,level),
     $     sourcesort(1,box(14)),box(15),
     $     ifpot,pot(box(14)),
     $     iffld,fld(1,box(14)),
     $     wlege,nlege,ier)

        call l3dtaevalall_trunc(scale(level),center0,
     $     rmlexp(iaddr(2,ibox)),
     $     nterms(level),nterms_eval(1,level),
     $     targetsort(1,box(16)),box(17),
     $     ifpottarg,pottarg(box(16)),
     $     iffldtarg,fldtarg(1,box(16)),
     $     wlege,nlege,ier)
        
        endif
c
        endif
c
 6201   continue
C$OMP END PARALLEL DO
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(7)=t2-t1
c
c
 8000   continue
c
c
c
        if( ifevalloc .eq. 0 ) goto 9000
c 
        if (ifprint.ge.1) 
     $     call prinf('=== STEP 8 (direct) =====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 8, evaluate direct interactions 
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,nkids,list,nlist,npts)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect2,radius)
C$OMP$PRIVATE(lused,ier,i,j,ptemp,ftemp,cd,ilist,itype) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
        do 6202 ibox=1,nboxes
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
c
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
            npts=box(15)
            if (ifprint .ge. 2) then
               call prinf('npts=*',npts,1)
               call prinf('isource=*',isource(box(14)),box(15))
            endif
        endif
c
c
        if (nkids .eq. 0) then
c
c       ... evaluate self interactions
c
        call lfmm3dtria_direct_self(nqtri,box,
     $     triaflatsort,trianormsort,sourcesort,
     $     ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $     ifpot,pot,iffld,fld,
     $     targetsort,ifpottarg,pottarg,iffldtarg,fldtarg)
c
c
c       ... retrieve list #1
c
c       ... evaluate interactions with the nearest neighbours
c
        itype=1
        call d3tgetl(ier,ibox,itype,list,nlist,wlists)
        if (ifprint .ge. 2) call prinf('list1=*',list,nlist)
c
c       ... for all pairs in list #1, 
c       evaluate the potentials and fields directly
c
            do 6203 ilist=1,nlist
               jbox=list(ilist)
               call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c
c       ... prune all sourceless boxes
c
          if( box1(15) .eq. 0 ) goto 6203
c    
            call lfmm3dtria_direct(nqtri,box1,box,
     $         triaflatsort,trianormsort,sourcesort,
     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     $         ifpot,pot,iffld,fld,
     $         targetsort,ifpottarg,pottarg,iffldtarg,fldtarg)
c
 6203       continue
        endif
c
 6202   continue
C$OMP END PARALLEL DO
c
ccc        call prin2('inside fmm, pot=*',pot,2*nsource)
c
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(8)=t2-t1
c
 9000   continue
c
ccc        call prinf('=== DOWNWARD PASS COMPLETE ===*',i,0)
c
        if (ifprint.ge.1) then
          call prin2('timeinfo=*',timeinfo,8)
          call prinf('nboxes=*',nboxes,1)
          call prinf('nsource=*',nsource,1)
          call prinf('ntarget=*',ntarget,1)
        endif
c       
        return
        end
c
c
c
c
c
        subroutine lfmm3dtria_direct_self(nqtri,box,
     $     triaflat,trianorm,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,
     $     target,ifpottarg,pottarg,iffldtarg,fldtarg)
        implicit real *8 (a-h,o-z)
c
        integer box(20),box1(20)
c
        dimension triaflat(3,3,1),trianorm(3,1)
c
        dimension source(3,1),dipvec(3,1)
        complex *16 charge(1),dipstr(1)
        dimension target(3,1)
c
        complex *16 pot(1),fld(3,1),pottarg(1),fldtarg(3,1)
        complex *16 ptemp,ftemp(3)
c
        ione = 1
c
        if( ifpot .eq. 1 .or. iffld .eq. 1 ) then
c
ccC$OMP PARALLEL DO DEFAULT(SHARED)
ccC$OMP$PRIVATE(i,j,ptemp,ftemp)
ccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=box(14),box(14)+box(15)-1
        do i=box(14),box(14)+box(15)-1
c
        if (ifcharge.ne.0) then
        if( i .eq. j ) then
        call direct3dtrilaps2(ione,ione,
     1     source(1,j),charge(j),triaflat(1,1,j),
     $     ifpot,ptemp,iffld,ftemp)
        else
        call direct3dtritarglaps2(ione,source(1,j),
     $     charge(i),triaflat(1,1,i),
     $     ifpot,ptemp,iffld,ftemp)
        endif
        if (ifpot.eq.1) pot(j)=pot(j)+ptemp
        if (iffld.eq.1) then
        fld(1,j)=fld(1,j)+ftemp(1)
        fld(2,j)=fld(2,j)+ftemp(2)
        fld(3,j)=fld(3,j)+ftemp(3)
        endif
        endif
c
        if (ifdipole.ne.0) then
        if( i .eq. j ) then
        call direct3dtrilapd2(ione,ione,
     $     source(1,j),dipstr(j),triaflat(1,1,j),
     $     trianorm(1,j),ifpot,ptemp,iffld,ftemp)
        else
        call direct3dtritarglapd2(ione,source(1,j),
     $     dipstr(i),triaflat(1,1,i),
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
ccC$OMP END PARALLEL DO
c
        endif
c
        if( ifpottarg .eq. 1 .or. iffldtarg .eq. 1 ) then
c
ccC$OMP PARALLEL DO DEFAULT(SHARED)
ccC$OMP$PRIVATE(i,j,ptemp,ftemp)
ccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=box(16),box(16)+box(17)-1
        do i=box(14),box(14)+box(15)-1
c
        ione = 1
        if (ifcharge.ne.0) then
        call direct3dtritarglaps2(ione,target(1,j),
     $     charge(i),triaflat(1,1,i),
     $     ifpottarg,ptemp,iffldtarg,ftemp)
        if (ifpottarg.eq.1) pottarg(j)=pottarg(j)+ptemp
        if (iffldtarg.eq.1) then
        fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
        fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
        fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
        endif
        endif
c
        if (ifdipole.ne.0) then
        call direct3dtritarglapd2(ione,target(1,j),
     $     dipstr(i),triaflat(1,1,i),
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
ccC$OMP END PARALLEL DO
c
        endif
c
        return
        end
c
c
c
c
c
        subroutine lfmm3dtria_direct(nqtri,box,box1,
     $     triaflat,trianorm,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,
     $     target,ifpottarg,pottarg,iffldtarg,fldtarg)
        implicit real *8 (a-h,o-z)
c
        integer box(20),box1(20)
c
        dimension triaflat(3,3,1),trianorm(3,1)
c
        dimension source(3,1),dipvec(3,1)
        complex *16 charge(1),dipstr(1)
        dimension target(3,1)
c
        complex *16 pot(1),fld(3,1),pottarg(1),fldtarg(3,1)
        complex *16 ptemp,ftemp(3)
c
c
        if( ifpot .eq. 1 .or. iffld .eq. 1 ) then
c
ccC$OMP PARALLEL DO DEFAULT(SHARED)
ccC$OMP$PRIVATE(i,j,ptemp,ftemp)
ccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=box1(14),box1(14)+box1(15)-1
c
        if (ifcharge.ne.0) then
        call direct3dtritarglaps2(box(15),source(1,j),
     $     charge(box(14)),triaflat(1,1,box(14)),
     $     ifpot,ptemp,iffld,ftemp)
        if (ifpot.eq.1) pot(j)=pot(j)+ptemp
        if (iffld.eq.1) then
        fld(1,j)=fld(1,j)+ftemp(1)
        fld(2,j)=fld(2,j)+ftemp(2)
        fld(3,j)=fld(3,j)+ftemp(3)
        endif
        endif
c
        if (ifdipole.ne.0) then
        call direct3dtritarglapd2(box(15),source(1,j),
     $     dipstr(box(14)),triaflat(1,1,box(14)),
     $     trianorm(1,box(14)),ifpot,ptemp,iffld,ftemp)
        if (ifpot.eq.1) pot(j)=pot(j)+ptemp
        if (iffld.eq.1) then
        fld(1,j)=fld(1,j)+ftemp(1)
        fld(2,j)=fld(2,j)+ftemp(2)
        fld(3,j)=fld(3,j)+ftemp(3)
        endif
        endif
c
        enddo
ccC$OMP END PARALLEL DO
c
        endif
c
        if( ifpottarg .eq. 1 .or. iffldtarg .eq. 1 ) then
c
ccC$OMP PARALLEL DO DEFAULT(SHARED)
ccC$OMP$PRIVATE(i,j,ptemp,ftemp)
ccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do j=box1(16),box1(16)+box1(17)-1
c
        if (ifcharge.ne.0) then
        call direct3dtritarglaps2(box(15),target(1,j),
     $     charge(box(14)),triaflat(1,1,box(14)),
     $     ifpottarg,ptemp,iffldtarg,ftemp)
        if (ifpottarg.eq.1) pottarg(j)=pottarg(j)+ptemp
        if (iffldtarg.eq.1) then
        fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
        fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
        fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
        endif
        endif
c
        if (ifdipole.ne.0) then
        call direct3dtritarglapd2(box(15),target(1,j),
     $     dipstr(box(14)),triaflat(1,1,box(14)),
     $     trianorm(1,box(14)),ifpottarg,ptemp,iffldtarg,ftemp)
        if (ifpottarg.eq.1) pottarg(j)=pottarg(j)+ptemp
        if (iffldtarg.eq.1) then
        fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
        fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
        fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
        endif
        endif
c
        enddo
ccC$OMP END PARALLEL DO
c
        endif
c
        return
        end
c
c
c
c
c
        subroutine l3dtriadirect(nsource,
     $     triaflat,trianorm,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,ntarget,
     $     target,ifpottarg,pottarg,iffldtarg,fldtarg)
        implicit real *8 (a-h,o-z)
c
c       Laplace interactions in R^3: evaluate all pairwise triangle
c       interactions (including self-interaction) + interactions with targets
c       via direct O(N^2) algorithm.
c
c       This is the principal subroutine for evaluating 
c       harmonic layer potentials on (flat) triangulated surfaces.
c       It permits the evaluation of a single layer potential
c       with piecewise constant density defined by <<charge>>
c       and a dipole layer potential with piecewise constant density 
c       and dipole orientation defined by <<dipstr,dipvec>>.
c
c       We use (1/r) for the Green's function,
c       without the (1/4 pi) scaling.  Self-interactions are included.
c   
c       It is capable of evaluating the layer potentials either on 
c       or off the surface (or both).            
c
c       INPUT PARAMETERS:
c
c       nsource: integer:  number of triangles
c       triaflat: real *8 (3,3,nsource): triangle coordinate array
c       trianorm: real *8 (3,nsource): triangle normals
c       source: real *8 (3,nsource):  triangle centroids
c       ifcharge:  single layer potential (SLP) flag
c                  ifcharge = 1   =>  include SLP contribution
c                                     otherwise do not
c       charge: complex *16 (nsource): piecewise constant SLP strength
c       ifdipole:  dipole layer potential (DLP) flag
c                  ifdipole = 1   =>  include DLP contribution
c                                     otherwise do not
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
        dimension triaflat(3,3,1),trianorm(3,1)
c
        dimension source(3,1),dipvec(3,1)
        complex *16 charge(1),dipstr(1)
        dimension target(3,1)
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
C$OMP$PRIVATE(i,j,ptemp,ftemp) 
        do j=1,nsource
        do i=1,nsource
c
        if (ifcharge.ne.0) then
        if( i .eq. j ) then
        call direct3dtrilaps2(ione,ione,
     1     source(1,j),charge(j),triaflat(1,1,j),
     $     ifpot,ptemp,iffld,ftemp)
        else
        call direct3dtritarglaps2(ione,source(1,j),
     $     charge(i),triaflat(1,1,i),
     $     ifpot,ptemp,iffld,ftemp)
        endif
        if (ifpot.eq.1) pot(j)=pot(j)+ptemp
        if (iffld.eq.1) then
        fld(1,j)=fld(1,j)+ftemp(1)
        fld(2,j)=fld(2,j)+ftemp(2)
        fld(3,j)=fld(3,j)+ftemp(3)
        endif
        endif
c
        if (ifdipole.ne.0) then
        if( i .eq. j ) then
        call direct3dtrilapd2(ione,ione,
     $     source(1,j),dipstr(j),triaflat(1,1,j),
     $     trianorm(1,j),ifpot,ptemp,iffld,ftemp)
        else
        call direct3dtritarglapd2(ione,source(1,j),
     $     dipstr(i),triaflat(1,1,i),
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
C$OMP$PRIVATE(i,j,ptemp,ftemp) 
        do j=1,ntarget
        do i=1,nsource
c
        if (ifcharge.ne.0) then
        call direct3dtritarglaps2(ione,target(1,j),
     $     charge(i),triaflat(1,1,i),
     $     ifpottarg,ptemp,iffldtarg,ftemp)
        if (ifpottarg.eq.1) pottarg(j)=pottarg(j)+ptemp
        if (iffldtarg.eq.1) then
        fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
        fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
        fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
        endif
        endif
c
        if (ifdipole.ne.0) then
        call direct3dtritarglapd2(ione,target(1,j),
     $     dipstr(i),triaflat(1,1,i),
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
        subroutine lfmm3dtriampftarg(ier,iprec,nsource,
     $     triaflat,trianorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ntarget,target,ifpottarg,pottarg,iffldtarg,fldtarg)
        implicit real *8 (a-h,o-z)
c
c       Laplace FMM in R^3: evaluate interactions with targets (off surface)
c
c       Harmonic layer potentials due to flat triangulated surfaces.
c
c       It permits the evaluation of a single layer potential
c       with piecewise constant density defined by <<charge>>
c       and a dipole layer potential with piecewise constant density 
c       and dipole orientation defined by <<dipstr,dipvec>>.
c
c       It is capable of evaluating the layer potentials off the surface only,
c       the far away targets are evaluated via an outgoing multipole expansion.
c       
c       We use (1/r) for the Green's function,
c       without the (1/4 pi) scaling. 
c   
c       This is primarily a memory management code. 
c       The actual work is carried out in subroutine lfmm3dtriampftarg0.
c
c       INPUT PARAMETERS:
c
c       iprec:  FMM precision flag
c
c                 -2 => tolerance =.5d0
c                 -1 => tolerance =.5d-1
c                  0 => tolerance =.5d-2
c                  1 => tolerance =.5d-3
c                  2 => tolerance =.5d-6
c                  3 => tolerance =.5d-9
c                  4 => tolerance =.5d-12
c                  5 => tolerance =.5d-15
c
c       nsource: integer:  number of triangles
c       triaflat: real *8 (3,3,nsource): triangle coordinate array
c       trianorm: real *8 (3,nsource): triangle normals
c       source: real *8 (3,nsource):  triangle centroids
c       ifcharge:  single layer potential (SLP) flag
c                  ifcharge = 1   =>  include SLP contribution
c                                     otherwise do not
c       charge: complex *16 (nsource): piecewise constant SLP strength
c       ifdipole:  dipole layer potential (DLP) flag
c                  ifdipole = 1   =>  include DLP contribution
c                                     otherwise do not
c       dipstr: complex *16 (nsource): piecewise constant DLP strengths
c       dipvec: real *8 (3,nsource): piecewise constant dipole orientation 
c                                    vectors. 
c
c           NOTE: In the present version, dipvec MUST BE SET EQUAL
c                 to the triangle normal. It is here as an additional
c                 parameter for future use, where an arbitrarily 
c                 oriented dipole vector is permitted.
c
c       ntarget: integer:  number of targets
c       target: real *8 (3,ntarget):  target locations
c       ifpottarg:  target potential flag 
c                   (1=compute potential, otherwise no)
c       iffldtarg:  target field flag 
c                   (1=compute field, otherwise no)
c
c       OUTPUT PARAMETERS:
c
c       ier   =  error return code
c                ier=0     =>  normal execution
c                ier=4     =>  cannot allocate tree workspace
c                ier=8     =>  cannot allocate bulk FMM  workspace
c                ier=16    =>  cannot allocate mpole expansion
c                              workspace in FMM
c
c       pottarg: complex *16 (ntarget): potential at target locations 
c       fldtarg: complex *16 (3,ntarget): field (-gradient) at target locations 
c
        dimension triaflat(3,3,1)
        dimension trianorm(3,1)
        dimension source(3,1)
        complex *16 charge(1)
        complex *16 dipstr(1)
        dimension dipvec(3,1)
        dimension target(3,1)
        complex *16 pottarg(1)
        complex *16 fldtarg(3,1)        
c
        real *8, allocatable :: w(:)
        dimension center(3),corners(3,8)
c
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        lused=0
c
        iiseptarg=1
        liseptarg=ntarget
        lused=lused+liseptarg
c
        itarget1=lused+1
        ltarget1=3*ntarget
        lused=lused+ltarget1
c
        ipottarg1=lused+1
        lpottarg1=2*ntarget
        lused=lused+lpottarg1
c
        ifldtarg1=lused+1
        lfldtarg1=2*3*ntarget
        lused=lused+lfldtarg1
        allocate (w(lused))
c        
        call lfmm3dtriampftarg0(ier,iprec,nsource,
     $     triaflat,trianorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ntarget,target,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     w(iiseptarg),w(itarget1),w(ipottarg1),w(ifldtarg1))
c
        return
        end
c
c
c
c
c
        subroutine lfmm3dtriampftarg0(ier,iprec,nsource,
     $     triaflat,trianorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ntarget,target,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     iseptarg,target1,pottarg1,fldtarg1)
        implicit real *8 (a-h,o-z)
        dimension triaflat(3,3,1)
        dimension trianorm(3,1)
        dimension source(3,1)
        complex *16 charge(1)
        complex *16 dipstr(1)
        dimension dipvec(3,1)
        dimension target(3,1)
        complex *16 pottarg(1)
        complex *16 fldtarg(3,1)        
c
        integer iseptarg(1)
        dimension target1(3,1)
        complex *16 pottarg1(1)
        complex *16 fldtarg1(3,1)        
c
        dimension center(3),corners(3,8)
c
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        lused7=0
c       
c       ... first, check if all targets are well-separated from the
c       computational box
c
        call d3tgetbbox(nsource,source,center,size,corners)
        do i=1,ntarget
        ifsep = 0
        if( abs(target(1,i)-center(1)) .gt. 1.5d0*size ) ifsep = 1
        if( abs(target(2,i)-center(2)) .gt. 1.5d0*size ) ifsep = 1
        if( abs(target(3,i)-center(3)) .gt. 1.5d0*size ) ifsep = 1
ccc        if( ifsep .eq. 0 ) call prinf('near field target, i=*',i,1)
        iseptarg(i)=ifsep
        enddo
c
ccc        call prinf('ntarget=*',ntarget,1)
c        
c
        ntarget1=0
        do i=1,ntarget
        if( iseptarg(i) .eq. 1 ) then
        ntarget1=ntarget1+1
        target1(1,ntarget1)=target(1,i)
        target1(2,ntarget1)=target(2,i)
        target1(3,ntarget1)=target(3,i)
        endif
        enddo
c
ccc        call prinf('ntarget1=*',ntarget1,1)
c
        if( ntarget1 .gt. 0 ) then
c
        call lfmm3dtriamptarg(ier,iprec,nsource,
     $     triaflat,trianorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ntarget1,target1,ifpottarg,pottarg1,iffldtarg,fldtarg1)
c
        itarget1=0
        do i=1,ntarget
        if( iseptarg(i) .eq. 1 ) then
        itarget1=itarget1+1
        if(ifpottarg .eq. 1) pottarg(i)=pottarg1(itarget1)
        if(iffldtarg .eq. 1) then
        fldtarg(1,i)=fldtarg1(1,itarget1)
        fldtarg(2,i)=fldtarg1(2,itarget1)
        fldtarg(3,i)=fldtarg1(3,itarget1)
        endif
        endif
        enddo
c
        endif
c

        ntarget2=0
        do i=1,ntarget
        if( iseptarg(i) .eq. 0 ) then
        ntarget2=ntarget2+1
        target1(1,ntarget2)=target(1,i)
        target1(2,ntarget2)=target(2,i)
        target1(3,ntarget2)=target(3,i)
        endif
        enddo
c
ccc        call prinf('ntarget2=*',ntarget2,1)
c
        if( ntarget2 .gt. 0 ) then
c
        ifpot=0
        iffld=0

        call lfmm3dtriatarg(ier,iprec,nsource,
     $     triaflat,trianorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ifpot,pot,iffld,fld,
     $     ntarget2,target1,ifpottarg,pottarg1,iffldtarg,fldtarg1)
c
        itarget2=0
        do i=1,ntarget
        if( iseptarg(i) .eq. 0 ) then
        itarget2=itarget2+1
        if(ifpottarg .eq. 1) pottarg(i)=pottarg1(itarget2)
        if(iffldtarg .eq. 1) then
        fldtarg(1,i)=fldtarg1(1,itarget2)
        fldtarg(2,i)=fldtarg1(2,itarget2)
        fldtarg(3,i)=fldtarg1(3,itarget2)
        endif
        endif
        enddo
c
        endif
c
        return
        end
c
c
c
c
c
        subroutine lfmm3dtriamptarg(ier,iprec,nsource,
     $     triaflat,trianorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ntarget,target,ifpottarg,pottarg,iffldtarg,fldtarg)
        implicit real *8 (a-h,o-z)
        dimension triaflat(3,3,1)
        dimension trianorm(3,1)
        dimension source(3,1)
        complex *16 charge(1),zk
        complex *16 dipstr(1)
        dimension dipvec(3,1)
        dimension target(3,1)
        complex *16 pottarg(1)
        complex *16 fldtarg(3,1)        
c
        real *8, allocatable :: w(:)
        dimension center(3),corners(3,8)
c
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        call d3tgetbbox(nsource,source,center,size,corners)
        scale=1
c
cc        call prin2('center=*',center,3)
cc        call prin2('size=*',size,1)
cc        call prin2('corners=*',corners,3*8)
c
        if( iprec .eq. -2 ) epsfmm=.5d-0 
        if( iprec .eq. -1 ) epsfmm=.5d-1
        if( iprec .eq. 0 ) epsfmm=.5d-2
        if( iprec .eq. 1 ) epsfmm=.5d-3
        if( iprec .eq. 2 ) epsfmm=.5d-6
        if( iprec .eq. 3 ) epsfmm=.5d-9
        if( iprec .eq. 4 ) epsfmm=.5d-12
        if( iprec .eq. 5 ) epsfmm=.5d-15
        if( iprec .eq. 6 ) epsfmm=0
c
        call l3dterms(epsfmm,nterms,ier)
ccc        call prinf('after l3dterms, nterms=*',nterms,1)
c
        lused=0
c
        impole=1
        lmpole=2*(nterms+1)*(2*nterms+1)
        lused=lused+lmpole
c
        imptemp=lused+1
        lmptemp=2*(nterms+1)*(2*nterms+1)
        lused=lused+lmptemp
        allocate (w(lused))
c
        call lfmm3dtriamptarg0(ier,iprec,nsource,
     $     triaflat,trianorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ntarget,target,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     center,size,corners,
     $     w(impole),w(imptemp),nterms,scale)

c
        return
        end
c
c
c
c
c
        subroutine lfmm3dtriamptarg0(ier,iprec,nsource,
     $     triaflat,trianorm,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     ntarget,target,ifpottarg,pottarg,iffldtarg,fldtarg,
     $     center0,size0,corners0,
     $     mpole,mptemp,nterms,scale)
c
        implicit real *8 (a-h,o-z)
        dimension triaflat(3,3,1)
        dimension trianorm(3,1)
        dimension source(3,1)
        complex *16 charge(1)
        complex *16 dipstr(1)
        dimension dipvec(3,1)
        dimension target(3,1)
        complex *16 pottarg(1)
        complex *16 fldtarg(3,1)        
c
        complex *16 mpole(1), mptemp(1)
        dimension center0(3), corners0(3,8)
c
        complex *16 ptemp,ftemp(3)
c
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
        call l3dzero(mpole,nterms)
c
c
c        call prinf('ntarget=*',ntarget,1)
c        call prinf('ifpottarg=*',ifpottarg,1)
c        call prinf('iffldtarg=*',iffldtarg,1)
c        if(ifpottarg.eq.1) call prin2('pottarg=*',pottarg,2*12)
c        if(iffldtarg.eq.1) call prin2('fldtarg=*',fldtarg,6*12)
c
c
        norder=1
        nqtri=1
c
        if( iprec .eq. -2 ) then
        norder=1
        nqtri=1
        endif
c
        if( iprec .eq. -1 ) then
        norder=2
        nqtri=2
        endif
c
        if( iprec .eq. 0 ) then
        norder=3
        nqtri=3
        endif
c
        if( iprec .ge. 1 ) then
        norder=4
        nqtri=6
        endif
c
c
        npts=nsource
c
c       
c       ... form outgoing multipole expansion
c
        radius = (corners0(1,1) - center0(1))**2
        radius = radius + (corners0(2,1) - center0(2))**2
        radius = radius + (corners0(3,1) - center0(3))**2
        radius = sqrt(radius)
c       
        if (ifcharge .eq. 1) then
        call l3dformmptris_add(ier,scale,
     1     triaflat,charge,npts,
     2     center0,norder,nterms,mpole)        
        endif
c
        if (ifdipole .eq. 1 ) then
        call l3dformmptrid_add(ier,scale,
     1     triaflat,trianorm,
     2     dipstr,npts,center0,norder,
     3     nterms,mpole)
        endif
c
c
        do j=1,ntarget

        call l3dmpeval(scale,center0,
     $     mpole,nterms,target(1,j),
     $     ptemp,iffldtarg,ftemp,ier)

        if( ifpottarg .eq. 1 ) pottarg(j)=pottarg(j)+ptemp
        if( iffldtarg .eq. 1 ) then
        fldtarg(1,j)=fldtarg(1,j)+ftemp(1)
        fldtarg(2,j)=fldtarg(2,j)+ftemp(2)
        fldtarg(3,j)=fldtarg(3,j)+ftemp(3)
        endif

        enddo
c
c
        return
        end
c
c
c
c
c
