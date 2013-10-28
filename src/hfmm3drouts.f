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
c    $Date: 2010-07-05 09:56:35 -0400 (Mon, 05 Jul 2010) $
c    $Revision: 1048 $
c
c       
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the routines for Helmholtz FMM in R^3
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine hfmm3dparttree(ier,iprec,zk,
     $     nsource,source,ntarget,target,
     $     nbox,epsfmm,iisource,iitarget,iwlists,lwlists,
     $     nboxes,laddr,nlev,center,size,
     $     w,lw,lused7)
        implicit real *8 (a-h,o-z)
c       
c       Helmholtz FMM in R^3: build the oct-tree
c
c     INPUT PARAMETERS:
c
c       zk: complex *16: Helmholtz parameter
c       nsource: integer:  number of sources
c       source: real *8 (3,n):  source locations
c       w: real *8 (lw): workspace
c       lw:  length of workspace
c
c     OUTPUT PARAMETERS:
c
c       ier   =  error return code
c       lused7 = the amount of workspace used
c
c
        dimension source(3,1),target(3,1)
c       
        dimension center(3)
c       
        dimension laddr(2,200)
        integer box(20)
        dimension center0(3),corners0(3,8)
c       
        integer box1(20)
        dimension center1(3),corners1(3,8)
c
        complex *16 zk
        dimension w(1)
c       
        ier=0
c       
        done=1
        pi=4*atan(done)
c       
        lused7=0
        ifprint=0
c       
        iisource=1
        lused7=lused7+nsource
        if (ifprint.eq.1) call prinf('lused7=*',lused7,1)
        if (lused7 .ge. lw) ier=128
        if( ier .ne. 0 ) return
c
        iitarget=iisource+nsource
        lused7=lused7+ntarget
        if (ifprint.eq.1) call prinf('lused7=*',lused7,1)
        if (lused7 .ge. lw) ier=128
        if( ier .ne. 0 ) return
c
        iwlists=iisource+lused7+10
c
c       ... construct the adaptive FMM oct-tree structure
c
        call d3tstrcr(ier,source,nsource,nbox,
     $     nboxes,w(iisource),laddr,nlev,center,size,
     $     target,ntarget,w(iitarget),w(iwlists),lw-lused7,lused)
c
        if( ier .ne. 0 ) return
c
        lwlists=lused
        lused7=lused7+lwlists
        if (lused7 .ge. lw) ier=128
        if( ier .ne. 0 ) return
c       
        if (ifprint.eq.1) 
     $       call prin2('after d3tstrcr, center=*',center,3)
        if (ifprint.eq.1) 
     $       call prin2('after d3tstrcr, size=*',size,1)
        if (ifprint.eq.1) 
     $       call prinf('after d3tstrcr, nlev=*',nlev,1)
        if (ifprint.eq.1) 
     $       call prinf('after d3tstrcr, nbox=*',nbox,1)
        if (ifprint.eq.1) 
     $       call prinf('after d3tstrcr, laddr=*',laddr,2*(nlev+1))
c
ccc        call prinf('after d3tstrcr, isource=*',w(iisource),nsource)
ccc        call prinf('after d3tstrcr, itarget=*',w(iitarget),ntarget)
c
ccc        call d3tprint(w(iwlists),lwlists)
c
c       ... optional, plot the oct-tree in gnuplot compatible format
c
        ifplot = 0
        if (ifplot .eq. 1 .and. nsource .lt. 10000 ) then
c
c       ... plot the boxes
c
        iw=51
        call plot_box3d(iw,center,size)
c       
        itag=1
        iw=52
        call plot_label3d(iw,center,size,itag,itag)
c
        iw=60
        call plot_points3d(iw,source,nsource)
c       
        iw=63
        call plot_points3d(iw,target,ntarget)
c       
        do ibox=1,nboxes
           call d3tgetb(ier,ibox,box,center0,corners0,w(iwlists))
           level=box(1)
           size0=size/2**level
c
           iw=61
           call plot_box3d(iw,center0,size0)
c       
           itag=ibox
           iw=62
           call plot_label3d(iw,center0,size0,itag,itag)
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
        subroutine hfmm3d_list2
     $     (zk,bsize,nlev,laddr,scale,nterms,rmlexp,iaddr,epsfmm,
     $     timeinfo,wlists,mptemp,lmptemp,xnodes,wts,nquad,
     $     ifprune_list2)
        implicit real *8 (a-h,o-z)
c
        integer iaddr(2,1),laddr(2,1),nterms(0:1)
        dimension rmlexp(1),scale(0:1),itable(-3:3,-3:3,-3:3)
c
        integer list(10 000)
c
        integer box(20)
        dimension bsize(0:200)
        dimension xnodes(nquad)
        dimension wts(nquad)
        dimension center0(3),corners0(3,8)
c       
        integer box1(20)
        dimension center1(3),corners1(3,8)
c
        dimension wlists(1)
        complex *16 mptemp(lmptemp)
        complex *16 zk
c
        complex *16 pot(1),fld(3,1),pottarg(1),fldtarg(3,1)
        complex *16 ptemp,ftemp(3)
c
        real *8, allocatable :: rotmatf(:,:,:,:)
        real *8, allocatable :: rotmatb(:,:,:,:)
        real *8, allocatable :: thetas(:,:,:)
        real *8 rvec(3)
c
        dimension timeinfo(10)

        real *8, allocatable :: xnodes2(:), wts2(:)

c       
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c       
        ifprint=0
c
        max_nodes = 10000
        allocate( xnodes2(max_nodes) )
        allocate( wts2(max_nodes) )
c
c
         if (ifprint .ge. 1) 
     $     call prinf('=== STEP 3 (merge mp) ===*',i,0)
         t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 3, merge all multipole expansions
c       
ccc         do 2200 ibox=nboxes,1,-1
         do 2300 ilev=nlev,3,-1
        nquad2=nterms(ilev-1)*2.5
        nquad2=max(6,nquad2)
ccc        write(*,*) 'in mpmp', ilev, nterms(ilev-1), nquad2
        ifinit2=1
        call legewhts(nquad2,xnodes2,wts2,ifinit2)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,level,npts,nkids,radius)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1)
C$OMP$PRIVATE(lused,ier,i,j,ptemp,ftemp,cd) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
         do 2200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
c
         call d3tgetb(ier,ibox,box,center0,corners0,wlists)
         call d3tnkids(box,nkids)
c
c       ... prune all sourceless boxes
c
         if( box(15) .eq. 0 ) goto 2200
c
         if (nkids .ne. 0) then
c
         level0=box(1)
         if( level0 .ge. 2 ) then
ccc         if (level0 .ge. 0) then
            radius = (corners0(1,1) - center0(1))**2
            radius = radius + (corners0(2,1) - center0(2))**2
            radius = radius + (corners0(3,1) - center0(3))**2
            radius = sqrt(radius)
c       
            if( ifprint .ge. 2 ) then
               call prin2('radius=*',radius,1)
               call prinf('ibox=*',ibox,1)
               call prinf('box=*',box,20)
               call prinf('nkids=*',nkids,1)
            endif
c
c       ... merge multipole expansions of the kids 
c
            call h3dzero(rmlexp(iaddr(1,ibox)),nterms(level0))
            if (ifprint .ge. 2) then
               call prin2('center0=*',center0,3)
            endif
c
            do 2100 i = 1,8
               jbox = box(5+i)
               if (jbox.eq.0) goto 2100
               call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
               if (ifprint .ge. 2) then
               call prinf('jbox=*',jbox,1)
               call prin2('center1=*',center1,3)
               endif
               level1=box1(1)
               call h3dmpmpquadu_add(zk,scale(level1),center1,
     1            rmlexp(iaddr(1,jbox)),nterms(level1),scale(level0),
     1            center0,rmlexp(iaddr(1,ibox)),
     $            nterms(level0),nterms(level0),
ccc     1            radius,xnodes,wts,nquad,ier)
     1            radius,xnodes2,wts2,nquad2,ier)
 2100       continue
            if (ifprint .ge. 2) then
            call prinf('=============*',x,0)
            endif
c       ... mark the local expansion of all kids and the parent
c
            endif
         endif
 2200    continue
C$OMP END PARALLEL DO
 2300    continue
c
c
ccc        call prinf('=== UPWARD PASS COMPLETE ===*',i,0)
c
c------------------------------------------------------------
c      DEBUGGING SEGMENT - once all multipole expansions are merged
c      to top level, one can compare it to the direct formation of the
c      expansion at the top level from the source locations.
c
ccc        call prinm(rmlexp(iaddr(1,1)),nterms(0))
c
ccc        call h3dformmp(ier,scale(0),source,charge,n,
ccc     1  	center,nterms(0),mptemp)
c
ccc        call prinm(mptemp,nterms(0))
c
ccc        call h3dmperr(rmlexp(iaddr(1,1)),mptemp,nterms(0),d)
ccc        call prin2('error in upward pass=*',d,1)
c
ccc        pause
ccc        stop
c      END DEBUGGING SEGMENT
c------------------------------------------------------------
c
         t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
         timeinfo(3)=t2-t1
c
        if (ifprint .ge. 1) 
     $     call prinf('=== STEP 4 (mp to lo) ===*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... precompute rotation matrices, useful up to order 10 or so
c       (approximately 30kB of storage for ldm=10)
c       (approximately 40MB of storage for ldm=30)
c
        ldm = 1
        do i=2,nlev
        if( nterms(i) .gt. ldm) ldm = nterms(i)
        enddo
        if( ldm .gt. 10 ) ldm = 10

        if( ifprint .ge. 2 ) call prinf('ldm=*',ldm,1)

        allocate(rotmatf((ldm+1)*(ldm+1)*(2*ldm+1),-3:3,-3:3,-3:3))
        allocate(rotmatb((ldm+1)*(ldm+1)*(2*ldm+1),-3:3,-3:3,-3:3))
        allocate(thetas(-3:3,-3:3,-3:3))
        nstor = (ldm+1)*(ldm+1)*(2*ldm+1)*7*7*7  * 2

        if( ifprint .ge. 2 ) then
        call prinf('nstor=*',nstor,1)
        call prinf('nstor=*',nstor/1000,1)
        call prinf('nstor=*',nstor/1000000,1)
        endif
c
        thetas(0,0,0)=0
        do i=-3,3
        do j=-3,3
        do k=-3,3
c        
        if( abs(i).gt.0 .or. abs(j).gt.0 .or. abs(k).gt.0 ) then
        rvec(1) = i
        rvec(2) = j
        rvec(3) = k
c
        call cart2polar(rvec,d,theta,phi)
        thetas(i,j,k)=theta
c
        call rotviarecur3p_init(ier,rotmatf(1,i,j,k),ldm,+theta)
        call rotviarecur3p_init(ier,rotmatb(1,i,j,k),ldm,-theta)
c
        endif
c
        enddo
        enddo
        enddo
        if( ifprint .ge. 2 ) call prin2('thetas=*',thetas,7*7*7)
c
c       ... step 4, convert multipole expansions into the local ones
c
cc        call prinf('laddr=*',laddr,2*(nlev+1))
cc        call prin2('bsize=*',bsize,(nlev+1))
cc        do 4200 ibox=1,nboxes
ccc        ntops=0
        do 4300 ilev=3,nlev+1
        call h3dterms_list2(bsize(ilev-1),zk,epsfmm, itable, ier)
ccc        call prinf('itable=*',itable,7*7*7)
ccc        t3=second()
cccC$        t3=omp_get_wtime()
c
c       ... truncate projection quadratures, spectrally accurate
c       for iprec = 0 and above, needs to be at least 6 for iprec=-2 or -1
        nquad2=nterms(ilev-1)*1.2
        nquad2=max(6,nquad2)
ccc        write(*,*) 'in mpta', ilev, nterms(ilev-1), nquad2
        ifinit2=1
        call legewhts(nquad2,xnodes2,wts2,ifinit2)
ccc        call prin2('xnodes=*',xnodes,nquad)
ccc        call prin2('xnodes2=*',xnodes2,nquad2)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,list,nlist)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect2,radius)
C$OMP$PRIVATE(lused,ier,i,j,ptemp,ftemp,cd,ilist,itype)
C$OMP$PRIVATE(if_use_trunc,nterms_trunc,ii,jj,kk) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 4200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
        endif
        level0=box(1)
        if (level0 .ge. 2) then
c
c       ... retrieve list #2
c
           itype=2
           call d3tgetl(ier,ibox,itype,list,nlist,wlists)
           if (ifprint .ge. 2) then
              call prinf('list2=*',list,nlist)
           endif
c
c       ... prune all sourceless boxes
c
ccc           if( box(15) .eq. 0 ) nlist=0
c
c       ... for all pairs in list #2, apply the translation operator 
c
           do 4150 ilist=1,nlist
              jbox=list(ilist)
              call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
              if (box1(15).eq.0) goto 4150
              if ((box(17).eq.0).and.(ifprune_list2.eq.1))
     $           goto 4150
              radius = (corners1(1,1) - center1(1))**2
              radius = radius + (corners1(2,1) - center1(2))**2
              radius = radius + (corners1(3,1) - center1(3))**2
              radius = sqrt(radius)
c
c       ... convert multipole expansions for all boxes in list 2 in local exp
c       ... if source is childless, evaluate directly (if cheaper)
c
              level1=box1(1)
c
              ifdirect2 = 0
c
              ii=box1(2)-box(2)
              jj=box1(3)-box(3)
              kk=box1(4)-box(4)
              nterms_trunc=itable(ii,jj,kk)
              nterms_trunc=min(nterms(level0),nterms_trunc)
              nterms_trunc=min(nterms(level1),nterms_trunc)
c
ccc              if_use_trunc=1
ccc              ntops=ntops+1
c
              if (ifdirect2 .eq. 0) then

              if_use_rotmatfb = 1
              if( nterms(level0) .gt. ldm ) if_use_rotmatfb = 0
              if( nterms(level1) .gt. ldm ) if_use_rotmatfb = 0
              if( nterms_trunc   .gt. ldm ) if_use_rotmatfb = 0
ccc              call prin2('theta=*',thetas(-ii,-jj,-kk),1)
              if( if_use_rotmatfb .eq. 1 ) then
              call h3dmplocquadu2_add_trunc(zk,scale(level1),center1,
     1           rmlexp(iaddr(1,jbox)),nterms(level1),nterms_trunc,
     $           scale(level0),center0,
     $           rmlexp(iaddr(2,ibox)),nterms(level0),nterms_trunc,
ccc     2            radius,xnodes,wts,nquad,ier)
     2           radius,xnodes2,wts2,nquad2,ier,
     $           rotmatf(1,-ii,-jj,-kk),rotmatb(1,-ii,-jj,-kk),ldm)
              else
              call h3dmplocquadu_add_trunc(zk,scale(level1),center1,
     1           rmlexp(iaddr(1,jbox)),nterms(level1),nterms_trunc,
     $           scale(level0),center0,
     $           rmlexp(iaddr(2,ibox)),nterms(level0),nterms_trunc,
ccc     2            radius,xnodes,wts,nquad,ier)
     2           radius,xnodes2,wts2,nquad2,ier)
              endif
              endif

 4150       continue
        endif
 4200   continue
C$OMP END PARALLEL DO
ccc        t4=second()
cccC$        t4=omp_get_wtime()
ccc        write(*,*) 'level ', ilev, ' time in list2:', t4-t3
ccc        write(*,*) 'time in list2:', second()-t1
ccc        write(*,*) 'ntops:', ntops
ccc        write(*,*) 'speed:', ntops/(second()-t1)
 4300   continue
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(4)=t2-t1
c       
        if (ifprint .ge. 1) 
     $     call prinf('=== STEP 5 (split lo) ===*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 5, split all local expansions
c
ccc        do 5200 ibox=1,nboxes
        do 5300 ilev=3,nlev
        nquad2=nterms(ilev-1)*2
        nquad2=max(6,nquad2)
ccc        write(*,*) 'in tata', ilev, nterms(ilev-1), nquad2
        ifinit2=1
        call legewhts(nquad2,xnodes2,wts2,ifinit2)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level0,level,npts,nkids,radius)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1)
C$OMP$PRIVATE(lused,ier,i,j,ptemp,ftemp,cd) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 5200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
c       
        if (nkids .ne. 0) then
            level0=box(1)
            if (level0 .ge. 2) then
               if (ifprint .ge. 2) then
                  call prinf('ibox=*',ibox,1)
                  call prinf('box=*',box,20)
                  call prinf('nkids=*',nkids,1)
                  call prin2('center0=*',center0,3)
               endif
c
c       ... split local expansion of the parent box
c
               do 5100 i = 1,8
	          jbox = box(5+i)
	          if (jbox.eq.0) goto 5100
                  call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
                  radius = (corners1(1,1) - center1(1))**2
                  radius = radius + (corners1(2,1) - center1(2))**2
                  radius = radius + (corners1(3,1) - center1(3))**2
                  radius = sqrt(radius)
                  if (ifprint .ge. 2) then
                     call prinf('jbox=*',jbox,1)
                     call prin2('radius=*',radius,1)
                     call prin2('center1=*',center1,3)
                  endif
                  level1=box1(1)
	          call h3dloclocquadu_add(zk,scale(level0),center0,
     1    	      rmlexp(iaddr(2,ibox)),nterms(level0),
     1                scale(level1),center1,rmlexp(iaddr(2,jbox)),
     $                nterms(level1),nterms(level1),
ccc     1    	      radius,xnodes,wts,nquad,ier)
     1    	      radius,xnodes2,wts2,nquad2,ier)
 5100          continue
               if (ifprint .ge. 2) call prinf('=============*',x,0)
            endif
        endif
c
        if (nkids .ne. 0) then
            level=box(1)
            if (level .ge. 2) then
               if( ifprint .ge. 2 ) then
                  call prinf('ibox=*',ibox,1)
                  call prinf('box=*',box,20)
                  call prinf('nkids=*',nkids,1)
               endif
            endif
        endif
 5200   continue
C$OMP END PARALLEL DO
 5300   continue
c       
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(5)=t2-t1
c
        return
        end
c
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the auxiliary routines
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine h3dpsort(n,isource,psort,pot)
        implicit real *8 (a-h,o-z)
        dimension isource(1)
        complex *16 pot(1),psort(1)
c
ccc        call prinf('isource=*',isource,n)
c        
        do i=1,n
        pot(isource(i))=psort(i)
        enddo
c
        return
        end
c
c
c
c
c
        subroutine h3dfsort(n,isource,fldsort,fld)
        implicit real *8 (a-h,o-z)
        dimension isource(1)
        complex *16 fld(3,1),fldsort(3,1)
c        
ccc        call prinf('isource=*',isource,n)
c
        do i=1,n
        fld(1,isource(i))=fldsort(1,i)
        fld(2,isource(i))=fldsort(2,i)
        fld(3,isource(i))=fldsort(3,i)
        enddo
c
        return
        end
c
c
c
c
c
        subroutine h3dpsortsub(n,isource,psort,pot)
        implicit real *8 (a-h,o-z)
        dimension isource(1)
        complex *16 pot(1),psort(1)
c
ccc        call prinf('isource=*',isource,n)
c        
        do i=1,n
        pot(isource(i))=pot(isource(i))-psort(i)
        enddo
c
        return
        end
c
c
c
c
c
        subroutine h3dfsortsub(n,isource,fldsort,fld)
        implicit real *8 (a-h,o-z)
        dimension isource(1)
        complex *16 fld(3,1),fldsort(3,1)
c        
ccc        call prinf('isource=*',isource,n)
c
        do i=1,n
        fld(1,isource(i))=fld(1,isource(i))-fldsort(1,i)
        fld(2,isource(i))=fld(2,isource(i))-fldsort(2,i)
        fld(3,isource(i))=fld(3,isource(i))-fldsort(3,i)
        enddo
c
        return
        end
c
c
c
c
c
        subroutine h3dreorder(nsource,source,
     $     ifcharge,charge,isource,ifdipole,
     1     dipstr,dipvec,sourcesort,chargesort,dipvecsort,dipstrsort) 
        implicit real *8 (a-h,o-z)
        dimension source(3,1),sourcesort(3,1),isource(1)
        dimension dipvec(3,1),dipvecsort(3,1)
        complex *16 charge(1),chargesort(1),dipstr(1),dipstrsort(1)
c       
ccc        call prinf('nsource=*',nsource,1)
        do i = 1,nsource
        sourcesort(1,i) = source(1,isource(i))
        sourcesort(2,i) = source(2,isource(i))
        sourcesort(3,i) = source(3,isource(i))
        if( ifcharge .ge. 1 ) then
        chargesort(i) = charge(isource(i))
        endif
        if (ifdipole .ge. 1) then
        dipstrsort(i) = dipstr(isource(i))
        dipvecsort(1,i) = dipvec(1,isource(i))
        dipvecsort(2,i) = dipvec(2,isource(i))
        dipvecsort(3,i) = dipvec(3,isource(i))
        endif
        enddo
        return
        end
c
c
c
c
c
        subroutine h3dreordertarg(ntarget,target,itarget,targetsort)
        implicit real *8 (a-h,o-z)
        dimension target(3,1),targetsort(3,1),itarget(1)
c       
ccc        call prinf('ntarget=*',ntarget,1)
        do i = 1,ntarget
        targetsort(1,i) = target(1,itarget(i))
        targetsort(2,i) = target(2,itarget(i))
        targetsort(3,i) = target(3,itarget(i))
        enddo
        return
        end
c
c
c
c
c
        subroutine h3dreordertria(nsource,isource,
     $     triaflat,triaflatsort,trianorm,trianormsort)
c
        implicit real *8 (a-h,o-z)
        dimension isource(1)
        dimension triaflat(3,3,1),triaflatsort(3,3,1)
        dimension trianorm(3,1),trianormsort(3,1)
c
        do i = 1,nsource
        triaflatsort(1,1,i) = triaflat(1,1,isource(i))
        triaflatsort(2,1,i) = triaflat(2,1,isource(i))
        triaflatsort(3,1,i) = triaflat(3,1,isource(i))
        triaflatsort(1,2,i) = triaflat(1,2,isource(i))
        triaflatsort(2,2,i) = triaflat(2,2,isource(i))
        triaflatsort(3,2,i) = triaflat(3,2,isource(i))
        triaflatsort(1,3,i) = triaflat(1,3,isource(i))
        triaflatsort(2,3,i) = triaflat(2,3,isource(i))
        triaflatsort(3,3,i) = triaflat(3,3,isource(i))
        trianormsort(1,i) = trianorm(1,isource(i))
        trianormsort(2,i) = trianorm(2,isource(i))
        trianormsort(3,i) = trianorm(3,isource(i))
        enddo
c
        return
        end
c
c
c
c
c
        subroutine h3dzero(mpole,nterms)
        implicit real *8 (a-h,o-z)
c
c       ... set multipole to zero
c
        complex *16 mpole(0:nterms,-nterms:nterms)
c       
        do n=0,nterms
        do m=-n,n
        mpole(n,m)=0
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine h3dmpalloc(wlists,iaddr,nboxes,lmptot,nterms)
        implicit real *8 (a-h,o-z)
        integer box(20)
        dimension nterms(0:1)
        dimension iaddr(2,nboxes)
        dimension center0(3),corners0(3,8)
        dimension wlists(1)
c
c       ... construct pointer array iaddr for addressing multipole and
c       local expansion
c
        iptr=1
        do ibox=1,nboxes
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        level=box(1)
c
c       ... first, allocate memory for the multipole expansion
c       
        iaddr(1,ibox)=iptr
        iptr=iptr+(nterms(level)+1)*(2*nterms(level)+1)*2
c
c       ... then, allocate memory for the local expansion
c       
        iaddr(2,ibox)=iptr
        iptr=iptr+(nterms(level)+1)*(2*nterms(level)+1)*2
c       
        enddo
        lmptot = iptr
        return
        end
c
c
c
c
c
