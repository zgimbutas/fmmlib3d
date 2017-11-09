cc Copyright (C) 2009-2012: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This software is being released under a modified FreeBSD license
cc (see COPYING in home directory). 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date: 2011-05-24 20:42:58 -0400 (Tue, 24 May 2011) $
c    $Revision: 2001 $
c
c       
c        this file contains various FMM tree plotting routines
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine plot_points3d(iw,z,n)
        implicit real *8 (a-h,o-z)
        dimension z(3,1)
c
c       write coordinates of n points to file with uit number iw.
c
        do i=1,n
        write(iw,1000) z(1,i),z(2,i),z(3,i)
        enddo
 1000   format(6(1x,e11.5))
        return
        end
c
c
c
        subroutine plot_box3d(iw,center,size)
        implicit real *8 (a-h,o-z)
        dimension center(3)
c
c       write out data to plot box boundaries.
c       
        write(iw,1000) 
     $     center(1)-size/2,center(2)-size/2,center(3)-size/2
        write(iw,1000) 
     $     center(1)-size/2,center(2)+size/2,center(3)-size/2
        write(iw,1000) 
     $     center(1)+size/2,center(2)+size/2,center(3)-size/2
        write(iw,1000) 
     $     center(1)+size/2,center(2)-size/2,center(3)-size/2
        write(iw,1000) 
     $     center(1)-size/2,center(2)-size/2,center(3)-size/2
        write(iw,1200)
c       
        write(iw,1000) 
     $     center(1)-size/2,center(2)-size/2,center(3)-size/2
        write(iw,1000) 
     $     center(1)-size/2,center(2)-size/2,center(3)+size/2
        write(iw,1200) 
c       
        write(iw,1000) 
     $     center(1)-size/2,center(2)+size/2,center(3)-size/2
        write(iw,1000) 
     $     center(1)-size/2,center(2)+size/2,center(3)+size/2
        write(iw,1200) 
c       
        write(iw,1000) 
     $     center(1)+size/2,center(2)+size/2,center(3)-size/2
        write(iw,1000) 
     $     center(1)+size/2,center(2)+size/2,center(3)+size/2
        write(iw,1200) 
c       
        write(iw,1000) 
     $     center(1)+size/2,center(2)-size/2,center(3)-size/2
        write(iw,1000) 
     $     center(1)+size/2,center(2)-size/2,center(3)+size/2
        write(iw,1200) 
c
        write(iw,1000) 
     $     center(1)-size/2,center(2)-size/2,center(3)+size/2
        write(iw,1000) 
     $     center(1)-size/2,center(2)+size/2,center(3)+size/2
        write(iw,1000) 
     $     center(1)+size/2,center(2)+size/2,center(3)+size/2
        write(iw,1000) 
     $     center(1)+size/2,center(2)-size/2,center(3)+size/2
        write(iw,1000) 
     $     center(1)-size/2,center(2)-size/2,center(3)+size/2
        write(iw,1200)
        write(iw,1000) 
c       
 1000   format(6(1x,e11.5))
 1200   format(80a1)
        return
        end 
c
c
c
        subroutine plot_label3d(iw,center,size,itag,label)
        implicit real *8 (a-h,o-z)
        dimension center(3)
c       
c        write(iw,*) 'set label ', itag, ' "', label, 
c     $     '" at ', center(1), ' , ', center(2), ' , ', 
c     $     center(3), ' center'  
c       
        write(iw,*) 'set label ', itag, ' sprintf("%d",', label, 
     $     ') at ', center(1), ' , ', center(2), ' , ', 
     $     center(3), ' center'  
c
        return
        end
c
