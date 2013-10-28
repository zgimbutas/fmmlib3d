function [U]=d3tstrcrem(nsource,source,nbox,ifempty,minlev,maxlev,ntarget,target)
%D3TSTRCREM See D3TSTRCR, include empty boxes, min and max level restriction.
%
%  [U]=D3TSTRCREM(NSOURCE,SOURCE,NBOX);
%
%  [U]=D3TSTRCREM(NSOURCE,SOURCE,NBOX,IFEMPTY,MINLEV,MAXLEV);
%
%  [U]=D3TSTRCREM(NSOURCE,SOURCE,NBOX,IFEMPTY,MINLEV,MAXLEV,NTARGET,TARGET);
%
%  This subroutine constructs the logical structure for the 
%  fully adaptive FMM in three dimensions and stores it in the
%  structure U.  It is capable of constructing the oct-tree on 
%  both sources and targets. 
%
%  After that, the user can obtain the information about various boxes 
%  and lists in it by calling D3TGETB, D3TGETL.
%
%  Input parameters:
%  
%  nsource - number of sources
%  source - real (3,nsource) - source locations
%  nbox - maximum number of points in a box on the finest level
%  ntarget - number of targets
%  target - real (3,ntarget) - target locations
%
%  ifempty - ifempty=0 - remove empty boxes, ifempty=1 - keep empty boxes 
%  minlevel - minimum level of refinement
%  maxlevel - minimum level of refinement
%
%  Output parameters:
%
%  U.nboxes - the number of boxes created
%  U.isource - the integer array, addressing the particles in boxes
%       Explanation: for a box ibox, the particles living in it are:
%         (source(1,j),source(2,j),source(3,j)),
%         (source(1,j+1),source(2,j+1),source(3,j+1)),
%         (source(1,j+2),source(2,j+2),source(3,j+3)), ... 
%         (source(1,j+nj-1),source(2,j+nj-1),source(3,j+nj-1)),
%         (source(1,j+nj),source(2,j+nj),source(3,j+nj)),
%         with j=boxes(14,ibox), and nj=boxes(15,ibox)
%  U.nlev - the maximum level number on which any boxes have been created
%  U.laddr - an integer array dimensioned (2,nlev), describing the
%         numbers of boxes on various levels of sybdivision, so that
%         the first box on level (i-1) has sequence number laddr(1,i),
%         and there are laddr(2,i) boxes on level i-1
%  U.center - the center of the box on the level 0, 
%             containing the whole simulation
%  U.size - the side of the box on the level 0
%  U.lists -  the array containing all tables describing boxes, lists, etc. 
%         it is a link-list (for the most part), and can only be accessed
%         via the entries D3TGETB, D3TGETL.
%  U.itarget - the integer array, addressing the target particles in boxes
%       Explanation: for a box ibox, the targets living in it are:
%         (target(1,j),target(2,j),target(3,j)),
%         (target(1,j+1),target(2,j+1),target(3,j+1)),
%         (target(1,j+2),target(2,j+2),target(3,j+3)), ... 
%         (target(1,j+nj-1),target(2,j+nj-1),target(3,j+nj-1)),
%         (target(1,j+nj),target(2,j+nj),target(3,j+nj)),
%         with j=boxes(16,ibox), and nj=boxes(17,ibox)
%  U.lused - the amount of workspace used for storing U.lists
%
%  U.lw - the amount of workspace used for creating U.lists
%  U.ier - error return code
%    ier=0  - successful execution.
%    ier=32 - the amount lw of space in array w is insufficient.
%    ier=16 - the subroutine attempted to construct more 
%        than 199 levels of subdivision; indicates bad trouble.
%    ier=64 - the amount lw of space in array w is severely insufficient.
%

if( nargin == 3 )
  ifempty = 0;
  minlev = 0;
  maxlev = 200;
  ntarget = 0;
  target = zeros(3,1);
  itarget = zeros(1,1);
end

if( nargin == 6 )
  ntarget = 0;
  target = zeros(3,1);
  itarget = zeros(1,1);
end


isource = zeros(1,nsource);
if( ntarget > 0 ), itarget = zeros(1,ntarget); end

nlev = 0;
nboxes = 0;
laddr = zeros(2,400);

center = zeros(3,1);
size = 0;

lw = 100*(nsource+ntarget)+10000;
for i=1:10

w = zeros(1,lw);
lused=0;

ier = 0;

mex_id_ = 'd3tstrcrem(io int[x], i double[], i int[x], io int[x], io int[x], io int[], io int[], io int[x], io double[], io double[], i double[], i int[x], io int[], io double[], i int[x], io int[x], i int[x], i int[x], i int[x])';
[ier, nbox, nboxes, isource, laddr, nlev, center, size, itarget, w, lused] = fmm3d_r2012a(mex_id_, ier, source, nsource, nbox, nboxes, isource, laddr, nlev, center, size, target, ntarget, itarget, w, lw, lused, ifempty, minvel, maxlev, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);

if( ier > 0 ) 
% increase memory allocation
  lw = lw*1.5;
else
  break;
end

end


U.nbox = nbox;
U.ifempty = ifempty;
U.minlev = minlev;
U.maxlev = maxlev;

U.nboxes = nboxes;
U.nlev = nlev;
U.laddr = laddr;
U.center = center;
U.size = size;

if( ier == 0 ), U.lists = w(1,1:lused); end

U.isource = isource;
if( ntarget > 0 ), U.itarget = itarget; end

U.lw = lw;
U.lused = lused;
U.ier=ier;



