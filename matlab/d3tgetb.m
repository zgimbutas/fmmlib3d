function [ier,box,center,corners]=d3tgetb(ibox,lists)
%D3TGETB Retrieve box information.
%
% [IER,BOX,CENTER,CORNERS]=D3TGETB(IBOX,LISTS);
%
% Input parameters:
%
% ibox - the box number for which the information is desired
% lists - storage area U.lists as created be D3TSTRCR or D3TSTRCREM
%
% Output parameters:
%
% ier - error return code
%    ier=0 - successful execution
%    ier=4 - ibox is either greater than the number of boxes 
%            in the structure or less than 1.
%
% box - an integer array dimensioned box(20). its elements describe 
%        the box number ibox, as follows:
%
%       1. level - the level of subdivision on which this box 
%             was constructed; 
%       2, 3, 4  - the coordinates of this box among  all
%             boxes on this level
%       5 - the daddy of this box, identified by it address
%             in array boxes
%       6,7,8,9,10,11,12,13 - the  list of children of this box 
%             (eight of them, and the child is identified by its address
%             in the array boxes; if a box has only one child, only the
%             first of the four child entries is non-zero, etc.)
%       14 - the location in the array iz of the particles 
%             living in this box
%       15 - the number of particles living in this box
%       16 - the location in the array iztarg of the targets
%             living in this box
%       17 - the number of targets living in this box
%       18 - source box type: 0 - empty, 1 - leaf node, 2 - sub-divided
%       19 - target box type: 0 - empty, 1 - leaf node, 2 - sub-divided
%       20 - reserved for future use
%
% center - real (3) - the center of the box number ibox 
% corners - real (3,8) - the corners of the box number ibox 
%

ier = 0;
center = zeros(3,1);
corners = zeros(3,8);
box = zeros(1,20);

mex_id_ = 'd3tgetb(io int[x], i int[x], io int[], io double[], io double[], i double[])';
[ier, box, center, corners] = fmm3d_r2012a(mex_id_, ier, ibox, box, center, corners, lists, 1, 1);



