function [triaflat,trianorm,triaarea,centroid,triatang1,triatang2]=atriproc(verts,ifaces)
%ATRIPROC Process triangulations in Cart3d format (flat).
%
%  [triaflat]=atriproc(verts,ifaces);
%
%  [triaflat,trianorm,triaarea,centroid]=atriproc(verts,ifaces);
%
%  [triaflat,trianorm,triaarea,centroid,triatang1,triatang2]=
%           atriproc(verts,ifaces);
%
%  Input parameters:
%
%  verts - real(3,nverts): array of triangulation vertices 
%  ifaces - integer(3,nfaces): indices of triangle vertices
%
%  Output parameters:
%
%  triaflat - real(3,3,ntri): array of triangle vertex coordinates 
%  trianorm - real(3,nsource): triangle normals at centroids
%  triaarea - real(nsource): triangle area elements at centroids
%  centroid - real(3,nsource): triangle centroids
%  triatang1 - real(3,nsource): triangle tangents at centroids (first set)
%  triatang2 - real(3,nsource): triangle tangents at centroids (second set)
%
%  Note: the first set of tangent vectors is (\partial xyz/\partial u).
%

%
%  Construct triangle vertex array
%
nverts=size(verts,2);
nfaces=size(ifaces,2);

ntri = nfaces;
triaflat = zeros(3,3,ntri);

%for i=1:ntri	
%  triaflat(1:3,1:3,i) = verts(1:3,ifaces(1:3,i));
%end

triaflat(1:3,1:3,1:ntri) = reshape(verts(1:3,ifaces(1:3,1:ntri)),3,3,ntri);

if( nargout > 1 ),
%
%  Parametrization constants
%
%       ... setup a flat triangle in R^3
%
%
%            2     
%          .   . 
%         .     .
%        .       .
%       0 .. . .. 1
%
%
x0=squeeze(triaflat(1:3,1,:));
x1=squeeze(triaflat(1:3,2,:));
x2=squeeze(triaflat(1:3,3,:));

xu=x1-x0;
xv=x2-x0;

%
%  Triangle centroids
%
u=1/3;
v=1/3;

centroid=x0+u*xu+v*xv;

%  Construct triangle normals
trianorm = zeros(3,ntri);
triaarea = zeros(1,ntri);

vec1 = xu;
vec2 = xv;

trianorm = cross(vec1,vec2);
ds = sqrt(sum(trianorm.^2,1));

trianorm = trianorm ./ repmat(ds,3,1);

%  Triangle area element at the centroid
triaarea = ds/2;
end

if( nargout > 4 ),
%
%  Construct tangent vectors
%
dt = sqrt(sum(vec1.^2,1));
triatang1 = vec1 ./ repmat(dt,3,1);
triatang2 = cross(trianorm,triatang1);
end
