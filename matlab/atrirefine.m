function [verts1,ifaces1]=atrirefine(verts,ifaces)
%ATRIREFINE Refine triangulations in Cart3d format (flat).
%
%  [verts1,ifaces1]=atrirefine(verts,ifaces);
%
%
%  Input parameters:
%
%  verts - real(3,nverts): array of triangulation vertices 
%  ifaces - integer(3,nfaces): indices of triangle vertices
%
%  Output parameters:
%
%  verts1 - real(3,nverts1): array of triangulation vertices (refined)
%  ifaces1 - integer(3,nfaces1): indices of triangle vertices (refined)
%

noversamp = 2;


nverts = size(verts,2);

verts1 = verts;

verts1 = [verts1 (verts(:,ifaces(1,:))+verts(:,ifaces(2,:)))/2];
verts1 = [verts1 (verts(:,ifaces(2,:))+verts(:,ifaces(3,:)))/2];
verts1 = [verts1 (verts(:,ifaces(3,:))+verts(:,ifaces(1,:)))/2];


nfaces = size(ifaces,2);

ifaces2 = [(nverts+1:nverts+nfaces); (nverts+nfaces+1:nverts+nfaces*2); (nverts+nfaces*2+1:nverts+nfaces*3)];

ifaces1 = [[ifaces(1,:); ifaces2(1,:); ifaces2(3,:)] ...
	   [ifaces(2,:); ifaces2(2,:); ifaces2(1,:)] ...
	   [ifaces(3,:); ifaces2(3,:); ifaces2(2,:)] ...
	   [ifaces2(1,:); ifaces2(2,:); ifaces2(3,:)]];