function atriplot(verts,ifaces)
%ATRIPLOT Plot flat triangulation.
%
%  atriplot(verts,ifaces)
%
%  Input parameters:
%
%  verts - real(3,nverts): array of triangulation vertices (flat)
%  ifaces - integer(3,nfaces): indices of triangle vertices (flat)
%

%%%trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)')

nfaces=size(ifaces,2);
C=zeros(nfaces,1);
trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)',C')

axis equal
axis auto
%%%alpha(.5)
