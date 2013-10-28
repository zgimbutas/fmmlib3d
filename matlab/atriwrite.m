function [nverts,nfaces]=atriwrite(filename,verts,ifaces)
%ATRIWRITE Store triangulations in Cart3d to a file.  (flat)
%
%  [nverts,nfaces]=atriwrite(filename,verts,ifaces);
%
%  Input parameters:
%
%  filename - output file name.
%  verts - real(3,nverts): array of triangulation vertices
%  ifaces - integer(3,nfaces): indices of triangle vertices
%
%  Output parameters:
%
%  nverts - number of triangulation vertices
%  nfaces - number of triangle vertex indices
%

%
%  Store flat triangulation
%

geom_type = 2;

fid = fopen(filename,'w');

nverts=size(verts,2);
nfaces=size(ifaces,2);

fprintf(fid,'%d %d\n',nverts,nfaces);

fprintf(fid,'%20.15e %20.15e %20.15e\n',verts);
fprintf(fid,'%d %d %d\n',ifaces);

fclose(fid);

%%%nverts,nfaces

