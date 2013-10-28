function [verts,ifaces,nverts,nfaces]=atriread(filename)
%ATRIREAD Retrieve triangulations in Cart3d from a file.  (flat)
%
%  [verts,ifaces,nverts,nfaces]=atriread(filename);
%
%  Input parameters:
%
%  filename - input file name.
%
%  Output parameters:
%
%  verts - real(3,nverts): array of triangulation vertices
%  ifaces - integer(3,nfaces): indices of triangle vertices
%

%
%  Retrieve flat triangulation
%

geom_type = 2;

fid = fopen(filename,'r');

nverts=0;
nfaces=0;
[nverts] = fscanf(fid,'%d',1);
[nfaces] = fscanf(fid,'%d',1);

verts=zeros(3,nverts);
ifaces=zeros(3,nfaces);

[verts] = fscanf(fid,'%f',[3,nverts]);
[ifaces] = fscanf(fid,'%d',[3,nfaces]);

fclose(fid);

%%%nverts,nfaces

