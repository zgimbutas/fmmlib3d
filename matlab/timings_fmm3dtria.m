%
%  Test Helmholtz and Laplace triangle FMMs in R^3
%

%
%  Retrieve flat triangulation
%

geom_type = 2;
filename_geo = 'sphere320.a.tri';
filename_geo = 'sphere1280.a.tri';
filename_geo = 'psm1.a.tri';
%%%filename_geo = 'psm2.a.tri';
%%%filename_geo = 'sphere11520.a.tri';
%%%filename_geo = 'cube768.a.tri';
%%%filename_geo = 'cube3072.a.tri';

[verts,ifaces,nverts,nfaces] = atriread(filename_geo);
nverts,nfaces

%
%  Construct triangle vertex and normal arrays, triangle area vector, 
%  and source locations at triangle centroids
%

ntri = nfaces;
[triangles,trianorm,triaarea,source]=atriproc(verts,ifaces);

%plot3(source(1,:),source(2,:),source(3,:))

%
%  Timings
%

zk=1;
ifcharge=1;
charge = ones(1,ntri)+1i*zeros(1,ntri);
ifdipole=1;
dipstr = ones(1,ntri)+1i*zeros(1,ntri);
dipvec = trianorm;

zk
ifcharge
ifdipole
ifpot = 1
iffld = 1

'Helmholtz triangle FMM in R^3'

tic
iprec=0
U=hfmm3dtrif(iprec,zk,ntri,triangles,trianorm,source,ifcharge,charge,ifdipole,dipstr,dipvec);
total_time=toc
speed=ntri/total_time

'Laplace triangle FMM in R^3'

tic
iprec=0
U=lfmm3dtria(iprec,ntri,triangles,trianorm,source,ifcharge,charge,ifdipole,dipstr,dipvec);
total_time=toc
speed=ntri/total_time

zk=1;
ifcharge=1;
charge = ones(1,ntri)+1i*zeros(1,ntri);
ifdipole=0;
dipstr = ones(1,ntri)+1i*zeros(1,ntri);
dipvec = trianorm;

zk
ifcharge
ifdipole
ifpot = 1
iffld = 0

'Helmholtz triangle FMM in R^3'

tic
iprec=0
U=hfmm3dtrif(iprec,zk,ntri,triangles,trianorm,source,ifcharge,charge,ifdipole,dipstr,dipvec);
total_time=toc
speed=ntri/total_time

'Laplace triangle FMM in R^3'

tic
iprec=0
U=lfmm3dtria(iprec,ntri,triangles,trianorm,source,ifcharge,charge,ifdipole,dipstr,dipvec);
total_time=toc
speed=ntri/total_time

