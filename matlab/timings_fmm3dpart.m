%
%  Test Helmholtz and Laplace particle FMMs in R^3
%

nsource = 100000

source = zeros(3,nsource);

  theta=rand(1,nsource)*pi;
  phi=rand(1,nsource)*2*pi;
  source(1,:)=.5*cos(phi).*sin(theta);
  source(2,:)=.5*sin(phi).*sin(theta);
  source(3,:)=.5*cos(theta);


%plot3(source(1,:),source(2,:),source(3,:))

%
%  Timings
%

zk=1;
ifcharge=1;
charge = ones(1,nsource)+1i*zeros(1,nsource);
ifdipole=1;
dipstr = ones(1,nsource)+1i*zeros(1,nsource);
dipvec = ones(3,nsource);

zk
ifcharge
ifdipole
ifpot = 1
iffld = 1

'Helmholtz particle FMM in R^3'

tic
iprec=0
[U]=hfmm3dpart(iprec,zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec);
total_time=toc
speed=nsource/total_time

'Laplace particle FMM in R^3'

tic
iprec=0
[U]=lfmm3dpart(iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec);
total_time=toc
speed=nsource/total_time


zk=1;
ifcharge=1;
charge = ones(1,nsource)+1i*zeros(1,nsource);
ifdipole=0;
dipstr = ones(1,nsource)+1i*zeros(1,nsource);
dipvec = ones(3,nsource);


zk
ifcharge
ifdipole
ifpot = 1
iffld = 0

'Helmholtz particle FMM in R^3'

tic
iprec=0
[U]=hfmm3dpart(iprec,zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec);
total_time=toc
speed=nsource/total_time

'Laplace particle FMM in R^3'

tic
iprec=0
[U]=lfmm3dpart(iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec);
total_time=toc
speed=nsource/total_time


