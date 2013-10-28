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
%  build FMM tree structures in R^3
%
'build FMM tree in R^3'
tic
nbox = 100;
U = d3tstrcr(nsource,source,nbox);
time_fmmtree=toc


lused=U.lused
nboxes=U.nboxes
nlev=U.nlev
laddr=U.laddr(1:2,1:U.nlev)

