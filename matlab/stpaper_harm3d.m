%
%  Test Laplace particle FMMs in R^3
%

nsources = [5000 5000 50000 500000]
iprecs = [-1 0 1 2]

for i=1:size(nsources,2)
    
nsource = nsources(i);

source = zeros(3,nsource);

idist=3;

if( idist == 1 ),
theta=rand(1,nsource)*pi;
phi=rand(1,nsource)*2*pi;
source(1,:)=.5*cos(phi).*sin(theta);
source(2,:)=.5*sin(phi).*sin(theta);
source(3,:)=.5*cos(theta);
end

if( idist == 2 ),
source(1,:)=rand(1,nsource);
source(2,:)=rand(1,nsource);
source(3,:)=rand(1,nsource);
end

if( idist == 3 ),
phi=rand(1,nsource)*2*pi;
source(1,:)=.5*cos(phi);
source(2,:)=.5*sin(phi);
source(3,:)=rand(1,nsource);
end

%scatter3(source(1,:),source(2,:),source(3,:))
%plot3(source(1,:),source(2,:),source(3,:))

for j=1:size(iprecs,2)

iprec = iprecs(j);

%
%  Timings
%

ifcharge=1;
charge = ones(1,nsource)+1i*zeros(1,nsource);
ifdipole=0;
dipstr = ones(1,nsource)+1i*zeros(1,nsource);
dipvec = ones(3,nsource);

ifcharge;
ifdipole;
ifpot = 1;
iffld = 1;

'Laplace particle FMM in R^3';

tic
[U]=lfmm3dpart(iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld);
total_time=toc;
total_time_fmm=total_time;
speed=nsource/total_time;

'Laplace particle direct evaluation in R^3';

k = 100;
if( nsource < 10000), k = 1; end;
if( nsource > 100000), k = 1000; end;

ms = nsource / k;
tic
[F]=l3dpartdirecttime(ms,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld);
total_time=toc;
total_time_direct=total_time*k;
speed=nsource/total_time;


if( ifpot ), U.pot=U.pot/(4*pi); end
if( iffld ), U.fld=U.fld/(4*pi); end

if( ifpot ), F.pot=F.pot/(4*pi); end
if( iffld ), F.fld=F.fld/(4*pi); end



if( ifpot ), U.pot=U.pot(1:ms); end
if( iffld ), U.fld=U.fld(1:ms); end

if( ifpot ), F.pot=F.pot(1:ms); end
if( iffld ), F.fld=F.fld(1:ms); end


if( ifpot ),
rel_error_pot = norm((U.pot - F.pot),2)/norm((F.pot),2);
end

if( iffld ),
rel_error_fld = norm(U.fld - F.fld,2)/norm(F.fld,2);
end


  fprintf('%5d & %2d & %5.3f & %5.3f & %6.3e \\\\ \n', nsource, iprec, total_time_fmm, total_time_direct, rel_error_pot );


end
end
