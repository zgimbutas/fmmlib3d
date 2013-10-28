%
%  Test Laplace triangle FMMs in R^3
%

%
%  Retrieve flat triangulation
%

geom_type = 2;
filename_geo = 'sphere180.a.tri';
filename_geo = 'sphere720.a.tri';
%filename_geo = 'sphere2880.a.tri';
%filename_geo = 'sphere11520.a.tri';
%filename_geo = 'sphere20480.a.tri';

fid = fopen(filename_geo,'r');

nverts=0;
nfaces=0;
[nverts] = fscanf(fid,'%d',1);
[nfaces] = fscanf(fid,'%d',1);

verts=zeros(3,nverts);
ifaces=zeros(3,nfaces);

[verts] = fscanf(fid,'%f',[3,nverts]);
[ifaces] = fscanf(fid,'%d',[3,nfaces]);

fclose(fid);


[triaflat,trianorm,triaarea,centroid]=atriproc(verts,ifaces);

triangles = triaflat;
dipvec = trianorm;
source = centroid;

ntri = nfaces;
nsource = ntri

%
%  timings
%

ifcharge=1;
charge = rand(3,ntri);
ifdipole=0;
dipstr = rand(3,ntri);
dipvec = trianorm;

%%dipstr = zeros(3,ntri);
%%dipstr(2,:) = 1;
%%dipstr(2,:) = 1*cos(.1*source(3,:));
%%dipstr = cross(dipstr,cross(dipstr,trianorm));

ifcharge
ifdipole
ifpot = 1
iffld = 1

ntarget = ntri
target = source(:,1:ntri);
target(1,:) = target(1,:) + 10;

[ndim,ntarget] = size(target)
ifpottarg = 1
iffldtarg = 0
ntarget
%ntarget = 0


disp('')
'Laplace triangle target FMM in R^3'

tic
iprec=0
[U]=lfmm3dtria(iprec,ntri,triangles,trianorm,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarget,target,ifpottarg,iffldtarg);
total_time=toc


'Laplace triangle target direct evaluation in R^3'
tic
[F]=l3dtriadirect(ntri,triangles,trianorm,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarget,target,ifpottarg,iffldtarg);
total_time=toc


if( ifpot ), U.pot=U.pot/(4*pi); end
if( iffld ), U.fld=U.fld/(4*pi); end

if( ifpot ), F.pot=F.pot/(4*pi); end
if( iffld ), F.fld=F.fld/(4*pi); end

if( ifpot ),
%rms_pot = norm((F.pot),2)/sqrt(nsource)
rms_error_pot = norm((U.pot - F.pot),2)/sqrt(nsource)
rel_error_pot = norm((U.pot - F.pot),2)/norm((F.pot),2)
end

if( iffld ),
%rms_fld = norm(reshape(F.fld,9,nsource),2)/sqrt(nsource)
rms_error_fld = norm(U.fld - F.fld,2)/sqrt(nsource)
rel_error_fld = norm(U.fld - F.fld,2)/norm(F.fld,2)
end
%%%break;

if( ifpottarg ), U.pottarg=U.pottarg/(4*pi); end
if( iffldtarg ), U.fldtarg=U.fldtarg/(4*pi); end

if( ifpottarg ), F.pottarg=F.pottarg/(4*pi); end
if( iffldtarg ), F.fldtarg=F.fldtarg/(4*pi); end

if( ifpottarg ),
%rms_pottarg = norm((F.pottarg),2)/sqrt(nsource)
rms_error_pottarg = norm((U.pottarg - F.pottarg),2)/sqrt(ntarget)
norm_pottarg = norm((F.pottarg),2)
rel_error_pottarg = norm((U.pottarg - F.pottarg),2)/norm((F.pottarg),2)
end

if( iffldtarg ),
rms_fldtarg = norm(F.fldtarg,2)/sqrt(ntarget)
rms_error_fldtarg = ...
    norm(U.fldtarg - F.fldtarg,2)/sqrt(ntarget)
rel_error_fldtarg = ...
    norm(U.fldtarg - F.fldtarg,2)/ ...
    norm(F.fldtarg,2)
end
%%%break;


disp('')
'Laplace triangle FMM in R^3'

ifpottarg = 0
iffldtarg = 0


tic
iprec=0
[U]=lfmm3dtria(iprec,ntri,triangles,trianorm,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld);
total_time=toc

'Laplace triangle direct evaluation in R^3'

tic
[F]=l3dtriadirect(ntri,triangles,trianorm,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld);
total_time=toc



if( ifpot ), U.pot=U.pot/(4*pi); end
if( iffld ), U.fld=U.fld/(4*pi); end

if( ifpot ), F.pot=F.pot/(4*pi); end
if( iffld ), F.fld=F.fld/(4*pi); end

if( ifpot ),
%rms_pot = norm((F.pot),2)/sqrt(nsource)
rms_error_pot = norm((U.pot - F.pot),2)/sqrt(nsource)
end

if( iffld ),
%rms_fld = norm(F.fld,2)/sqrt(nsource)
rms_error_fld = norm(U.fld - F.fld,2)/sqrt(nsource)
end
%%%break;

