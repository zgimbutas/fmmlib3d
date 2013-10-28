function [U]=l3dtriadirect(nsource,triaflat,trianorm,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarget,target,ifpottarg,iffldtarg)
%LFMM3DTRIADIRECT Laplace interactions in R^3, direct evaluation.
%
% [U]=L3DTRIADIRECT(NSOURCE,TRIAFLAT,TRIANORM,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC);
%
% [U]=L3DTRIADIRECT(NSOURCE,TRIAFLAT,TRIANORM,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFPOT,IFFLD);
%
% [U]=L3DTRIADIRECT(NSOURCE,TRIAFLAT,TRIANORM,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFPOT,IFFLD,...
%         NTARGET,TARGET);
%
% [U]=L3DTRIADIRECT(NSOURCE,TRIAFLAT,TRIANORM,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFPOT,IFFLD,...
%         NTARGET,TARGET,IFPOTTARG,IFFLDTARG);
%
%
% This subroutine evaluates the Laplace potential and field due
% to a collection of flat triangles with constant single and/or
% double layer densities. We use (1/r) for the Green's function,
% without the (1/4 pi) scaling.  Self-interactions are included.
%
% It is capable of evaluating the layer potentials either on or 
% off the surface (or both).            
%
%
% Input parameters:
% 
% nsource - number of triangles
% triaflat - real (3,3,nsource): array of triangle vertex coordinates
% trianorm - real (3,nsource): triangle normals
% source - real (3,nsource): triangle centroids
% ifcharge - single layer computation flag
%
%         0 => do not compute
%         1 => include Laplace SLP contribution
% 
% charge - complex (nsource): piecewise constant SLP (charge) strength 
% ifdipole - double layer computation flag
%
%         0 => do not compute
%         1 => include Laplace DLP contribution
% 
% dipole - complex (nsource): piecewise constant DLP (dipole) strength 
% dipvec - real (3,nsource): piecewise constant dipole orientation vectors
%
%        In the present version, dipvec MUST BE SET EQUAL
%        to the triangle normal. It is here as an additional
%        parameter for future use, where an arbitrarily 
%        oriented dipole vector is permitted.
%
% ifpot - potential computation flag, 1 => compute the potential, otherwise no
% iffld - field computation flag, 1 => compute the field, otherwise no
%
% ntarget - number of targets
% target - real (3,ntarget): target locations
%
% ifpottarg - target potential computation flag, 
%         1 => compute the potential, otherwise no
% iffldtarg - target field computation flag, 
%         1 => compute the field, otherwise no
%
%
% Output parameters: 
%
% U.pot - complex (nsource) - potential at triangle centroids
% U.fld - complex (3,nsource) - field (i.e. -gradient) at triangle centroids
% U.pottarg - complex (ntarget) - potential at targets
% U.fldtarg - complex (3,ntarget) - field (i.e. -gradient) at targets
%
% U.ier - error return code
%
%             ier=0     =>  normal execution
%

if( nargin == 9 ) 
  ifpot = 1;
  iffld = 1;
  ntarget = 0;
  target = zeros(3,1);
  ifpottarg = 0;
  iffldtarg = 0;
end

if( nargin == 11 ) 
  ntarget = 0;
  target = zeros(3,1);
  ifpottarg = 0;
  iffldtarg = 0;
end

if( nargin == 13 ) 
  ifpottarg = 1;
  iffldtarg = 1;
end

ifcharge = double(ifcharge); ifdipole = double(ifdipole);
ifpot = double(ifpot); iffld = double(iffld);
ifpottarg = double(ifpottarg); iffldtarg = double(iffldtarg);

pot=0;
fld=zeros(3,1);
pottarg=0;
fldtarg=zeros(3,1);

if( ifpot == 1 ), pot=zeros(1,nsource)+1i*zeros(1,nsource); end;
if( iffld == 1 ), fld=zeros(3,nsource)+1i*zeros(3,nsource); end;
if( ifpottarg == 1 ), pottarg=zeros(1,ntarget)+1i*zeros(1,ntarget); end;
if( iffldtarg == 1 ), fldtarg=zeros(3,ntarget)+1i*zeros(3,ntarget); end;

ier=0;

mex_id_ = 'l3dtriadirect(i int[x], i double[], i double[xx], i double[xx], i int[x], i dcomplex[], i int[x], i dcomplex[], i double[xx], i int[x], io dcomplex[], i int[x], io dcomplex[], i int[x], i double[], i int[x], io dcomplex[], i int[x], io dcomplex[])';
[pot, fld, pottarg, fldtarg] = fmm3d_r2012a(mex_id_, nsource, triaflat, trianorm, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, pot, iffld, fld, ntarget, target, ifpottarg, pottarg, iffldtarg, fldtarg, 1, 3, nsource, 3, nsource, 1, 1, 3, nsource, 1, 1, 1, 1, 1);


if( ifpot == 1 ), U.pot=pot; end
if( iffld == 1 ), U.fld=fld; end
if( ifpottarg == 1 ), U.pottarg=pottarg; end
if( iffldtarg == 1 ), U.fldtarg=fldtarg; end
U.ier=ier;



