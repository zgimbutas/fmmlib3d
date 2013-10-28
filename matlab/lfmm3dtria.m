function [U]=lfmm3dtria(iprec,nsource,triaflat,trianorm,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarget,target,ifpottarg,iffldtarg)
%LFMM3DTRIA Laplace triangle FMM in R^3.
%
% [U]=LFMM3DTRIA(IPREC,NSOURCE,TRIAFLAT,TRIANORM,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC);
%
% [U]=LFMM3DTRIA(IPREC,NSOURCE,TRIAFLAT,TRIANORM,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFPOT,IFFLD);
%
% [U]=LFMM3DTRIA(IPREC,NSOURCE,TRIAFLAT,TRIANORM,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFPOT,IFFLD,...
%         NTARGET,TARGET);
%
% [U]=LFMM3DTRIA(IPREC,NSOURCE,TRIAFLAT,TRIANORM,SOURCE,...
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
% iprec - FMM precision flag
%
%             -2 => tolerance =.5d0   =>  
%             -1 => tolerance =.5d-1  =>  1 digit 
%              0 => tolerance =.5d-2  =>  2 digits
%              1 => tolerance =.5d-3  =>  3 digits
%              2 => tolerance =.5d-6  =>  6 digits
%              3 => tolerance =.5d-9  =>  9 digits
%              4 => tolerance =.5d-12 => 12 digits
%              5 => tolerance =.5d-15 => 15 digits
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
%             ier=4     =>  cannot allocate tree workspace
%             ier=8     =>  cannot allocate bulk FMM  workspace
%             ier=16    =>  cannot allocate mpole expansion workspace in FMM
%

if( nargin == 10 ) 
  ifpot = 1;
  iffld = 1;
  ntarget = 0;
  target = zeros(3,1);
  ifpottarg = 0;
  iffldtarg = 0;
end

if( nargin == 12 ) 
  ntarget = 0;
  target = zeros(3,1);
  ifpottarg = 0;
  iffldtarg = 0;
end

if( nargin == 14 ) 
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

if( ntarget == 0 ) 
mex_id_ = 'lfmm3dtriaself(io int[x], i int[x], i int[x], i double[], i double[xx], i double[xx], i int[x], i dcomplex[], i int[x], i dcomplex[], i double[xx], i int[x], io dcomplex[], i int[x], io dcomplex[])';
[ier, pot, fld] = fmm3d_r2012a(mex_id_, ier, iprec, nsource, triaflat, trianorm, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, pot, iffld, fld, 1, 1, 1, 3, nsource, 3, nsource, 1, 1, 3, nsource, 1, 1);
else
mex_id_ = 'lfmm3dtriatarg(io int[x], i int[x], i int[x], i double[], i double[xx], i double[xx], i int[x], i dcomplex[], i int[x], i dcomplex[], i double[xx], i int[x], io dcomplex[], i int[x], io dcomplex[], i int[x], i double[], i int[x], io dcomplex[], i int[x], io dcomplex[])';
[ier, pot, fld, pottarg, fldtarg] = fmm3d_r2012a(mex_id_, ier, iprec, nsource, triaflat, trianorm, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, pot, iffld, fld, ntarget, target, ifpottarg, pottarg, iffldtarg, fldtarg, 1, 1, 1, 3, nsource, 3, nsource, 1, 1, 3, nsource, 1, 1, 1, 1, 1);
end

if( ifpot == 1 ), U.pot=pot; end
if( iffld == 1 ), U.fld=fld; end
if( ifpottarg == 1 ), U.pottarg=pottarg; end
if( iffldtarg == 1 ), U.fldtarg=fldtarg; end
U.ier=ier;


