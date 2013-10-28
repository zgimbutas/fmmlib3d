function [U]=hfmm3dtriampf(iprec,zk,nsource,triaflat,trianorm,source,ifcharge,charge,ifdipole,dipstr,dipvec,ntarget,target,ifpottarg,iffldtarg)
%HFMM3DTRIAMPF Helmholtz triangle FMM in R^3, target evaluation only.
%
% [U]=HFMM3DTRIAMPF(IPREC,ZK,NSOURCE,TRIAFLAT,TRIANORM,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,...
%         NTARGET,TARGET);
%
% [U]=HFMM3DTRIAMPF(IPREC,ZK,NSOURCE,TRIAFLAT,TRIANORM,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,...
%         NTARGET,TARGET,IFPOTTARG,IFFLDTARG);
%
%
% This subroutine evaluates the Helmholtz potential and field due
% to a collection of flat triangles with constant single and/or
% double layer densities. We use (exp(ikr)/r) for the Green's function,
% without the (1/4 pi) scaling. 
%
% It is capable of evaluating the layer potentials off the surface only.
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
% zk - complex Helmholtz parameter
% nsource - number of triangles
% triaflat - real (3,3,nsource): array of triangle vertex coordinates
% trianorm - real (3,nsource): triangle normals
% source - real (3,nsource): triangle centroids
% ifcharge - single layer computation flag
%
%         0 => do not compute
%         1 => include Helmholtz SLP contribution
%         2 => include Helmholtz SLP contribution and subtract Laplace SLP
% 
% charge - complex (nsource): piecewise constant SLP (charge) strength 
% ifdipole - double layer computation flag
%
%         0 => do not compute
%         1 => include Helmholtz DLP contribution
%         2 => include Helmholtz DLP contribution and subtract Laplace DLP
% 
% dipole - complex (nsource): piecewise constant DLP (dipole) strength 
% dipvec - real (3,nsource): piecewise constant dipole orientation vectors
%
%        In the present version, dipvec MUST BE SET EQUAL
%        to the triangle normal. It is here as an additional
%        parameter for future use, where an arbitrarily 
%        oriented dipole vector is permitted.
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

if( nargin == 13 ) 
  ifpottarg = 1;
  iffldtarg = 1;
end

pottarg=0;
fldtarg=zeros(3,1);

if( ifpottarg == 1 ), pottarg=zeros(1,ntarget)+1i*zeros(1,ntarget); end;
if( iffldtarg == 1 ), fldtarg=zeros(3,ntarget)+1i*zeros(3,ntarget); end;

ier=0;

mex_id_ = 'hfmm3dtriampftarg(io int[x], i int[x], i dcomplex[x], i int[x], i double[], i double[xx], i double[xx], i int[x], i dcomplex[], i int[x], i dcomplex[], i double[xx], i int[x], i double[], i int[x], io dcomplex[], i int[x], io dcomplex[])';
[ier, pottarg, fldtarg] = fmm3d_r2012a(mex_id_, ier, iprec, zk, nsource, triaflat, trianorm, source, ifcharge, charge, ifdipole, dipstr, dipvec, ntarget, target, ifpottarg, pottarg, iffldtarg, fldtarg, 1, 1, 1, 1, 3, nsource, 3, nsource, 1, 1, 3, nsource, 1, 1, 1);


if( ifpottarg == 1 ), U.pottarg=pottarg; end
if( iffldtarg == 1 ), U.fldtarg=fldtarg; end
U.ier=ier;


