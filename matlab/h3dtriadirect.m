function [U]=h3dtriadirect(iprec,zk,nsource,triaflat,trianorm,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarget,target,ifpottarg,iffldtarg)
%H3DTRIADIRECT Helmholtz interactions in R^3, direct evaluation.
%
% [U]=H3DTRIADIRECT(IPREC,ZK,NSOURCE,TRIAFLAT,TRIANORM,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC);
%
% [U]=H3DTRIADIRECT(IPREC,ZK,NSOURCE,TRIAFLAT,TRIANORM,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFPOT,IFFLD);
%
% [U]=H3DTRIADIRECT(IPREC,ZK,NSOURCE,TRIAFLAT,TRIANORM,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFPOT,IFFLD,...
%         NTARGET,TARGET);
%
% [U]=H3DTRIADIRECT(IPREC,ZK,NSOURCE,TRIAFLAT,TRIANORM,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFPOT,IFFLD,...
%         NTARGET,TARGET,IFPOTTARG,IFFLDTARG);
%
%
% This subroutine evaluates the Helmholtz potential and field due
% to a collection of flat triangles with constant single and/or
% double layer densities. We use (exp(ikr)/r) for the Green's function,
% without the (1/4 pi) scaling.  Self-interactions are included.
%
% It is capable of evaluating the layer potentials either on or 
% off the surface (or both).            
%
%
% Input parameters:
% 
% iprec - precision flag, as used in FMM
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

if( nargin == 11 ) 
  ifpot = 1;
  iffld = 1;
  ntarget = 0;
  target = zeros(3,1);
  ifpottarg = 0;
  iffldtarg = 0;
end

if( nargin == 13 ) 
  ntarget = 0;
  target = zeros(3,1);
  ifpottarg = 0;
  iffldtarg = 0;
end

if( nargin == 15 ) 
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

%
% The following assumes reasonably regular triangulation.
%
%       nqtri: integer: number of quadrature nodes.
%          Suggested values for nqtri
%          iprec:  FMM precision flag
%                 -2 => tolerance =.5d0
%                 -1 => tolerance =.5d-1
%                  0 => tolerance =.5d-2
%                  1 => tolerance =.5d-3
%                  2 => tolerance =.5d-6
%                  3 => tolerance =.5d-9
%                  4 => tolerance =.5d-12
%                  5 => tolerance =.5d-15
%          if( iprec .eq. -2 ) nqtri=1
%          if( iprec .eq. -1 ) nqtri=2
%          if( iprec .eq.  0 ) nqtri=3
%          if( iprec .ge.  1 ) nqtri=6
%
nqtri=1;
if( iprec == -2), nqtri = 1; end;
if( iprec == -1), nqtri = 2; end;
if( iprec ==  0), nqtri = 3; end;
if( iprec >=  1), nqtri = 6; end;


mex_id_ = 'h3dtriadirect(i int[x], i dcomplex[x], i int[x], i double[], i double[xx], i double[xx], i int[x], i dcomplex[], i int[x], i dcomplex[], i double[xx], i int[x], io dcomplex[], i int[x], io dcomplex[], i int[x], i double[], i int[x], io dcomplex[], i int[x], io dcomplex[])';
[pot, fld, pottarg, fldtarg] = fmm3d(mex_id_, nqtri, zk, nsource, triaflat, trianorm, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, pot, iffld, fld, ntarget, target, ifpottarg, pottarg, iffldtarg, fldtarg, 1, 1, 1, 3, nsource, 3, nsource, 1, 1, 3, nsource, 1, 1, 1, 1, 1);


if( ifpot == 1 ), U.pot=pot; end
if( iffld == 1 ), U.fld=fld; end
if( ifpottarg == 1 ), U.pottarg=pottarg; end
if( iffldtarg == 1 ), U.fldtarg=fldtarg; end
U.ier=ier;


