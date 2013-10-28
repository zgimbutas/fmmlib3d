function fmm3dprini(unit1,unit2)
%FMM3DPRINI Initialize simple printing routines.
%
% Calling FMM3DPRINI(6,13) causes printing to screen and file fort.13.     
%

if (nargin == 1 )
unit2=0;
end

mex_id_ = 'prini(i int[x], i int[x])';
fmm3d_r2012a(mex_id_, unit1, unit2, 1, 1);


