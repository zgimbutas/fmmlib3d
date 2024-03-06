% test_mesh3.m
%   Test MWrap interface to mesh data structure.
%
% Copyright (c) 2007  David Bindel
% See the file COPYING for copying permissions

init;

mobj = Mesh_create(2, 2, 4);
try
  fprintf('-- Four element elastic patch test --\n');

  % Nodes:
  %       1     2     3     4     5     6     7     8     9
  x = [ 0.0,  4.0, 10.0,  0.0,  5.5, 10.0,  0.0,  4.2, 10.0;
        0.0,  0.0,  0.0,  4.5,  5.5,  5.0, 10.0, 10.0, 10.0];

  % Boundary codes and values
  bc =[   1,    0,    0,    1,    0,    0,    1,    0,    0;
          1,    0,    0,    0,    0,    0,    0,    0,    0];
  bv =[   0,    0,  2.5,    0,    0,  5.0,    0,    0,  2.5,
          0,    0,    0,    0,    0,    0,    0,    0,    0];

  % Elements:
  %      1  2  3  4
  ix = [ 1, 2, 4, 5 ;
         2, 3, 5, 6 ;
         5, 6, 8, 9 ;
         4, 5, 7, 8 ];

  % Set up problem mesh
  material = Mesh_add_elastic2d(mobj, 1000.0, 0.25, 'plane strain');
  Mesh_load(mobj, material, x, ix);
  Mesh_initialize(mobj);

  % Set boundary conditions
  Mesh_set_bc(mobj, bc);
  Mesh_set_bv(mobj, bv);
  Mesh_assign_ids(mobj);
  
  % Assemble stiffness and force
  K = Mesh_assemble_K(mobj);
  F = Mesh_assemble_F(mobj);

  % Solve for the reduced displacement
  u = -K\F;

  % Patch test should recover linear fields -- check that it does
  uu = Mesh_u(mobj);
  resid = uu - (uu/x)*x;
  fprintf('Patch test residual: %g\n', norm(resid));

  % Assemble residual and make sure it's zero
  Mesh_set_ur(mobj, u);
  RR = Mesh_assemble_F(mobj);
  fprintf('Reported residual force: %g\n', norm(RR));

catch
  
  fprintf('Error: %s\n', lasterr);
  
end
Mesh_delete(mobj);
