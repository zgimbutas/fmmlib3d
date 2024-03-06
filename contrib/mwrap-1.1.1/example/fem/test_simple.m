% test_mesh.m
%   Test MWrap interface to mesh data structure.
%
% Copyright (c) 2007  David Bindel
% See the file COPYING for copying permissions

init;

mobj = Mesh_create(1, 1, 2);
try

  % -- Set up a simple mesh
  i1 = Mesh_add_node(mobj, 0);
  i2 = Mesh_add_node(mobj, 1);
  i3 = Mesh_add_node(mobj, 2);
  i4 = Mesh_add_node(mobj, 3);
  s1 = Mesh_add_scalar1d(mobj, 1);
  e1 = Mesh_add_element(mobj, s1, [i1, i2]);
  e2 = Mesh_add_element(mobj, s1, [i2, i3]);
  e3 = Mesh_add_element(mobj, s1, [i3, i4]);
  Mesh_initialize(mobj);

  % -- Assign Dirichlet BC to first dof at node 1
  Mesh_set_bc(mobj, 1,  1, 1);

  fprintf('-- Mesh nodes and connectivities --\n');
  x  = Mesh_x(mobj)
  ix = Mesh_ix(mobj)

  fprintf('-- Mesh dof assignments --\n');
  numid = Mesh_initialize(mobj)
  id = Mesh_id(mobj)

  fprintf('-- Assembled stiffness (standard three-point stencil) --\n');
  K = Mesh_assemble_K(mobj);
  K = full(K)
  
catch
  
  fprintf('Error: %s\n', lasterr);
  
end
Mesh_delete(mobj);
