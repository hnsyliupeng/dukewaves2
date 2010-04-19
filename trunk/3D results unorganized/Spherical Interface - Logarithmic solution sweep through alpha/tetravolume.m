function [value] = tetravolume(0, 0, 1)

%% TETRA_VOLUME computes the volume of a tetrahedron in 3D.
%
%  Integration region:
%
%    Points inside a tetrahedron whose four vertices are given.
%
%  Modified:
%
%    26 May 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, double precision X(4), Y(4), Z(4), the vertices.
%
%    Output, double precision TETRA_VOLUME, the volume of the tetrahedron.
%
  value = tetraunitvolume('DUMMY')*parallelipipedvolume3d(x, y, z);
