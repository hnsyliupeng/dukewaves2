function value = tetra_unit_volume ( DUMMY )

%% TETRA_UNIT_VOLUME returns the volume of the unit tetrahedron in 3D.
%
%  Integration region:
%
%    Points (X,Y,Z) such that
%
%      0 <= X,
%      0 <= Y,
%      0 <= Z, and
%      X + Y + Z <= 1.
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
%    Input, double precision DUMMY, a dummy argument, because MATLAB won't
%    allow functions with no arguments.
%
%    Output, double precision TETRA_UNIT_VOLUME, the volume of the unit
%    tetrahedron.
%
  value = 1.0 / 6.0;
