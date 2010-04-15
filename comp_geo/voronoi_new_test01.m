%% VORONOI_NEW_TEST01 compares VORONOI and VORONOI_NEW.
%
%  Modified:
%
%    21 April 2005
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'VORONOI_NEW_TEST01:\n' );
  fprintf ( 1, '  This script compares the behavior of MATLAB''s standard\n' );
  fprintf ( 1, '  "voronoi" command against that of an unreleased version,\n' );
  fprintf ( 1, '  here called "voronoi_new".\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  The problem with the "voronoi" command is that it does\n' );
  fprintf ( 1, '  not return information that can be used to plot the\n' );
  fprintf ( 1, '  semi-infinite sides of the Voronoi regions that lie on\n' );
  fprintf ( 1, '  the boundary.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  This means that a Voronoi diagram drawn from this\n' );
  fprintf ( 1, '  information is incomplete.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  An unreleased version of "voronoi" corrects this\n' );
  fprintf ( 1, '  problem, as demonstrated here.\n' );

  p = [ 0.0  0.0;
        0.0  1.0;
        0.2  0.5;
        0.3  0.6;
        0.4  0.5;
        0.6  0.3;
        0.6  0.5;
        1.0  0.0;
        1.0  1.0 ];
%
%  First plot with VORONOI.
%
  [ vx, vy ] = voronoi ( p(:,1), p(:,2) );

  scatter ( p(:,1), p(:,2), 'b', 'filled' );
  hold on
  plot ( vx, vy, '-', 'LineWidth', 3, 'Color', 'r' );
  axis square
  axis ( [ -1, 2, -1, 2 ] )
  title ( 'Diagram using "Voronoi"' );
  hold off
  
  input ( 'Press return to continue' )
%
%  Now plot with VORONOI_NEW.
%
  [ vx, vy ] = voronoi_new ( p(:,1), p(:,2) );

  scatter ( p(:,1), p(:,2), 'b', 'filled' );
  hold on
  plot ( vx, vy, '-', 'LineWidth', 3, 'Color', 'r' );
  axis square
  axis ( [ -1, 2, -1, 2 ] )
  title ( 'Diagram using "Voronoi\_New"' );
  hold off

  input ( 'Press return to continue' )
