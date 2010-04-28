% colorelements.m
%
% Function call:    colorelements(node,x,y,eleIDs,color)
%
% This method colors the elements, listed in the vector 'eleIDs', with the
% color, given in 'color'.
%
% Input parameters
%   node        mapping between elements and nodes
%   x           x-coordinates of all nodes
%   y           y-coordinates of all nodes
%   eleIDs      vector with global element IDs of those elements, that
%               shall be colored
%   color       color for elements (same syntax as for plot-command)
%               example 1: 'g' for green
%               example 2: [1 0 0] for red via RGB-scheme
%
% Returned variables
%   none
%

% Author: Matthias Mayr (04/2010)

function [] = colorelements(node,x,y,eleIDs,color)

hold on;

if ischar(color)
    patch(x(node(:,eleIDs)),y(node(:,eleIDs)),color);
elseif isnumeric
    patch(x(node(:,eleIDs)),y(node(:,eleIDs)),'Color',color);
end

