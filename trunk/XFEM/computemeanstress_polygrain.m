% computemeanstress_polygrain.m
%
% CALL:
%
% Computes the mean stress as in 'Kanit2003' eq. (3).
%
% Input parameters:
%
%
% Returned variables:
%   meanstress        vector, containing mean stress in xx, yy, xy and
%                     von-Mises
%

% Author: Matthias Mayr (09/2010)

function [meanstress] = computemeanstress_polygrain(stress, ...
  stressvonmises,x,y,node,X,Y,CONN,SUBELEM_INFO, ...
  SUBELEMENT_GRAIN_MAP,cutlist,elemgrainmap)

%% INITIALIZE
numele = size(node,2);  % number of elments
area_tot = 0;           % total area
ms_xx = 0;  % mean stress in xx
ms_yy = 0;  % mean stress in yy
ms_xy = 0;  % mean stress in xy
ms_vM = 0;  % mean stress von-Mises
% ----------------------------------------------------------------------- %
%% LOOP OVER ELEMENTS
for e=1:numele
  if cutlist(e) == 0
    % element is uncut
    % get some element data
    elenodes = node(:,e);
    xcoords = x(elenodes);
    ycoords = y(elenodes);
    area_ele = det([[1 1 1]' xcoords' ycoords'])/2;
    grain = SUBELEMENT_GRAIN_MAP(e);
    
    % get stress in current element
    s_xx = stress(e,4,grain);
    s_yy = stress(e,5,grain);
    s_xy = stress(e,6,grain);
    s_vM = stressvonmises(e,4,grain);
    
    % update meanstresses
    ms_xx = ms_xx + s_xx * area_ele;
    ms_yy = ms_yy + s_yy * area_ele;
    ms_xy = ms_xy + s_xy * area_ele;
    ms_vM = ms_vM + s_vM * area_ele;
    
    % add element area to total area
    area_tot = area_tot + area_ele;
  else
    % element is cut
    % loop over subelements
    for s=1:SUBELEM_INFO(e).no_kids
      % get some data from subelement
      sID = SUBELEM_INFO(e).kids(s);
      subelenodes = CONN(:,sID);
      Xcoords = X(subelenodes);
      Ycoords = Y(subelenodes);
      area_subele = det([[1 1 1]' Xcoords' Ycoords'])/2;
      grain = SUBELEMENT_GRAIN_MAP(sID);
      
      % get stress in current element
      s_xx = stress(e,4,grain);
      s_yy = stress(e,5,grain);
      s_xy = stress(e,6,grain);
      s_vM = stressvonmises(e,4,grain);

      % update meanstresses
      ms_xx = ms_xx + s_xx * area_subele;
      ms_yy = ms_yy + s_yy * area_subele;
      ms_xy = ms_xy + s_xy * area_subele;
      ms_vM = ms_vM + s_vM * area_subele;

      % add element area to total area
      area_tot = area_tot + area_subele;
    end;
  end;
end;

% compute mean stress vector
meanstress = [ms_xx ms_yy ms_xy ms_vM] / area_tot;
% ----------------------------------------------------------------------- %
end