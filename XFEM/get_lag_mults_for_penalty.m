% get_lag_mults_for_penalty.m
%
% This method computes the Lagrange multipliers (internal forces at the
% interface) via the penalty method. They are obtained by integrating over
% the gap, so they will be piecewise constant.
%
% Input arguments:
%   node            mapping between elements and nodes
%   x               x-coordinates of all nodes
%   y               y-coordinates of all nodes
%   parent          global element ID of current element
%   id_eqns         mapping between nodes and global DOFs, respecting
%                   enriched nodes
%   id_dof          shows, which base DOFs have enriched DOFs, too.
%   pn_nodes        
%   pos_g           positive grain at this interface
%   neg_g           negative grain at this interface
%   intersection    points on elements edges, that are cut by an interface
%   endpoints       points, that define the interface
%   penalty         penalty paramter
%   fdisp           entire solution vector
%
% Returned variables
%   lagmult         Vector of Lagrange multipliers
%

% Author: Matthias Mayr (04/2010)

function [lagmult] = get_lag_mults_for_penalty(node,x,y,parent,id_eqns, ...
    id_dof,pn_nodes,pos_g,neg_g,intersection,endpoints,penalty,fdisp)

% Initialize
xep = zeros(1,3);
yep = zeros(1,3);
xes = zeros(1,3);
yes = zeros(1,3);

enrich1 = zeros(1,3);
enrich2 = zeros(1,3);

% get nodes of current element 'parent'
nodes = node(:,parent);

% Establish a set of flags

% Get coordinates of parent element
for m=1:3
    jep = node(m,parent); xep(m) = x(jep); yep(m) = y(jep);
end

flg = [0 0 0 0 0 0];

% First enrichment
for n = 1:3     % loop over nodes
    
    % Get the "first" enrichment of node
    
    enrich1(n) = id_dof(nodes(n),3);
    
    if enrich1(n) == pos_g
    
        if (pn_nodes(n,1) == 1)
            flg(n) = 1;
        else
            flg(n) = 0;
        end
      
    elseif enrich1(n) == neg_g
        
        if (pn_nodes(n,2) == 1)
            flg(n) = -1;
        else
            flg(n) = 0;
        end        
    end
end

% Second Enrichment
for n = 1:3     % loop over nodes
    
    % Get the "second" enrichment of nodes
    
    enrich2(n) = id_dof(nodes(n),5);  
    
    if enrich2(n) == pos_g  % If this enrichment corresponds 
                            % to the positive grain
    
        if (pn_nodes(n,1) == 1)
            flg(3 + n) = 1;
        else
            flg(3 + n) = 0;
        end
      
    elseif enrich2(n) == neg_g
        
        if (pn_nodes(n,2) == 1)
            flg(3 + n) = -1;
        else
            flg(3 + n) = 0;
        end        
    end
end


% end points of intersection - direction doesn't matter - this is for the
% segment jacobian calculation

if all(size(intersection) == [2 2])     % no triple junction in the element
    % get the two intersection points
    p1 = intersection(1,:);
    p2 = intersection(2,:);
    
elseif all(size(intersection) == [1 2]) % triple junction in the element
    % get the two points, that define the subsegment
    p1 = intersection(1,:);
    
    % Second endpoint of segment is also end point of interface
    endpoint = endpoints(1,:);
    
    inside = polygon_contains_point_2d ( 3, [xep;yep], endpoint );
    
    if inside
        p2 = endpoint;
    else
        p2 = endpoints(2,:);
    end

end

% initialize
N = zeros(2,12);

% jacobian of segment to global
he = sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2);
seg_jcob = he/2;

% Gauss points on segments
gauss = [-sqrt(3)/3 sqrt(3)/3];
weights = [1 1];

% loop over Gauss points to assemble N
for g = 1:2
    
    % Get real coordinates of gauss points
    xn = 0.5*(1-gauss(g))*p1(1)+0.5*(1+gauss(g))*p2(1);
    yn = 0.5*(1-gauss(g))*p1(2)+0.5*(1+gauss(g))*p2(2);
    

    for b = 1:3     % Evaluate shape functions
        
        % Get coorindates of area opposite node of concern
        for m=1:3
            jes = node(m,parent); xes(m) = x(jes); yes(m) = y(jes);
        end

        xes(b) = xn; yes(b) = yn;

        Area = det([[1 1 1]' xep' yep'])/2;
        Larea = det([[1 1 1]' xes' yes'])/2;
    
        % Evaluate shape function for node 'b'
        N(1,2*b-1) = N(1,2*b-1) + Larea/Area*seg_jcob*weights(g);    % First enrichment
        N(2,2*b)   = N(2,2*b)   + Larea/Area*seg_jcob*weights(g);
        N(1,2*b+5) = N(1,2*b+5) + Larea/Area*seg_jcob*weights(g);    % Second enrichment
        N(2,2*b+6) = N(2,2*b+6) + Larea/Area*seg_jcob*weights(g);
    end
    
end

% set some values to zero depending on, whether they are "positively" or
% "negatively" enriched.
for c = 1:6
    N(:,2*c-1:2*c) = N(:,2*c-1:2*c)*flg(c);
end
% here N contains the evatluated shape functions for two possible
% enrichments. If there is only one enrichment in this element, the
% corresponding entries in N are set to zero by the latter for-loop.

%--------------------------------------------------------------------------

% To decide, if you need one or two enrichments, you cannot look at the
% element, you have to look at each node of the element, because the node
% can be enriched two times, also if this element is only cut by one
% interface. This is the fact, if the double-enriched node belongs to a
% neighboured element, that contains a triple junction.
%
% Check every node, if it is enriched or not. If it is enriched only once, 
% add two entries to 'localdis' and fill up with two zeros. If it is 
% enriched twice, fill in the extra DOFs.

% Check first node
index1 = id_eqns(node(1,parent),3:6);
if all(index1)      % node 1 of element 'parent' is enriched twice
    localdis1 = [fdisp(index1(1:2)) zeros(1,4) fdisp(index1(3:4)) zeros(1,4)];
else
    localdis1 = [fdisp(index1(1:2)) zeros(1,10)];
end;

% Check second node
index2 = id_eqns(node(2,parent),3:6);
if all(index2)      % node 2 of element 'parent' is enriched twice
    localdis2 = [zeros(1,2) fdisp(index2(1:2)) zeros(1,4) ...
        fdisp(index2(3:4)) zeros(1,2)];
else
    localdis2 = [zeros(1,2) fdisp(index2(1:2)) zeros(1,8)];
end;

% Check third node
index3 = id_eqns(node(3,parent),3:6);
if all(index3)      % node 3 of element 'parent' is enriched twice
    localdis3 = [zeros(1,4) fdisp(index3(1:2)) zeros(1,4) ...
        fdisp(index3(3:4))];
else
    localdis3 = [zeros(1,4) fdisp(index3(1:2)) zeros(1,6)];
end;

% Assemble 'localdis'
localdis = localdis1' + localdis2' + localdis3';

%--------------------------------------------------------------------------

% compute lagrange multipliers
lagmult = penalty * N * localdis;

end

