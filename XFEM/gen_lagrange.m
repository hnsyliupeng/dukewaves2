function [ke,id_node,id_lag] =... 
    gen_lagrange(node,x,y,parent,id_eqns,id_dof,...
                 pn_nodes,pos_g,neg_g,numeqns,ex_dof,intersection,endpoints)

% Initialize
ke = zeros(12,2);

xep = [];
yep = [];
xes = [];
yes = [];

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

if all(size(intersection) == [2 2])
    p1 = intersection(1,:);
    p2 = intersection(2,:);
elseif all(size(intersection) == [1 2])
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
      

% jacobian of segment to global
he = sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2);
seg_jcob = he/2;

% Gauss points on segments
gauss = [-sqrt(3)/3 sqrt(3)/3];
weights = [1 1];

N = zeros(2,12);

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
    
        % Evaluate shape function
            N(1,2*b-1) = N(1,2*b-1) + Larea/Area*seg_jcob*weights(g);    % First enrichment
            N(2,2*b)   = N(2,2*b)   + Larea/Area*seg_jcob*weights(g);
            N(1,2*b+5) = N(1,2*b+5) + Larea/Area*seg_jcob*weights(g);    % Second enrichment
            N(2,2*b+6) = N(2,2*b+6) + Larea/Area*seg_jcob*weights(g);
    end
    
end

for c = 1:6
    N(:,2*c-1:2*c) = N(:,2*c-1:2*c)*flg(c);
end
    
ke = N';


% Build id array
nodes = node(:,parent);
id_node(1) = id_eqns(nodes(1),3);  % 1st extra x dof
id_node(2) = id_eqns(nodes(1),4);  % 1st extra y dof
id_node(3) = id_eqns(nodes(2),3);  % 1st extra x dof
id_node(4) = id_eqns(nodes(2),4);  % 1st extra y dof
id_node(5) = id_eqns(nodes(3),3);  % 1st extra x dof
id_node(6) = id_eqns(nodes(3),4);  % 1st extra y dof

id_node(7)  = id_eqns(nodes(1),5);  % 2nd extra x dof
id_node(8)  = id_eqns(nodes(1),6);  % 2nd extra y dof
id_node(9)  = id_eqns(nodes(2),5);  % 2nd extra x dof
id_node(10) = id_eqns(nodes(2),6);  % 2nd extra y dof
id_node(11) = id_eqns(nodes(3),5);  % 2nd extra x dof
id_node(12) = id_eqns(nodes(3),6);  % 2nd extra y dof

id_lag(1) = numeqns + ex_dof - 1;
id_lag(2) = numeqns + ex_dof;



