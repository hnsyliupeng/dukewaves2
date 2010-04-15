function [pos_g,neg_g,pn_nodes] = get_positive(i,j,nodegrainmap,p)

global NODAL_ENRICH
global SUBELEMENT_GRAIN_MAP
global INT_INTERFACE 
global CONN X Y

% "get" first subelement associated with this subsegment

parent = i;

sub1 = INT_INTERFACE(i).pairings(j,1);
sub2 = INT_INTERFACE(i).pairings(j,2);

grain_1 = SUBELEMENT_GRAIN_MAP(sub1);
grain_2 = SUBELEMENT_GRAIN_MAP(sub2);

p1 = [INT_INTERFACE(i).shared(2*j-1,1) INT_INTERFACE(i).shared(2*j-1,2)];
p2 = [INT_INTERFACE(i).shared(2*j,1) INT_INTERFACE(i).shared(2*j,2)];

center_of_grain1 = p(grain_1,:);

if halfplane_contains_point_2d(p1,p2,center_of_grain1)
    pos_sub = sub2;
    neg_sub = sub1;
else
    pos_sub = sub1;
    neg_sub = sub2;
end

pos_g = SUBELEMENT_GRAIN_MAP(pos_sub);
neg_g = SUBELEMENT_GRAIN_MAP(neg_sub);

nodes = CONN(:,parent);

pn_nodes = zeros(3,2);

for m = 1:3
    if find(NODAL_ENRICH(nodes(m)).enrichment == pos_g)
        pn_nodes(m,1) = 1;
    end
    if find(NODAL_ENRICH(nodes(m)).enrichment == neg_g)
        pn_nodes(m,2) = 1;
    end        
end
