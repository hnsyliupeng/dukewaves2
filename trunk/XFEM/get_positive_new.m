function [pn_nodes] = get_positive_new(parent,pos_g,neg_g)

global NODAL_ENRICH 
global CONN

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
