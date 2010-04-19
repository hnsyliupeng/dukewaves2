function [volumep,kenit,fe_nit,kepen,fe_pen] = volcalc(node,e,x,y,z)

xe = x(node(1:4,e));
ye = y(node(1:4,e));
ze = z(node(1:4,e));
nlink = 4;
% if (ls(node(:,e))<0)        %% If level set value is less than zero, volume on the positive side is zero
%     volumep = 0;
%     kenit = zeros(nlink,nlink);
%     kepen = zeros(nlink,nlink);
%     fe_nit = zeros(nlink,1);
%     fe_pen = zeros(nlink,1);
% else
    volumep = tetraunitvolume('DUMMY')*parallelipipedvolume3d(xe, ye, ze);
    kenit = zeros(nlink,nlink);
    kepen = zeros(nlink,nlink);
    fe_nit = zeros(nlink,1);
    fe_pen = zeros(nlink,1);
% end
