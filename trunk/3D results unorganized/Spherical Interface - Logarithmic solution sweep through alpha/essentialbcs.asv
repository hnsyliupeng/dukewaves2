function [ifixu, u] = essentialbcs(x,y,z,ls,numnod)
ifixu = zeros(numnod,1);
u = zeros(numnod,1);
tol = 1.0e-13;
for nod=1:numnod
    %%%% Top Surface %%%%%%%
    if(abs(z(nod)-1)<tol)
        [uex] = exactsolution(x(nod),y(nod),z(nod));
        u(nod) = uex;      
        ifixu(nod) = 1;
    end
%      %%%% Bottom Surface %%%%%
%     if(z(nod)==0)
%         [uex] = exactsolution(x(nod),y(nod),z(nod));
%         u(nod) = uex;
%         ifixu(nod) = 1;
%     end
%     %%%% Left Surface %%%%%
%     if (y(nod)==0)
%         [uex] = exactsolution(x(nod),y(nod),z(nod));
%         u(nod) = uex;
%         ifixu(nod) = 1;
%     end
%     %%%% Right Surface %%%%
    if (abs(y(nod)-1)<tol)
        [uex] = exactsolution(x(nod),y(nod),z(nod));
        u(nod) = uex;
        ifixu(nod) = 1;
    end
    %%%% Back Surface %%%%%
    if (abs(x(nod)-1)<tol)
        [uex] = exactsolution(x(nod),y(nod),z(nod));
        u(nod) = uex;
        ifixu(nod) = 1;
    end
%     %%%% Front Surface %%%%%
%     if (x(nod)==0)
%         [uex] = exactsolution(x(nod),y(nod),z(nod));
%         u(nod) = uex;
%         ifixu(nod) = 1;
%     end
%     %%% Set all the nodes in the negative region to zero %%%
%     if (ls(nod)<0)
%         ifixu(nod) = 1;
%         u(nod) = 0;
%     end
end