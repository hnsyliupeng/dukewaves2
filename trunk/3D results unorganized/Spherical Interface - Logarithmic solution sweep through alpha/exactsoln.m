function [uexact] = exactsoln(x,y,z,numnod,ls,ifixu,u)

% for nod=1:numnod
%     if (ls(nod)>=0)
%         uexact(nod) = exactsolution(x(nod),y(nod),z(nod));
%     end
%     if (ls(nod)<0)
%         if (ifixu(nod)==1)
%             uexact(nod) = 0;
%         else
%             uexact(nod) = exactsolution(x(nod),y(nod),z(nod));
%         end
%     end
% end

for nod=1:numnod
    if(ifixu(nod)==1)
        uexact(nod) = u(nod);
    else
        uexact(nod) = exactsolution(x(nod),y(nod),z(nod));
    end
end
        