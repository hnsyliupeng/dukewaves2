function [ls] = calclevelset(x,y,z,numnod)
rad = 0.45;
for nod=1:numnod
    %     if (y(nod)<=0.5)
%             ls(nod) = z(nod) - 0.25;
    %     else
    %         ls(nod) = z(nod);
    %     end
%          rcirc(nod) = sqrt(x(nod)*x(nod)+y(nod)*y(nod));
         Rad(nod) = sqrt(x(nod)*x(nod) + y(nod)*y(nod) + z(nod)*z(nod));
%          if(rcirc(nod)<=rad)
%              ls(nod) = Rad(nod)-rad;
%          else
%              ls(nod) = z(nod);
%          end
    ls(nod) = Rad(nod)-rad;
% ls(nod) = z(nod)-0;
% ls(nod) = (1/(0.46))*z(nod)+(1/(23))*y(nod)-1;
end