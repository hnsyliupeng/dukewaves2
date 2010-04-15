function [flag,a,b] = shared_coord(el_1,el_2)

% shared_node takes two elements and checks to see if they share 2 (and
% only 2) nodes.  If they share 0 , 1, or 3, flag is returned as 0.  If
% they share 2, flag is returned as 1, and the nodal id's are returned.

global CONN X Y

flag = 0;
a = 0;
b = 0;
x = [];
y = [];
tol = 0.00001;

nodes1 = CONN(:,el_1);
nodes2 = CONN(:,el_2);
exes1 = X(nodes1);
exes2 = X(nodes2);
wyes1 = Y(nodes1);
wyes2 = Y(nodes2);

cnt = 0;

for i = 1:3
    for j = 1:3
        if (abs(exes1(i) - exes2(j)) < tol) &&...
                (abs(wyes1(i) - wyes2(j)) < tol);
            i; 
            j;
            cnt = cnt + 1;
            x = [x; exes1(i)];
            y = [y; wyes1(i)];
            break
            
        end
    end
end

if cnt == 2
    flag = 1;
    a = x;
    b = y;
end
    