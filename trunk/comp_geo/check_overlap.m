function flag = check_overlap(x1,y1,x2,y2,tol);

%checks if distance between two points is within a tolerance

dist = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );

if (dist < tol)
    flag = 1;
else
    flag = 0;
end