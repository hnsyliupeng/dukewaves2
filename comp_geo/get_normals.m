function [norm,tang] = get_normals(xcoords,ycoords)

tangx = xcoords(2) - xcoords(1);
tangy = ycoords(2) - ycoords(1);

tang = [tangx tangy]'./(sqrt(tangx^2 + tangy^2));

norm = [-tang(2) tang(1)]';

