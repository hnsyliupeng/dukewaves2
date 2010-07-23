% m-file-version of 'database.xlsx'. Data is copied manually. 
% This file is used to generate the convergence plots
%
% following matrices are filled with the following data:
% column    description
%   1       number of elements
%   2       number of base DOFs
%   3       number of all DOFs
%   4       L2-norm of displacement error
%   5       energy-norm
%   6       L2-norm of traction error
%

% Author: Matthias Mayr (05/2010)

lagrange = [252   308	  336	0.1519    0.27419	0.84239;
            820   924   968	0.043575	0.28494	0.9175;
					 3240	 3444	 3528	0.013543	0.28938	0.95969;
           7260	 7564	 7688	0.0077796	0.29049	0.97304;
          12880	13284	13448	0.0055011	0.29097	0.9797;
          51360	52164	52488	0.0026101	0.2916	0.98976];

penalty = [ 252	  308	  336	0.1514	  0.27483	0.084292;
            820	  924	  968	0.042927	0.28494	0.07652;
           3240	 3444	 3528	0.01344	  0.28938	0.014087;
           7260	 7564	 7688	0.0077651	0.29049	0.0084514;
          12880	13284	13448	0.0054999	0.29097	0.0036425;
          51360	52164	52488	0.002614	0.2916	0.00074396];

nitsche = [ 252	  308	  336	0.1502    0.27421	0.061199;
            820	  924	  968	0.04312   0.28494	0.010854;
					 3240	 3444	 3528	0.013974	0.28938	0.0060262;
           7260	 7564	 7688	0.0077744	0.29049	0.0034504;
          12880	13284	13448	0.0055    0.29097	0.0013597;
          51360	52164	52488	0.0025857	0.2916	0.00055273];
 
slope1 = [1.0e+3 1.0e-1;  % slope "2" in double logarithmic diagramm
  1.0e+4 1.0e-3];  

slope2 = [1.0e+3 1.0e-1;  % slope "1" in double logarithmic diagramm
  1.0e+4 1.0e-2]; 
        
% generate a plot of the displacement error
figure(1);
% hold on;
% plot(lagrange(:,2),lagrange(:,4),'b');
% plot(penalty(:,2),penalty(:,4),'r');
% plot(nitsche(:,2),nitsche(:,4),'g');
loglog(lagrange(:,2),lagrange(:,4),'-vb',...
       penalty(:,2),penalty(:,4),'-or',...
       nitsche(:,2),nitsche(:,4),'-sg',...
       slope1(:,1),slope1(:,2),'-k');
% grid on;
% axis([1e2 1e5 1e-4 1e1]);

% title('convergence plot for displacement')
xlabel('Number of base degrees of freedom')
ylabel('error in displacement field, measured in L_2 -norm');
legend('Lagrange mulitpliers','penalty method','Nitsche´s method',...
  'Location','SouthWest');

        
% generate a plot of the traction error
figure(2);
% hold on;
% plot(lagrange(:,2),lagrange(:,6),'b');
% plot(penalty(:,2),penalty(:,6),'r');
% plot(nitsche(:,2),nitsche(:,6),'g');
loglog(lagrange(:,2),lagrange(:,6),'-vb',...
       penalty(:,2),penalty(:,6),'-or',...
       nitsche(:,2),nitsche(:,6),'-sg',...
       slope2(:,1),slope2(:,2),'-k');
% grid on;
axis([1e2 1e5 1e-4 1e1]);

% title('convergence plot for normal traction at the inteface')
xlabel('Number of base degrees of freedom')
ylabel('error in traction field, measured in L_2 -norm');
legend('Lagrange mulitpliers','penalty method','Nitsche´s method',...
  'Location','SouthWest');




