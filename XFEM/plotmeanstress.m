% plotmeanstress.m
%
% Plots the mean stress for the polygrain tensile test. 
%

% Author: Matthias Mayr (09/2010)

loadingsteps = 101;
unloadingsteps = loadingsteps - 1;
comp = 4;
g = atan(0.001);

figure(2);
% clf;
hold on;
plot(meanstress(1:loadingsteps,1)*2*g/(loadingsteps-1),meanstress(1:loadingsteps,comp),...
  'r-','LineWidth',4);
xlabel('shear angle','FontSize',36);
ylabel('spatially averaged shear stress','FontSize',36);

loadstep = [];
for i=0:unloadingsteps
  loadstep = [loadstep loadingsteps-1-i];
end;

plot(loadstep*2*g/unloadingsteps,meanstress(loadingsteps:end,comp),'r:','LineWidth',4);