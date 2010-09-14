% plotmeanstress.m
%
% Plots the mean stress for the polygrain tensile test. 
%

% Author: Matthias Mayr (09/2010)

loadingsteps = 21;
unloadingsteps = 20;
comp = 4;

figure(3);
clf;
hold on;
plot(meanstress(1:loadingsteps,1),meanstress(1:loadingsteps,comp),'*b-');

% loadstep = [];
% for i=0:unloadingsteps
%   loadstep = [loadstep loadingsteps-1-i];
% end;
% 
% plot(loadstep,meanstress(loadingsteps:end,comp),'*r-');