function [sens,sens_norm,Rsens,Isens] = model_sens(pars,data,plotmarkers,patient)
%function model_sens
%Input: model parameters, data structure, and figure setttings
%Output: sensitivity matrix (sens) and normalized sensitivities.
%Sensitivities wrt the model residuals Rsens, and the order of
%sensitivities Isens
%Uses: senseq.m using finite differences to compute sensitivities
%Description: Using nominal parameter values (pars) and initial conditions for the ODEs
%senseq finds the non-weighted sensitivities

% Figure properties
fs = plotmarkers.fs;
ms = plotmarkers.ms;
lw = plotmarkers.lwt;
figSize = plotmarkers.figSize;

% Patient name
% patient = data.patient;

% Computing sensitivities wrt model residual
disp('Sensitivity analysis: It takes about a minute ...');
sens = senseq(pars,data);  % dr/dlog(theta)

% ranked classical sensitivities
[M,N] = size(sens);
for i = 1:N
  sens_norm(i)=norm(sens(:,i),2);  % dr/dtheta 
end

% Sorting sensitivities from most to least sensitive
[Rsens,Isens] = sort(sens_norm,'descend'); % Isens is sensitivity ranking order, Rsens are the sorted sensitivities.

% Plot ranked sensitivities
f = figure(1); clf
set(gcf,'units','points','position',figSize)
h=semilogy(Rsens(1:end-1)./max(Rsens),'x');
title(patient,'interpreter','none')
set(h,'linewidth',lw);
set(h,'Markersize',ms);
set(gca,'Fontsize',fs);
xticks(1:length(pars)-1);
for i = 1:length(pars)
  parsNames{i} = data.pars_names{Isens(i)};
end
set(gca,'XTickLabel',parsNames(1:end-1),'TickLabelInterpreter','latex','fontsize',fs+5);
xtickangle(45)
grid on;
ylabel('Sensitivites');
xlabel('Parameters')
xlim([0.5 length(pars)-0.5]);

figFolder = plotmarkers.figFolder;
fileName = strcat(patient,'_sensitivities.png');
fullFileName = fullfile(figFolder,'Model_fits',fileName);
exportgraphics(f,fullFileName);

bbutton = uicontrol('Parent',f,'Style','pushbutton',...
    'Position',[25,10,figSize(3)*0.25,figSize(4)*0.05],...
    'String','Save and exit','fontsize',fs,'Callback',@closeGenButton);

uiwait(f);

end % function model_sens.m %
