% Calculates the coefficient of variation for multiple iterations
% of optimized parameters.
%
% pars = [D; B;                                  % Convex combination parameters 1-2
%    A;                                          % Wall strain parameter 3 
%    K; K_b; K_pb; K_pr; K_s;                    % Gains 4-8
%    tau; tau_b; tau_p; tau_r; tau_s; tau_H;     % Time Constants 9-14
%    q_w; q_p; q_r; q_s;                         % Sigmoid Steepnesses 15-18
%    xi_w; xi_p; xi_r; xi_s;                     % Sigmoid Shifts 19-22
%    HI; H_p; H_r; H_s;                          % Heart Rate Parameters 23-26
%    HIa; tchar];

clear; close all;
pars_names = {'$D$', '$B$', '$A$', ...
    '$K$','$K_b$','$K_p$','$K_r$','$K_s$', ...
    '$\tau$','$\tau_b$','$\tau_p$','$\tau_r$','$\tau_s$','$\tau_H$',...
    '$q_w$','$q_p$','$q_r$','$q_s$', ...
    '$\xi_w$','$\xi_p$','$\xi_r$','$\xi_s$', ...
    '$H_I$','$H_[$','$H_r$','$H_s$', ...
    '$H_{Ia}$','tchar'};

N = 28;

all_pars = zeros(N,1);
for iter=1:9
   s = strcat('Ex1_opt',num2str(iter),'.mat');
   cd ../Optimized
   load(s);
   parsnom = exp(opt_info.pars);
   parsopt = exp(optpars);
   all_pars_opt(:,iter) = parsopt; 
   all_pars_nom(:,iter) = parsnom; 
   all_cost_grad(:,iter)= opt_info.histout(end,1:2)';
end
for i = 1:N
  mean_pars_opt(i) = mean(all_pars_opt(i,:));
  sd_pars_opt(i)   = std(all_pars_opt(i,:));
end;
CoV = sd_pars_opt./mean_pars_opt;

x = [1:length(parsopt)];
f = figure;
h = plot(x,CoV,'x','Markersize',10,'linewidth',3);
set(gca,'FontSize',18);
title('All pars');
xticks(1:length(pars));
set(gca,'XTickLabel',pars_names,'TickLabelInterpreter','latex','fontsize',18);
xtickangle(45);
grid on;
ylabel('Coeff Var')
xlim([1 length(pars)]);

disp('optpars')
all_pars_opt
disp('nompars')
all_pars_nom
disp('cost & grad')
all_cost_grad

cd ../../Core
