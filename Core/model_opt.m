function [optpars,opt_info] = model_opt(pars,data)
%function model_opt
%Input: model parameters and data structure
%Output: optimized model parameters
%Uses: newlsq_v2.m gradient based levenberg Marquart to estimate parameters
%and model.m the right hand side of the ODEs
%Description: Estimates identifiable model parameters

% ODE solver tolerance
ODE_TOL = data.ODE_TOL; 
    
% lower and upper parameter bounds
lb = data.lb; 
ub = data.ub; 
    
% Global parameters
ALLPARS  = pars;
DIFF_INC = sqrt(ODE_TOL);

%pars = [D; B;                                   % Convex combination parameters 1-2
%    A;                                          % Wall strain parameter 3 
%    K_P; K_b; K_pb; K_pr; K_s;                  % Gains 4-8
%    tau; tau_b; taup_p tau_r; tau_s; tau_H;     % Time Constants 9-14
%    q_w; q_p; q_r; q_s;                         % Sigmoid Steepnesses 15-18
%    xi_w; xi_p; xi_r; xi_s;                     % Sigmoid Shifts 19-22
%    HI; H_p; H_r; H_s;                          % Heart Rate Parameters 23-26
%    HIa; tchar];                                % HR mean after VM 27          

% Parameters estimated
INDMAP = [2 11 12 13 21 22 23 24 25 ]; 
% Randomized initial parameters
%rng('Shuffle');
%rt = -1+2*rand(length(INDMAP),8);
%pars(INDMAP) = pars(INDMAP).*(1+0.1*rt(:,6));
%pars = log(pars);

% logscale parameters to be estimated
optx = pars(INDMAP);

% Global parameters
gpars.INDMAP   = INDMAP;
gpars.ALLPARS  = ALLPARS;
gpars.ODE_TOL  = ODE_TOL;
gpars.DIFF_INC = DIFF_INC;
data.gpars = gpars;

% Set optimization hyperparameters and run optimization 
optub = ub(INDMAP);
optlb = lb(INDMAP);
    
maxiter = 30; 
mode    = 2; 
nu0     = 2.d-1; 
    
% Run Levenberg-Marquardt optimization 
disp('Optimization: It takes a few minutes...')
[xopt, histout, costdata, jachist, xhist, rout, sc] = ...
         newlsq_v2(optx,'opt_wrap',1.d-3,maxiter,...
         mode,nu0,optub,optlb,data); 

% Outputs (estimated parameters)
optpars = pars;
optpars(INDMAP) = xopt;

opt_info.histout  = histout; 
opt_info.costdata = costdata; 
opt_info.jachist  = jachist; 
opt_info.xhist    = xhist; 
opt_info.rout     = rout; 
opt_info.sc       = sc; 
opt_info.ub       = ub(INDMAP);
opt_info.lb       = lb(INDMAP);
opt_info.INDMAP   = INDMAP;
opt_info.pars     = pars;
opt_info.optpars  = optpars;
    
end % function model_opt.m %
