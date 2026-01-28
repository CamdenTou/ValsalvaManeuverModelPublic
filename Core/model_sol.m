function [HR,rout,J,Outputs,Init] = model_sol(pars,data)
% function model_sol
% Input: model parameters and data structure
% Output: computed heart rate, model residual, least squares cost, outputs and initial conditions
% Uses: model.m via ODE15s and initialconditions to set the initial
% conditions for the ODEs
% Description: solves the differential equations model

% Exponentiates model parameters 
pars = exp(pars); 
    
% Unpack data structure 
Tdata   = data.Tdata;
SPdata  = data.SPdata;
PPdata  = data.PPdata;
Pthdata = data.Pthdata; 
Hdata   = data.HRdata; 
dSPdt   = data.dSBPdt; 
dPthdt  = data.dPthdt; 
dPPdt   = data.dPPdt; 
i_ts    = data.i_ts;
i_t3    = data.i_t3;
    
ODE_TOL = data.ODE_TOL; 
    
% Convex combination of SPdata and PP data 
% Model parameter
D = pars(1); 
Pdata = D * SPdata + (1 - D) * PPdata;
dPdt  = D * dSPdt  + (1 - D) * dPPdt; 

Texnt = 1;
dte   = 0.1;
data.Texnt = Texnt;
TdataOff   = [0:0.1:Texnt];
TdataL     = [TdataOff Tdata(2:end)'+Texnt]';
PdataL     = [Pdata(1)*ones(size(TdataOff)) Pdata(2:end)']';
dPdtL      = [dPdt(1)*ones(size(TdataOff)) dPdt(2:end)']';
PthdataL   = [Pthdata(1)*ones(size(TdataOff)) Pthdata(2:end)']';
dPthdtL    = [dPthdt(1)*ones(size(TdataOff)) dPthdt(2:end)']';
HdataL     = [Hdata(1)*ones(size(TdataOff)) Hdata(2:end)']';

data.P_spline      = griddedInterpolant(TdataL,PdataL); 
data.dPdt_spline   = griddedInterpolant(TdataL,dPdtL); 
data.Pth_spline    = griddedInterpolant(TdataL,PthdataL); 
data.dPthdt_spline = griddedInterpolant(TdataL,dPthdtL); 

%Respiration (NEW)
Rdata = data.Rdata; 
RdataL = [Rdata(1)*ones(size(TdataOff)) Rdata(2:end)'];

R_d = griddedInterpolant(TdataL, RdataL);

Tshift = length(TdataL) - length(Tdata) + 1; 
R_d.Values(i_ts+Tshift:i_t3+Tshift) = 0;
data.R_spline = R_d;

% Get initial conditions 
Init = initialconditions(pars,data);

% Solve the model 
opts  = odeset('RelTol',ODE_TOL,'AbsTol',ODE_TOL); 
sol   = ode15s(@model,[TdataL(1) TdataL(end)], Init, opts, pars, data); 
sols  = deval(TdataL,sol); 

% % Solve the model (0 to tchar)
% opts  = odeset('RelTol',ODE_TOL,'AbsTol',ODE_TOL); 
% sol   = ode15s(@model,TdataL, Init, opts, pars, data); 
% sols  = deval(TdataL,sol); 

% Solution
Pc  = sols(1,:)';
Pa  = sols(2,:)';
ebc = sols(3,:)';
eba = sols(4,:)';
Tp  = sols(5,:)';
Tr  = sols(6,:)';
Ts  = sols(7,:)';
HR  = sols(8,:)';

Ms    = length(HdataL)-length(Hdata)+1;
i_ts  = i_ts + Ms;
i_ts2 = round(0.9*i_ts);
i_t3 = i_t3 + Ms; 

% Calculate residual and cost 
R_out1 = (HdataL(Ms    :i_ts) - HR(Ms    :i_ts))./ HdataL(Ms   :i_ts); %max(HdataL); % HdataL(Ms:round(i_ts/2));
R_out2 = (HdataL(i_ts+1:i_t3) - HR(i_ts+1:i_t3))./ HdataL(i_ts+1:i_t3); %max(HdataL); % HdataL(i_ts:i_t3);
R_out3 = (HdataL(i_t3+1:end)  - HR(i_t3+1:end)) ./ HdataL(i_t3+1:end); %max(HdataL); %HdataL(i_t3+1:end);
R_out4 = (max(HdataL) - max(HR)) ./ max(HdataL) ./sqrt(length(HdataL(Ms:end)));

rout = [2*R_out1 / sqrt(length(R_out1)); 2*R_out2 / sqrt(length(R_out2)); R_out3 / sqrt(length(R_out3)); R_out4];
J = rout'*rout;

Outputs = [TdataL,Pc,Pa,ebc,eba,Tp,Tr,Ts]; 

end % function model_sol.m % 