function dxdt = model(t,x,pars,data)
% function model.m (called from ODE15.s) 
% Input: current time and state, model parameters and data structure
% Output: the rate of change of the model states at time t
% Description: Right hand side of ordinary differential equations

% Assign pressure values from data splines 
Pc_d     = data.P_spline; 
Pth_d    = data.Pth_spline; 
dPc_ddt  = data.dPdt_spline; 
dPth_ddt = data.dPthdt_spline; 
R_d = data.R_spline; % NEW
    
Pc  = Pc_d(t); 
Pa  = Pc_d(t) - Pth_d(t);
Pth = Pth_d(t); 

Rd = R_d(t);
    
dPcdt = dPc_ddt(t); 
dPadt = dPc_ddt(t) - dPth_ddt(t);
    
% Model parameters 
B   = pars(2); % Convex combination parameter 
A   = pars(3); % Wall strain parameter 
    
% Gains 
K_P    = pars(4);       % unitless
K_b    = pars(5);       % unitless
K_p    = pars(6);       % unitless
K_r    = pars(7);       % unitless
K_s    = pars(8);       % unitless
    
% Time scales
tau    = pars(9);       % s 
tau_b  = pars(10);      % s
tau_p  = pars(11);      % s 
tau_r  = pars(12);      % s 
tau_s  = pars(13);      % s
tau_H  = pars(14);      % s
    
% Sigmoid steepnesses
q_w    = pars(15);      % mmHg^{-1}
q_p    = pars(16);      % unitless
q_s    = pars(17);      % unitless
    
% Sigmoid shifts 
xi_w   = pars(18);      % mmHg
xi_p   = pars(19);      % unitless 
xi_s   = pars(20);      % unitless
    
% Heart rate parameters 
H_I    = pars(21);      % bpm
H_p    = pars(22);      % unitless
H_r    = pars(23);      % unitless
H_s    = pars(24);      % unitless 
    
% HR drop after phase 3 
H_Ia  = pars(25);        % bpm
tchar = pars(26)+data.Texnt; % Move since initial condition is shifted

% respiration amplitude change
respAmp = pars(27);

% Current states
Pc_m   = x(1);          % mmHg
Pa_m   = x(2);          % mmHg
eps_bc = x(3);          % unitless
eps_ba = x(4);          % unitless
T_p    = x(5);          % unitless
T_r    = x(6);          % unitless
T_s    = max(x(7),0);   % unitless 
H      = x(8);          % bpm
    
 %% Auxiliary equations 
    
 % Exponents in wall strain equation 
 y_c = -q_w * (Pc_m - xi_w);
 y_a = -q_w * (Pa_m - xi_w);
    
 % Wall strain
 eps_wc = 1 - sqrt((1 + exp(y_c)) / (A + exp(y_c)));
 eps_wa = 1 - sqrt((1 + exp(y_a)) / (A + exp(y_a)));
    
 % Convex combination of baroreceptor strain 
 n = B*(eps_wc - eps_bc) + (1 - B)*(eps_wa - eps_ba);
    
 % Sigmoidal function 
 G_p  = 1/(1 + exp(-q_p*(n - xi_p)));
 G_s  = 1/(1 + exp( q_s*(n - xi_s)));
 % G_r  = 1/(1 + exp( q_r*(Pth - xi_r)));
 G_r = respAmp * Rd; 
    
if t > tchar 
   H_I = H_Ia; 
end 
    
 % Heart rate (bpm) 
 H_tilde = H_I*(1 - H_p*T_p + H_r*T_r + H_s*T_s);
    
 %% ODEs     
 dPc_mdt = (-Pc_m + Pc - K_P * dPcdt)/tau;
 dPa_mdt = (-Pa_m + Pa - K_P * dPadt)/tau;
 debcdt  = (-eps_bc + K_b*eps_wc)/tau_b;
 debadt  = (-eps_ba + K_b*eps_wa)/tau_b;
 dTpdt   = (-T_p + K_p*G_p )/tau_p;
 dTrdt   = (-T_r + K_r*G_r)/tau_r;
 dTsdt   = (-T_s + K_s*G_s)/tau_s;
 dHdt    = (-H + H_tilde)/tau_H;
 
 % The vector of the right hand sides of the ODEs
 dxdt = [dPc_mdt; dPa_mdt; debcdt; debadt; dTpdt; dTrdt; dTsdt; dHdt;];

end % function model.m %