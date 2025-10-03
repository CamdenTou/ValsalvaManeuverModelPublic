function Init = initialconditions(pars,data)
% function initialconditions.m
% Input: model parameters and data structure
% Output: initial conditions for differential equations

% Model parameters
alpha  = pars(1); 
beta   = pars(2);
A      = pars(3);
    
K      = pars(4);
K_b    = pars(5);
K_p    = pars(6);
K_r    = pars(7);
K_s    = pars(8);
    
q_w    = pars(15);
q_p    = pars(16);
q_s    = pars(17);
    
xi_w   = pars(18);
xi_p   = pars(19);
xi_s   = pars(20);
    
HI     = pars(21);

respAmp = pars(27);
    
% Calculate initial conditions ensuring that the model starts in steady
% state
Pbar      = alpha * data.SPbar + (1 - alpha) * data.PPbar; 
dPdt_bar  = alpha * data.dSBPdt_bar + (1 - alpha) * data.dPPdt_bar; 
dPadt_bar = dPdt_bar - data.dPthdt_bar; 
    
Pc_m_0 = Pbar + K * dPdt_bar; 
Pa_m_0 = (Pbar - data.Pthbar) + K * dPadt_bar; 
    
y_c_0 = -q_w * (Pc_m_0 - xi_w); 
y_a_0 = -q_w * (Pa_m_0 - xi_w); 
    
e_wc_0 = 1 - sqrt((1 + exp(y_c_0))/(A + exp(y_c_0))); 
e_wa_0 = 1 - sqrt((1 + exp(y_a_0))/(A + exp(y_a_0))); 
    
e_bc_0 = K_b * e_wc_0; 
e_ba_0 = K_b * e_wa_0;  
    
n_0   = beta*(e_wc_0 - e_bc_0) + (1 - beta)*(e_wa_0 - e_ba_0);
    
G_p_0  = 1/(1 + exp(-q_p*(n_0 - xi_p)));
% G_r_0  = 1/(1 + exp( q_r*(data.Pthbar - xi_r)));
G_r_0 =  respAmp * data.Rbar;
G_s_0  = 1/(1 + exp( q_s*(n_0 - xi_s))); 
    
T_p_0  = K_p*G_p_0; 
T_r_0  = K_r*G_r_0;
T_s_0  = K_s*G_s_0; 
    

% Output vector with initial conditions
Init = [Pc_m_0; Pa_m_0; e_bc_0; e_ba_0; T_p_0; T_r_0; T_s_0; HI];

end % function initialconditions.m % 