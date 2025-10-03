function [data,pars] = parameters(data)
% function parameters.m
% Input: data structure
% Outout: list of model parameters and augmented data structure
% Description: sets nominal model parameters

% Integration tolerance
data.ODE_TOL  = 1e-6;
data.DIFF_INC = sqrt(1e-6);

% Initial mean values
SPbar  = data.SPbar;
PPbar  = data.PPbar;
Pthbar = data.Pthbar;
HminR  = data.HminR; 
HmaxR  = data.HmaxR; 
Hbar   = data.Hbar;
HRdata = data.HRdata;
Tdata  = data.Tdata;
i_HR4e_min = data.i_HR4e_min;
i_HR4e_Max = data.i_HR4e_Max; 

% Model parameters
D  = 0.75; 
B  = 0.4; 
A  = 5;         

% Gains
K    = 1; 
K_b  = .1;
K_p  = 1;
K_r  = 1; 
K_s  = 1; 

% Time scales
tau    = 1;  
tau_b  = 1;              
tau_p  = 0.5;          
tau_r  = 4;
tau_s  = 10;          
tau_H  = .1;

% Sigmoid shifts
q_w  = .04;         
q_p  = 10;          
q_s  = 10;   

% Patient specific parameters
xi_w  = D * SPbar + (1 - D) * PPbar;           

% Intrinsic HR
HI     = HminR; 
iH_I   = i_HR4e_Max;  
tchar  = Tdata(floor(iH_I));

% Maximal HR
HM  = 1.5*max(HRdata);    
H_s = (1/K_s)*(HM/HI - 1); 

% At end of expiration and inspiration
% Gr_ss = 1/(1 + exp(q_r*(Pthbar - xi_r)));

respAmp = 2*10^2;

Gr_ss = respAmp * data.Rbar;

Tr_ss  = K_r*Gr_ss; 
H_r    = (HmaxR - HminR)/HI/Tr_ss ;

% Calculate sigmoid shifts
Pc_ss  = xi_w; 
Pa_ss  = xi_w - Pthbar; 

ewc_ss = 1 - sqrt((1 + exp(-q_w*(Pc_ss - xi_w)))/(A + exp(-q_w*(Pc_ss - xi_w)))); 
ewa_ss = 1 - sqrt((1 + exp(-q_w*(Pa_ss - xi_w)))/(A + exp(-q_w*(Pa_ss - xi_w)))); 

ebc_ss = K_b*ewc_ss; 
eba_ss = K_b*ewa_ss;  

n_ss  = B*(ewc_ss - ebc_ss) + (1 - B)*(ewa_ss - eba_ss);

Tp_ss  = .8;
Ts_ss  = .2; 

% Steady-state sigmoid shifts 
xi_p  = n_ss + log(K_p/Tp_ss - 1)/q_p;  
xi_s  = n_ss - log(K_s/Ts_ss - 1)/q_s;   

H_p   = (H_r*Tr_ss + H_s*Ts_ss)/Tp_ss;
HIa   = min(HRdata(i_HR4e_min:end));

% delay value
Ds = 3;

% Parameter vector (output)
pars = [D; B;                                   % Convex combination parameters 1-2
    A;                                          % Wall strain parameter 3 
    K; K_b; K_p; K_r; K_s;                      % Gains 4-8
    tau; tau_b; tau_p; tau_r; tau_s; tau_H;     % Time Constants 9-14
    q_w; q_p; q_s;                              % Sigmoid Steepnesses 15-17
    xi_w; xi_p; xi_s;                           % Sigmoid Shifts 18-20
    HI; H_p; H_r; H_s;                          % Heart Rate Parameters 21-24
    HIa; tchar;                                 % Mean HR after VM and index for time at which it changes 25-26
    respAmp; Ds];                                % Respiration data amplification, delay 27-28
               
pars_names = {'$\alpha$', '$\beta$', '$A$', ...
    '$K$','$K_b$','$K_p$','$K_r$','$K_s$', ...
    '$\tau$','$\tau_b$','$\tau_p$','$\tau_r$','$\tau_s$','$\tau_H$',...
    '$q_w$','$q_p$','$q_s$', ...
    '$\xi_w$','$\xi_p$','$\xi_s$', ...
    '$H_I$','$H_p$','$H_r$','$H_s$', ...
    '$H_{Ia}$','$tchar$','$respAmp$','$Ds$'};

% Parameter bounds

% Upper and lower parameter bounds
lb      = pars/10; 
ub      = pars*10;
ub(3)   = pars(3) * 1.5;
lb(3)   = pars(3) / 2;
lb([5,6,7,8,11,12,14,23,24])   = 0.01;
ub([5,6,7,8,14])   = 10;
lb(10)  = 0.1;
ub(10)  = 1.5;
ub(11)  = 17.9;
ub(12)  = 31.2;
%lb([13,15,16,17,18])  = pars([13,15,16,17,18]) - 0.5*pars([13,15,16,17,18]);
lb([13,15,16,17])  = pars([13,15,16,17])/2;
ub([13,15,16,17])  = pars([13,15,16,17]) + 1.5*pars([13,15,16,17]); 
lb(18)  = 83;
ub(18)  = 163;
lb(19)  = 0.53;
ub(19)  = 0.55;
lb(20)  = 0.04;
ub(20)  = 0.06;
lb(21)  = 60;
ub(21)  = 114;
lb(22)  = 0.1;
ub(22)  = 0.9;
ub([23,24]) = 1.1;
ub(28)  = 15;
lb(28)  = 0.5;


% B - Convex combination bounds
lb([1,2])  = .01;                       
ub([1,2])  = 1;
%lb([11])   = 0.0001;

% Log scaled outputs
pars = log(pars);
lb   = log(lb);
ub   = log(ub); 
data.lb = lb;
data.ub = ub;
data.pars_names = pars_names; 

end  % function parameters.m %