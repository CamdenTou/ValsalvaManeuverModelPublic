function [] = operations(list,call_flag,plotmarkers)
% function operations.m
% Input: list of datasets, flag with operation, figure markers
% Output: None
% Description: Runs operations required to analyze data
scp = strcmp(call_flag,'12. Summary');

if scp ~= 1
    for i = 1:length(list)
      close all
      patient  = list{i}(1:end-4);
     
      switch call_flag
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case '1.  Patient Information'
            create_WS_info(patient,plotmarkers);
            close all;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case '2.  Electrocardiogram (ECG)'
            create_WS_ECG(patient,plotmarkers);
            close all;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case '3.  Heart Rate (HR)'
            create_WS_HR(patient,plotmarkers);  
            close all;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        case '4.  Respiration'
            create_WS_Resp(patient,plotmarkers);  
            close all;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        case '5.  Blood Pressure (BP)'
            create_WS_BP(patient,plotmarkers);
            close all;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        case '6.  VM Phases'
            create_WS_indices(patient,plotmarkers);
            close all;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        case '7.  Clinical Ratios'
            clinicalratios(patient,plotmarkers);
            close all;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case '8.  Model Prediction (Nominal)'
            load(strcat('../WS/',patient,'_WS.mat'), 'data');
            % Get nominal parameter values
            [data,pars] = parameters(data); 
            
            % Run forward model with nominal parameter values 
            [HR,~,~,Outputs] = model_sol(pars,data);  
            
            % Plot forward model with nominal parameter values
            plot_model(patient,Outputs,HR,data,pars,0,plotmarkers);
            close all;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case '9.  Sensitivity Analysis'
            
            load(strcat('../WS/',patient,'_WS.mat'), 'data');
           
            % Get nominal parameter values
            [data,pars] = parameters(data);
            [sens,sens_norm,Rsens,Isens] = model_sens(pars,data,plotmarkers,patient); 
            
            save(strcat('../Sensitivities/',patient,'_sens.mat'),'sens','sens_norm','Rsens','Isens');
            close all;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        case '10. Optimization'
            
            load(strcat('../WS/',patient,'_WS.mat'), 'data')
        
            % Get nominal parameter values
            [data, pars] = parameters(data);
            
            [optpars,opt_info] = model_opt(pars,data);
            
            % Run forward model with optimized parameters
            [HR_opt,~,~,Outputs_opt] = model_sol(optpars,data);
            
            save(strcat('../Optimized/',patient,'_opt.mat'), ...
                'HR_opt','Outputs_opt','optpars','pars','opt_info');
            close all;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        case '11. Plot Model Predictions'
            
            answer = questdlg('Select model prediction to view', ...
                'View Model Prediction', ...
                'Nominal','Optimized','Optimized');
            
            switch answer
                case 'Nominal'
                    load(strcat('../WS/',patient,'_WS.mat'), 'data');
                    
                    % Get nominal parameter values
                    [data,pars] = parameters(data);
                    
                    % Run forward model with nominal parameter values
                    [HR,~,~,Outputs] = model_sol(pars,data);
                    
                    % Plot forward model with nominal parameter values
                    plot_model(patient,Outputs,HR,data,pars,0,plotmarkers);
                    
                    % add save and exit button to plot
                case 'Optimized'
                    % Load workspace
                    load(strcat('../WS/',patient,'_WS.mat'), 'data');

                    % Get nominal parameter values
                    [data] = parameters(data);
                    
                    % Load optimized paramters
                    load(strcat('../Optimized/',patient,'_opt.mat'))
                    
                    [HR_opt,~,~,Outputs_opt] = model_sol(optpars,data);
                    
                    % Plot forward model with optimized parameters
                    plot_model(patient,Outputs_opt,HR_opt,data,optpars,1,plotmarkers); %plot results 
            end
            close all;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      end
    end
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   case '12. Summary'
        T = DataSummary(list);
        disp(T)
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

end % function Operations.m %

