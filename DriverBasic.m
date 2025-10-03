function DriverBasic
% Input: none
% Output: generated figures and excel file with analyzed data
% Description: Runs software

clear all
close all
warning off

% Select screen and figure properties. Recommendations for PCs and Macs
% listed. Default is set to a Mac.
prompt         = {'Type (1) for PC and (2) for MAC', ...
                  'Fontsize (PC: 12, Mac: 16)',...
                  'Markersize (PC: 10, Mac: 14)',...
                  'Linewidth (PC: 2, Mac: 3)'};
dlgtitle       = 'Select Figure Parameters'; 
definput       = {'2','16','14','3'};
dims           = [1 65];
answer         = inputdlg(prompt,dlgtitle,dims,definput);

PCcomp         = str2num(answer{1});
fontsize       = str2num(answer{2});
markersize     = str2num(answer{3});
linwidth       = str2num(answer{4});
cd 
if isempty(answer) == 0
    plotmarkers.fs = fontsize;
    plotmarkers.ms = markersize;
    plotmarkers.lwt= linwidth;
else
    return
end 
indx = 'notempty';

% Default location for generated figures
figFolder = fullfile(pwd,'Figures');

% Set figure size
screenSize   = get(0,'screensize');
width   = floor(screenSize(3)/2);
height  = floor(screenSize(4)/1.5);
left    = floor((screenSize(3)-width)/4);
bottom  = floor((screenSize(4)-height)/4);
if PCcomp == 1
    figSize = [left bottom width height];
    posCB = [figSize(1)+figSize(3)*0.8,10,figSize(3)*0.25,figSize(4)*0.05];
else
    figSize = 0.75.*screenSize;
    figSize(1:2) = 0.12*screenSize(3:4);
    figSize(4) = figSize(4)*0.90;
    posCB = [25+figSize(3)*0.7,10,figSize(3)*0.25,figSize(4)*0.05];
end
  
plotmarkers.screenSize = screenSize;
plotmarkers.figSize    = figSize;
plotmarkers.figFolder  = figFolder;
plotmarkers.posCB      = posCB;

save('Plotmarkers.mat','plotmarkers')

% Patient selection
while isempty(indx) == 0

prompt   = {'Select patients'};
dlgtitle = 'Patient List';
d        = dir('Labchart/*.mat');
liststr  = {d.name};
prompt  = {'Select a file',...
           'Hold control/command to slect multiple',''};

[indx,~] = listdlg('PromptString',prompt,...
    'SelectionMode','multiple','ListString',liststr,...
    'ListSize',[200,300],'Name','Patient Selection');

if isempty(indx) == 1
    return
end

list = liststr(indx);

% Perform operations to analyze data
prompt_op = {'Select an operation',...
    'Run operations in numerical order',''};
call_flag_list = {'1.  Patient Information', ...
    '2.  Electrocardiogram (ECG)',...
    '3.  Heart Rate (HR)',...
    '4.  Respiration',...
    '5.  Blood Pressure (BP)',...
    '6.  VM Phases',...
    '7.  Clinical Ratios',...
    '8.  Model Prediction (Nominal)',...
    '9.  Sensitivity Analysis',...
    '10. Optimization',...
    '11. Plot Model Predictions',...
    '12. Summary'};

[indx_op,tf] = listdlg('PromptString',prompt_op,...
    'SelectionMode','multiple','ListString',call_flag_list,...
    'ListSize',[200,300],'Name','Operation Selection',...
    'SelectionMode','single');

if isempty(indx_op) == 1
    notdone = false; 
else 
    notdone = true;
    cd Core 
end

while notdone
    call_flag = call_flag_list{indx_op};
    operations(list,call_flag,plotmarkers); % code used to select operations

    [indx_op,tf] = listdlg('PromptString',prompt_op,...
    'SelectionMode','multiple','ListString',call_flag_list,...
    'ListSize',[200,300],'Name','Operation Selection',...
    'SelectionMode','single');
    if isempty(indx_op) == 1
       notdone = false; 
    else
       call_flag = call_flag_list{indx_op};
    end
    
    if sum(indx_op) == 0 
        str = pwd;
        if strcmp('ValsalvaAnalyzer',str(end-17:end-2))==1
            return
        elseif strcmp('Core',str(end-3:end)) == 1
            cd ..
            return
        end
     end
  end

end % DriverBasic %