function [] = create_WS_HR(patient,plotmarkers)
% function create_WS_HR
% Input: patient name and figure settings
% Output: saved temporary workspace (saved in WS)
% Uses: correct_HR to correct heart rate values
% Description: loads temporary workspace with corrected ECG data calculates
% heart rate and saves temporary workspace with Heart rate output

% nested functions bc using user info (blue color is not global just bc of
% nested)

global HR_new % has to be global because of nested functions

%% Figure properties
fs = plotmarkers.fs;
lw = plotmarkers.lwt;
figSize = plotmarkers.figSize;

if lw > 2
    lw2 = lw-1;
else
    lw2 =1;
end

%% Load temporary workspace and assign vectors
filename = strcat(patient,'_WS_temp.mat');
load(strcat('../WS/',filename),'data','raw_data'); 

Traw   = raw_data.Traw0;    % time starts at zero
ECGraw = raw_data.ECGraw;

RRint   = data.RRint;
T_RRint = data.T_RRint;

HR      = 60./RRint;

%% ECG to HR Calculations
f = figure(1); clf; hold on;
set(gcf,'units','points','position',figSize)
f.Units = 'pixels';
sgtitle(patient,'Fontsize',fs+2,'FontWeight','bold','Interpreter','none')
sub1 = subplot(2,1,1); hold on;
    plot(T_RRint,HR,'bo-','linewidth',lw2);
    set(gca,'fontsize',fs)
    ylabel('HR (bpm)')
    xlim([Traw(1),Traw(end)])
    ylim([min(HR)-5 max(HR)+5]);

sub2 = subplot(2,1,2); hold on;
    plot(Traw,ECGraw,'b-')
    set(gca,'fontsize',fs)
    xlabel('Time (s)')
    ylabel('ECG (mV)')
    xlim([Traw(1),Traw(end)])
linkaxes([sub1, sub2],'x')

figFolder = plotmarkers.figFolder;
fileName = strcat(patient,'_HeartRateECG.png');
fullFileName = fullfile(figFolder,'Data',fileName);
saveas(f,fullFileName,'png');

bbutton = uicontrol('Parent',f,'Style','pushbutton',...
    'Position',[25,10,figSize(3)*0.25,figSize(4)*0.05],...
    'String','Save and exit','fontsize',fs,'Callback',@closeButton);

cbutton = uicontrol('Parent',f,'Style','pushbutton',...
    'Position',plotmarkers.posCB,...
    'String','Correct Heart Rate','fontsize',fs,'Callback',@correctButton);

uiwait(f)


function [] = closeButton(a,~)
    f = a.Parent;
    HR_new = HR;
    close(f)
end % function close botton %

function [] = correctButton(a,~) 
    
    [HR_new] = correct_HR(HR,Traw,T_RRint,patient,plotmarkers);
    
    f = figure(1); clf; hold on;
    set(gcf,'units','points','position',figSize)
    f.Units = 'pixels';
    sgtitle(patient,'Fontsize',fs+2,'FontWeight','bold','Interpreter','none')
    sub1 = subplot(2,1,1); hold on;
        plot(T_RRint,HR_new,'bo-','linewidth',lw2);
        set(gca,'fontsize',fs)
        ylabel('HR (bpm)')
        xlim([Traw(1),Traw(end)])
        ylim([min(HR_new)-5 max(HR_new)+5]);
        
    sub2 = subplot(2,1,2); hold on;
        plot(Traw,ECGraw,'b-')
        set(gca,'fontsize',fs)
        xlabel('Time (s)')
        ylabel('ECG (mV)')
        xlim([Traw(1),Traw(end)])
        %ylim([-1.25e-3 1.25e-3]);
        
    linkaxes([sub1, sub2],'x')

    figFolder = plotmarkers.figFolder;
    fileName = strcat(patient,'_HeartRateECG');
    fullFileName = fullfile(figFolder,'Data',fileName);
    saveas(f,fullFileName,'png');

    bbutton = uicontrol('Parent',f,'Style','pushbutton',...
    'Position',[25,10,figSize(3)*0.25,figSize(4)*0.05],...
    'String','Save and exit','fontsize',fs,'Callback',@closeGenButton);
   
end % function correct botton %

% Interpolate over step function and evaluate at Tdata
HRraw = interp1(T_RRint,HR_new,Traw,'pchip');

% Use HR to find RR intervals
RRraw = 60./HRraw;


%% Save temporary workspace
% data.HR     = HR;
% data.HR_new = HR_new;   
raw_data.HRraw = HRraw;   % user corrected HR data
data.RRraw     = RRraw;

plotmarkers.fs = fs;
plotmarkers.figSize = figSize;
plotmarkers.sub1 = sub1;
plotmarkers.sub2 = sub2;
plotmarkers.lw2 = lw2;

% Save temporary workspace
s = strcat('../WS/',patient,'_WS_temp.mat');
save(s,'data','raw_data','plotmarkers','-append');
end % function create_WS_HR %

