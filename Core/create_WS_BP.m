function data = create_WS_BP(patient,plotmarkers)
% function create_WS_BP
% Input: patient name and plot markers
% Output: augmented data, figures displaying systolic and diastolic
% pressure and workspace with systolic and diastolic pressures
% Uses: period_test needed to identify approximate length of cardiac cycle
% Description: generates spline with systolic and diastolic blood
% pressures. Saves figures, generates a workspace, and augments the data
% structure.

% global Traw Praw

%% Load temporary workspace and assign vectors
load(strcat('../WS/',patient,'_WS_temp.mat'));

Traw   = raw_data.Traw-raw_data.Traw(1);
Praw   = raw_data.Praw;
%% Figure properties
fs = plotmarkers.fs;
lw = plotmarkers.lwt;
if lw > 2
    lw2 = lw -1;
else
    lw2 = 1;
end
ms = plotmarkers.ms;
figSize = plotmarkers.figSize;

%% Correcting systolic and diastoloic blood pressure signals
[SPinds, DPinds] = period_test(Traw,Praw,patient,plotmarkers);
PrawOrg = Praw;
SPindsOrg = SPinds;
DPindsOrg = DPinds;

%% Plot results
f = figure(2); clf; hold on
set(gcf,'units','points','position',figSize)
f.Units = 'pixels';

p1 = plot(Traw,Praw,'b-','DisplayName','BP (mmHg)');
p2 = plot(Traw(DPinds),Praw(DPinds),'-o','Color',[0, 0.5, 0],'DisplayName','Diastolic','Linewidth',lw2);
p3 = plot(Traw(SPinds),Praw(SPinds),'-or','DisplayName','Systolic','Linewidth',lw2);
set(gca,'Fontsize',fs);
title(patient,'Scrolling (enter when done)', 'interpreter', 'none');
xlabel('Time (s)')
ylabel('BP (mmHg)');
legend([p1 p2 p3],{'BP','Diastolic BP','Systolic BP'})
xlim([0 max(Traw)])
grid on;

answer_SBP = questdlg('Do you want to correct systolic points?', ...
	'Systolic Blood Pressure','Yes','No ','Yes');

if strcmp(answer_SBP,'Yes') == 1 
    figure(2); hold on;
        pan on
        xlim([0 30]) % make window smaller so you can see the first 15 seconds
        ylim([min(Praw(DPinds))-5 max(Praw(SPinds))+5])
        title(patient,'Scrolling (enter when done)', 'interpreter', 'none');
        pause
        
    figure(2); hold on
        title(patient,'To finish press enter w/o clicking mouse', 'interpreter', 'none');
        xlabel('Click on SBP points to fix (enter when done)');
        [s,ys] = ginput;
        
    while isempty(s) == 0 % make sure there was an input
        i_s = zeros(1,length(s));
        for i = 1:length(s)
            i_st = find(Traw <= s(i),1,'last'); %find time closest to input
            if isempty(i_st) == 1
                i_st = 1;
            end
            if i_st+50 >= max(length(Praw))
                i_sp   = find(Praw(i_st-50:end) == max(Praw(i_st-50:end)),1); %if it misses the last peak
                i_sp   = i_sp + i_st-50;
            elseif i_st-50 <= 1
                i_sp   = find(Praw(1:i_st+50) == max(Praw(1:i_st+50)),1); %if it misses the first peak
                i_sp   = i_sp + i_st;
            else
                i_sp   = find(Praw(i_st-50:i_st+50) == max(Praw(i_st-50:i_st+50)),1); %find peak within one second of click
                i_sp   = i_sp + i_st-50;
            end
                
            try
                i_s(i) = i_sp;
            catch
                i_s(i);
                i_sp;
            end
        end
            
        [~,closestIndFirst] = min(abs(SPinds-i_s(1)));
        [~,closestIndLast]  = min(abs(SPinds-i_s(end)));
            
        if Traw(i_s(end)) > Traw(max(SPinds)) %last point is greater than current last sp ind
            if length(s) == 2
                SPinds = [SPinds(1:closestIndFirst) i_s(end)]; %might be able to get rid of this case logically
            else
                SPinds = [SPinds(1:closestIndFirst) i_s(2:end)];
            end
        elseif Traw(i_s(1)) < Traw(SPinds(1)) % first point is less than current sp ind
            if length(s) == 2
                SPinds = [i_s(1) SPinds(closestIndLast:end)];
            else
                SPinds = [SPinds(closestIndFirst) i_s(2:end-1) SPinds(closestIndLast:end)]; %might be able to eliminate this case logically
            end
        else
            if length(s) == 2
                SPinds = [SPinds(1:closestIndFirst) SPinds(closestIndLast:end)];
            else
                SPinds = [SPinds(1:closestIndFirst) i_s(2:end-1) SPinds(closestIndLast:end)];
            end
        end
        
        Praw(i_s) = ys;  
        plot(Traw(i_s),Praw(i_s),'Or--','linewidth',lw2,'DisplayName','Data Correction');
        set(gca,'Fontsize',fs);
        xlabel('Time (s)');
        ylabel('BP (mmHg)');
        legend([p1 p2 p3],{'BP','Diastolic BP','Systolic BP'})
        
        figure(2); hold on
        pan on
        title(patient,'Scrolling (enter when done)', 'interpreter', 'none');
        pause
            
        figure(2); hold on
        title(patient,'To finish press enter w/o clicking mouse', 'interpreter', 'none');
        [s,ys] = ginput;
    end
end
    
answer_DBP = questdlg('Do you want to correct diastolic points?', ...
                      'Diastolic Blood Pressure','Yes','No ','Yes');
    
if strcmp(answer_DBP,'Yes') == 1
    figure(2); hold on
    pan on
    xlim([0 15]) % make window smaller so you can see the first 15 seconds
    ylim([min(Praw(DPinds))-5 max(Praw(SPinds))+5])
    xlabel('Time (s)');
    pause
    
    figure(2); hold on
    title(patient,'To finish press enter w/o clicking mouse', 'interpreter', 'none');
    xlabel('Click on DBP points to fix (enter when done)');
    [d,yd] = ginput;
    
    while isempty(d) == 0 % make sure there was an input
        i_d = zeros(1,length(d));
        for i = 1:length(d)
            i_dt   = find(Traw <= d(i),1,'last'); %find time closest to input
            if isempty(i_dt) == 1
                i_dt =1;
            end
            if i_dt+50 >= max(length(Praw))
                i_dp   = find(Praw(i_dt-50:end) == min(Praw(i_dt-50:end)),1); %if it misses the last peak
                i_dp   = i_dp + i_dt-50;
            elseif i_dt-100 <= 1
                i_dp   = find(Praw(1:i_dt+50) == min(Praw(1:i_dt+50)),1); %if it misses the first peak
                i_dp   = i_dp + i_dt;
            else
                i_dp   = find(Praw(i_dt-50:i_dt+50) == min(Praw(i_dt-50:i_dt+50)),1); %find min within one second of click
                i_dp   = i_dp + i_dt-50;
            end
            i_d(i) = i_dp;
        end
        try
            i_d(i) = i_dp;
        catch
            i_d(i);
            i_dp;
        end
        
        if i_d(end) >= length(Traw)
            i_d(end) = length(Traw);
        end
        [~,DclosestIndFirst] = min(abs(DPinds-i_d(1)));
        [~,DclosestIndLast]  = min(abs(DPinds-i_d(end)));
        
        
        if Traw(i_d(end)) > Traw(max(DPinds)) %last point is greater than current last ind
            if length(d) == 2
                DPinds = [DPinds(1:DclosestIndFirst) i_d(end)]; %might be able to get rid of this case logically
            else
               DPinds = [DPinds(1:DclosestIndFirst) i_d(2:end)];
            end
        elseif Traw(i_d(1)) < Traw(DPinds(1)) % first point is less than current sp ind
            if length(d) == 2
                DPinds = [i_d(1) DPinds(DclosestIndLast:end)];
            else
                DPinds = [DPinds(DclosestIndFirst) i_d(2:end-1) DPinds(DclosestIndLast:end)]; %might be able to eliminate this case logically
            end
        else
            if length(d) == 2
                DPinds = [DPinds(1:DclosestIndFirst) DPinds(DclosestIndLast:end)];
            else
                DPinds = [DPinds(1:DclosestIndFirst) i_d(2:end-1) DPinds(DclosestIndLast:end)];
            end
        end
        
        Praw(i_d) = yd;  
        plot(Traw(i_d),Praw(i_d),'Og--','linewidth',lw2,'DisplayName','Data Correction');
        xlabel('time (sec)');
        ylabel('BP (mmHg)');
        legend([p1 p2 p3],{'BP','Diastolic BP','Systolic BP'})
        set(gca,'Fontsize',fs);
        
        figure(2); hold on
        pan on
        xlabel('Scrolling in progress (enter when done)');
        pause
        
        figure(2); hold on
        xlabel('Click on DBP points to fix (enter when done)');
        title(patient,'To finish press enter w/o clicking mouse', 'interpreter', 'none');
        [d,yd] = ginput;
    end
end

bbutton = uicontrol('Parent',f,'Style','pushbutton',...
    'Position',[25,10,figSize(3)*0.25,figSize(4)*0.05],...
    'String','Save and exit','fontsize',fs,'Callback',@closeGenButton);

uiwait(f)

%% Save temporary workspace
SPinds = unique(SPinds);
DPinds = unique(DPinds);
SPraw  = interp1(Traw(SPinds),Praw(SPinds),Traw,'pchip'); % extrapolating pressures from click on time vector from data 1000Hz
DPraw  = interp1(Traw(DPinds),Praw(DPinds),Traw,'pchip'); % extrapolating pressures from click on time vector from data 1000Hz
PPraw  = SPraw - DPraw;

data.PrawOrg = PrawOrg;
data.SPindsOrg = SPindsOrg;
data.SPinds = SPinds;
data.DPindsOrg = DPindsOrg;
data.DPinds = DPinds;
data.SPraw = SPraw;
data.DPraw = DPraw;
data.PPraw = PPraw;

s = strcat('../WS/',patient,'_WS_temp.mat');
save(s,'data','-append');

end % function create_WS_BP.m %
