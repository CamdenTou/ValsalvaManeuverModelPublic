function [phases] = set_indices(set_ind, patient, plotmarkers)
% function set_indices.m
% Input: time, blood pressure, heart rate, valsalva information, patient
% name, and figure properties
% Output: Valsalva maneuver phases
% Uses:
% Description: Determines the valsalva maneuver phases

%close all
fs  = plotmarkers.fs;
lw  = plotmarkers.lwt;
if lw > 2
    lw2 = lw-1;
else
    lw2 = 1;
end
figSize = plotmarkers.figSize;

% Load data in from structure
Traw        = set_ind.Traw;
Praw        = set_ind.Praw;
Rdata       = set_ind.Rdata;
Tdata       = set_ind.Tdata;
SPdata      = set_ind.SPdata;
DPdata      = set_ind.DPdata;
PPdata      = set_ind.PPdata;
HRdata      = set_ind.HRdata;
RRdata      = set_ind.RRdata;
val_start   = set_ind.val_start;
val_end     = set_ind.val_end;
PthAvail    = set_ind.PthAvail;

 if PthAvail ~= 0
     Pthdata = set_ind.Pthdata;
 end

% Time Indices - breath hold start (ts) & breath hold end (te) !!!
% Calculate mean distance between timepoints
dt = mean(diff(Tdata));

% Find the indices of the time data points that are closest to the VM start
% and end times
i_ts  = find(Tdata>val_start,1);
i_t2l = find(Tdata>val_end,1);
i_t1_ini = i_ts + (i_t2l-i_ts)/2;

% Set rest before (r1) and after (r2)
i_r1 = 1;
i_r2 = length(Tdata);

% Find steady-state baseline values up to VM start
baselineOffsetSeconds = 2;
i_baselineOffset = baselineOffsetSeconds * round(1/dt);
SPbar = trapz(SPdata(1:i_ts - i_baselineOffset - 1))/(i_ts - i_baselineOffset - 1); % sp
DPbar = trapz(DPdata(1:i_ts - i_baselineOffset - 1))/(i_ts - i_baselineOffset - 1); % dp

% Find time for end of phase I
% Max blood pressure before start of late phase II
[~,k_SPmax1] = max(SPdata(i_ts:i_t1_ini));
i_t1 = i_ts + k_SPmax1 - 1;

SPsmooth = sgolayfilt(SPdata,3,7);
SP1st = gradient(SPsmooth,dt);

t1Offset = round(2/dt);
i_t1linear = findchangepts(SP1st(i_t1:i_t1+t1Offset),'Statistic','mean') + i_t1;

% Find time for end early phase II
% Select the point that which the SBP reaches a minimum during phase II
[~,k_t2e] = min(SPdata(i_t1:min(i_t2l-round(1/dt)))); % subtract 1 second to avoid drop off from breath hold
if isempty(k_t2e)
    k_t2e = round((i_ts-i_t1)/2);
end
i_t2emin = i_t1 + k_t2e - 1; % minimum phase 2 early
if i_t1linear >= i_t2emin  % CAMDEN ISSUE
    i_t1linear = i_t1;
end
if i_t2emin > i_t1linear
    t2echangePoint = findchangepts(SPdata(i_t1linear:i_t2emin),'Statistic','linear') + i_t1linear;
    i_t2elinear    = t2echangePoint;
else
    i_t2elinear    = i_t2emin;  % CAMNDEN - ISSUE
end
normalT2eThreshold = 1.10; % Threshold to check "normal" phase 2
% if normal use minimum point for phase II early
if SPdata(i_t2l-1) > normalT2eThreshold*SPdata(i_t2emin)
    i_t2e = i_t2emin; 
else
    i_t2e = i_t2elinear;
end

% Find end of phase III/start of phase II (min 5 seconds after end of BH
[~,k_t3] = min(SPdata(i_t2l:min(i_t2l+round(5/dt),length(SPdata))));
if isempty(k_t3)
    i_t3 = i_t2e;
else
    i_t3   = i_t2l + k_t3 - 1;
end

i_t3linear = findchangepts(SP1st(i_t3:i_t3+round(2/dt)),'Statistic','mean') + i_t3;

if SPdata(i_t3 + round(1/dt)) > 1.1*SPdata(i_t3)
    i_t3linear = i_t3;
end
 
% Check if the pressure at the end of phase III is greater than the mean
N = length(SPdata);
if SPdata(i_t3) > SPbar
    k_PRT = find(SPdata(i_t3+1:i_t3+round(5/dt)) > SPdata(i_t2l),1,'first');
else
    k_PRT = find(SPdata(i_t3+1:end) > SPbar,1,'first');
end
if isempty(k_PRT) == 1
    i_PRT = N;
else
    i_PRT = k_PRT + i_t3;
end

% Find time for end of phase IV
% Select the point at which the SBP returns to baseline after the overshoot
returnBaseTol = 1.5; % return to baseline tolerance 
% k_t4  = find((SPdata(i_PRT+1:end) <= SPbar + returnBaseTol) & (DPdata(i_PRT+1:end) <= DPbar + returnBaseTol), 1, 'first');
k_t4SP = find(SPdata(i_PRT+1:end) <= SPbar + returnBaseTol, 1, 'first');
k_t4DP = find(DPdata(i_PRT+1:end) <= DPbar + returnBaseTol, 1, 'first');

if ~isempty(k_t4SP) && ~isempty(k_t4DP)
    k_t4 = max(k_t4SP,k_t4DP);
    i_t4  = i_PRT+1 + k_t4;
    T4    = i_t4*dt;
end
MinP4 = 10;
MaxP4 = 30/dt;
if isempty(k_t4SP) || isempty(k_t4DP)
    i_t4 = length(SPdata);
    k_t4 = i_t4-i_PRT;
else
    SPdata_i_t4_mean = mean(SPdata(i_t4-1:i_t4)) - mean(SPdata(i_t4:i_t4+1));
    while (SPdata_i_t4_mean > 0 && i_t4 < length(SPdata)-1)
        i_t4 = i_t4 + 1;
        SPdata_i_t4_mean = mean(SPdata(i_t4-1:i_t4)) - mean(SPdata(i_t4:i_t4+1));
        % is not exceeding 0 bc there aren't enough indicies to exceed 0
    end
end
T4D = (i_t4-i_t3)*dt;

if T4D < MinP4
    notdone = 0;
    i = round(MinP4/dt)+i_t3;
    while notdone == 0
        if mean(SPdata(i-1:i))-mean(SPdata(i:i+1)) > 0
            i = i+1;
            if i == N
                notdone = 1;
                T4 = i*dt;
            end
        else
            notdone = 1;
            T4 = i*dt;
            i_t4 = i;
        end
    end
end

t4change = findchangepts(SP1st(i_t3+round(1/dt):i_t3+round(5.5/dt)),'Statistic','mean');
t4change = t4change + i_t3+round(1/dt);
i_t4linear = t4change;

%% Plot results
xt(1) = mean(Tdata(i_ts:i_t1));
xt(2) = mean(Tdata(i_t1:i_t2l));
xt(3) = Tdata(i_t2l+2);
xt(4) = mean(Tdata(i_t3:i_t4));

yt_BP(1:4) = max(SPdata)+5;
yt_HR(1:4) = max(HRdata)+5;

if PthAvail ~= 0
    yt_Pth(1:4) = max(Pthdata)+5;
end

str = {'I','II','III','IV'};

f = figure(2);clf;
set(gcf,'units','points','position',figSize);
f.Units = 'pixels';
sgtitle(patient,'Fontsize',fs,'FontWeight','bold','Interpreter', 'none')

if PthAvail ~= 0 
    subplot(4,1,1); hold on
        p1 = plot(Traw,Praw,'Color',[0 0.4470 0.7410],'DisplayName','BP');
        p2 = plot(Tdata,SPdata,'b','LineWidth',lw2,'DisplayName','Systolic BP');
        xline(Tdata(i_ts),'-',' ','LineWidth',1,'LabelVerticalAlignment', 'top',...
            'LabelHorizontalAlignment','right','FontSize',fs)
        xline(Tdata(i_t2l),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t1),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t2e),'--',' ','LineWidth',1)
        xline(Tdata(i_t3),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t4),'-','','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_PRT),'--','','LineWidth',1,'FontSize',fs)
        yline(SPbar,'--','Baseline before','linewidth',1,'FontSize',fs)
        text(xt,yt_BP,str,'FontSize',fs)
        set(gca,'Fontsize',fs);
        axis([Tdata(i_r1) Tdata(i_r2) min(Praw)-5 max(SPdata)+10]);
        ylabel('BP (mmHg)');
        legend([p1 p2],{'BP','SBP'});
        xlim([Tdata(i_r1) Tdata(i_r2)]);
        ylim([min(Praw)-5 max(Praw)+5])
    
    subplot(4,1,2); hold on
        plot(Tdata,HRdata,'b','LineWidth',lw2)
        xline(Tdata(i_ts),'-',' ','LineWidth',1,'LabelVerticalAlignment', 'top',...
            'LabelHorizontalAlignment','right','FontSize',fs)
        xline(Tdata(i_t2l),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t1),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t2e),'--',' ','LineWidth',1)
        xline(Tdata(i_t3),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t4),'-','','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_PRT),'--','','LineWidth',1,'FontSize',fs)
        text(xt,yt_HR,str,'FontSize',fs)
        set(gca,'Fontsize',fs);
        xlim([Tdata(i_r1) Tdata(i_r2)]);
        ylim([min(HRdata)-2 max(HRdata)+2]);
        ylabel('HR (bpm)');

     subplot(4,1,3); hold on    % using Traw with Pthdata bc they are same length
        plot(Tdata,Pthdata,'b','LineWidth',lw2)
        xline(Tdata(i_ts),'-',' ','LineWidth',1,'LabelVerticalAlignment', 'top',...
            'LabelHorizontalAlignment','right','FontSize',fs)
        xline(Tdata(i_t2l),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t1),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t2e),'--',' ','LineWidth',1)
        xline(Tdata(i_t3),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t4),'-','','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_PRT),'--','','LineWidth',1,'FontSize',fs)
        text(xt,yt_Pth,str,'FontSize',fs)
        set(gca,'Fontsize',fs);
        xlim([Tdata(i_r1) Tdata(i_r2)]);
        ylim([min(Pthdata)-2 max(Pthdata)+2]);
        ylabel('Pth (mmHg)');
        
    subplot(4,1,4); hold on;
        xline(Tdata(i_ts),'-',' ','LineWidth',1,'LabelVerticalAlignment', 'top',...
            'LabelHorizontalAlignment','right','FontSize',fs)
        xline(Tdata(i_t2l),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t1),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t2e),'--',' ','LineWidth',1)
        xline(Tdata(i_t3),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t4),'-','','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_PRT),'--','','LineWidth',1,'FontSize',fs)
        plot(Tdata,Rdata,'b','LineWidth',lw2)
        set(gca,'Fontsize',fs); 
        xlabel('Time (s)')
        ylabel('Respiration')
        xlim([Tdata(i_r1) Tdata(i_r2)]);

else
    subplot(3,1,1); hold on
        p1 = plot(Traw,Praw,'Color',[0 0.4470 0.7410],'DisplayName','BP');
        p2 = plot(Tdata,SPdata,'b','LineWidth',lw2,'DisplayName','Systolic BP');
        xline(Tdata(i_ts),'-',' ','LineWidth',1,'LabelVerticalAlignment', 'top','LabelHorizontalAlignment','right','FontSize',fs)
        xline(Tdata(i_t2l),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t1),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t2e),'--',' ','LineWidth',1)
        xline(Tdata(i_t3),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t4),'-','','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_PRT),'--','','LineWidth',1,'FontSize',fs)
        yline(SPbar,'--','Baseline before','linewidth',1,'FontSize',fs)
        text(xt,yt_BP,str,'FontSize',fs)
        set(gca,'Fontsize',fs);
        axis([Tdata(i_r1) Tdata(i_r2) min(Praw)-5 max(SPdata)+10]);
        ylabel('BP (mmHg)');
        legend([p1 p2],{'BP','SBP'});
        xlim([Tdata(i_r1) Tdata(i_r2)]);
        ylim([min(Praw)-5 max(Praw)+5])
    
    subplot(3,1,2); hold on
        plot(Tdata,HRdata,'b','LineWidth',lw2)
        xline(Tdata(i_ts),'-',' ','LineWidth',1,'LabelVerticalAlignment', 'top','LabelHorizontalAlignment','right','FontSize',fs)
        xline(Tdata(i_t2l),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t1),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t2e),'--',' ','LineWidth',1)
        xline(Tdata(i_t3),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t4),'-','','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_PRT),'--','','LineWidth',1,'FontSize',fs)
        text(xt,yt_HR,str,'FontSize',fs)
        set(gca,'Fontsize',fs);
        xlim([Tdata(i_r1) Tdata(i_r2)]);
        ylim([min(HRdata)-2 max(HRdata)+2]);
        ylabel('HR (bpm)');
        
    subplot(3,1,3); hold on;
        xline(Tdata(i_ts),'-',' ','LineWidth',1,'LabelVerticalAlignment', 'top','LabelHorizontalAlignment','right','FontSize',fs)
        xline(Tdata(i_t2l),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t1),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t2e),'--',' ','LineWidth',1)
        xline(Tdata(i_t3),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t4),'-','','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_PRT),'--','','LineWidth',1,'FontSize',fs)
        plot(Tdata,Rdata,'b','LineWidth',lw2)
        set(gca,'Fontsize',fs); 
        xlabel('Time (s)')
        ylabel('Respiration')
        xlim([Tdata(i_r1) Tdata(i_r2)]);
end

% Correct other time indices
answer = questdlg('Accept indices?','Index Correction','Yes','No ','Yes');

% Handle response
while strcmp(answer,'No ')==1
    prompt = {'Choose index to correct.',...
    'Hold control/command to select multiple.',...
    'Cancel to accept indices.'};

    indices = {'Phase I start','Early Phase II start','Late Phase II start',...
    'Phase III start','Phase IV start','Phase IV end'};

    [indx,~] = listdlg('PromptString',prompt,...
    'SelectionMode','multiple','ListString',indices,...
    'ListSize',[200,300],'Name','Index Correction',...
    'CancelString',{'Cancel'});
    choice = indices(indx);
        
    for i = 1:length(choice)
        switch choice{i}
            case 'Phase I start'
                if PthAvail ~= 0
                    sgtitle('Click on start of breath hold (Phase I start)');
                    title('- Phase I starts at the minimum before the first peak');
                    [tsc, ~ ] = ginput(1);
                    if isempty(tsc) == 0
                        [~,i_tsc] = min(abs(Tdata - tsc));
                    end
                    subplot(4,1,1); hold on
                    xline(Tdata(i_tsc),'r-','LineWidth',2)

                    subplot(4,1,2); hold on
                    xline(Tdata(i_tsc),'r-','LineWidth',2)
                    
                    subplot(4,1,3); hold on
                    xline(Tdata(i_tsc),'r-','LineWidth',2)
                    
                    subplot(4,1,4); hold on
                    xline(Tdata(i_tsc),'r-','LineWidth',2)
                else
                    sgtitle('Click on start of breath hold (Phase I start)');
                    title('- Phase I starts at the minimum before the first peak');
                    [tsc, ~ ] = ginput(1);
                    if isempty(tsc) == 0
                        [~,i_tsc] = min(abs(Tdata - tsc));
                    end
                    subplot(3,1,1); hold on
                    xline(Tdata(i_tsc),'r-','LineWidth',2)
                    subplot(3,1,2); hold on
                    xline(Tdata(i_tsc),'r-','LineWidth',2)
                    subplot(3,1,3); hold on
                    xline(Tdata(i_tsc),'r-','LineWidth',2)
                end
                
                ts = tsc;
                i_ts = i_tsc;
                
            case 'Early Phase II start'
                if PthAvail ~= 0
                    sgtitle('Click on start of early phase II');
                    title(' - Phase II starts when SBP peaks after start of breath hold');
                    [t1c, ~ ] = ginput(1);
                    if isempty(t1c) == 0
                        [~,i_t1c] = min(abs(Tdata - t1c));
                    end
                    subplot(4,1,1); hold on
                    xline(Tdata(i_t1c),'r-','LineWidth',2)
                    subplot(4,1,2); hold on
                    xline(Tdata(i_t1c),'r-','LineWidth',2)
                    subplot(4,1,3); hold on
                    xline(Tdata(i_t1c),'r-','LineWidth',2)
                    subplot(4,1,4); hold on
                    xline(Tdata(i_t1c),'r-','LineWidth',2)
                
                else
                    sgtitle('Click on start of early phase II');
                    title(' - Phase II starts when SBP peaks after start of breath hold');
                    [t1c, ~ ] = ginput(1);
                    if isempty(t1c) == 0
                        [~,i_t1c] = min(abs(Tdata - t1c));
                    end
                    subplot(3,1,1); hold on
                    xline(Tdata(i_t1c),'r-','LineWidth',2)
                    subplot(3,1,2); hold on
                    xline(Tdata(i_t1c),'r-','LineWidth',2)
                    subplot(3,1,3); hold on
                    xline(Tdata(i_t1c),'r-','LineWidth',2)
                end 

                t1 = t1c;
                i_t1 = i_t1c;
                
            case 'Late Phase II start'
                if PthAvail ~= 0
                    sgtitle('Click on start of late phase II');
                    title(' - minimum SBP in phase II');
                    [t2ec, ~ ] = ginput(1);
                    if isempty(t2ec) == 0
                        [~,i_t2ec] = min(abs(Tdata - t2ec));
                    end
                    subplot(4,1,1); hold on
                    xline(Tdata(i_t2ec),'r--','LineWidth',2)
                    subplot(4,1,2); hold on 
                    xline(Tdata(i_t2ec),'r--','LineWidth',2)
                    subplot(4,1,3); hold on 
                    xline(Tdata(i_t2ec),'r--','LineWidth',2)
                    subplot(4,1,4); hold on
                    xline(Tdata(i_t2ec),'r--','LineWidth',2)
                else
                    sgtitle('Click on start of late phase II');
                    title(' - minimum SBP in phase II');
                    [t2ec, ~ ] = ginput(1);
                    if isempty(t2ec) == 0
                        [~,i_t2ec] = min(abs(Tdata - t2ec));
                    end
                    subplot(3,1,1); hold on
                    xline(Tdata(i_t2ec),'r--','LineWidth',2)
                    subplot(3,1,2); hold on 
                    xline(Tdata(i_t2ec),'r--','LineWidth',2)
                    subplot(3,1,3); hold on
                    xline(Tdata(i_t2ec),'r--','LineWidth',2)
                end
                
                t2e = t2ec;
                i_t2e = i_t2ec;
                
            case 'Phase III start'
                if PthAvail ~= 0
                    sgtitle('Click on end of breath hold (start of Phase III)');
                    title(' - Phase III starts at the peak before transient drop in BP (~15 seconds after start of breath hold)');
                    [t2lc, ~ ] = ginput(1);
                    if isempty(t2lc) == 0
                        [~,i_t2lc] = min(abs(Tdata - t2lc));
                    end
                    subplot(4,1,1); hold on
                    xline(Tdata(i_t2lc),'r-','LineWidth',lw)
                    subplot(4,1,2); hold on
                    xline(Tdata(i_t2lc),'r-','LineWidth',lw)
                    subplot(4,1,3); hold on
                    xline(Tdata(i_t2lc),'r-','LineWidth',lw)
                    subplot(4,1,4); hold on
                    xline(Tdata(i_t2lc),'r-','LineWidth',lw)
                else
                    sgtitle('Click on end of breath hold (start of Phase III)');
                    title(' - Phase III starts at the peak before transient drop in BP (~15 seconds after start of breath hold)');
                    [t2lc, ~ ] = ginput(1);
                    if isempty(t2lc) == 0
                        [~,i_t2lc] = min(abs(Tdata - t2lc));
                    end
                    subplot(3,1,1); hold on
                    xline(Tdata(i_t2lc),'r-','LineWidth',lw)
                    subplot(3,1,2); hold on
                    xline(Tdata(i_t2lc),'r-','LineWidth',lw)
                    subplot(3,1,3); hold on
                    xline(Tdata(i_t2lc),'r-','LineWidth',lw)
                end
                
                i_t2l = i_t2lc;
                
            case 'Phase IV start'
                if PthAvail ~= 0
                    sgtitle('Click on start of phase IV');
                    title(' - Phase IV starts at the minimum afer phase III');
                    [t3c, ~ ] = ginput(1);
                    if isempty(t3c) == 0
                        [~,i_t3c] = min(abs(Tdata - t3c));
                    end
                    subplot(4,1,1); hold on
                    xline(Tdata(i_t3c),'r-','LineWidth',2)
                    subplot(4,1,2); hold on
                    xline(Tdata(i_t3c),'r-','LineWidth',2) 
                    subplot(4,1,3); hold on
                    xline(Tdata(i_t3c),'r-','LineWidth',2) 
                    subplot(4,1,4); hold on
                    xline(Tdata(i_t3c),'r-','LineWidth',2) 
                else
                    sgtitle('Click on start of phase IV');
                    title(' - Phase IV starts at the minimum afer phase III');
                    [t3c, ~ ] = ginput(1);
                    if isempty(t3c) == 0
                        [~,i_t3c] = min(abs(Tdata - t3c));
                    end
                    subplot(3,1,1); hold on
                    xline(Tdata(i_t3c),'r-','LineWidth',2)
                    subplot(3,1,2); hold on
                    xline(Tdata(i_t3c),'r-','LineWidth',2) 
                    subplot(3,1,3); hold on
                    xline(Tdata(i_t3c),'r-','LineWidth',2) 
                end
                
                t3 = t3c;
                i_t3 = i_t3c;
                N = length(SPdata);
                if SPdata(i_t3) > SPbar
                   k_PRT = find(SPdata(i_t3+1:i_t3+round(5/dt)) > SPdata(i_t2l),1,'first');
                else
                   k_PRT = find(SPdata(i_t3+1:end) > SPbar,1,'first');
                end
                if isempty(k_PRT) == 1
                    i_PRT = N;
                else
                    i_PRT = k_PRT + i_t3;
                end

            case 'Phase IV end'
                if PthAvail ~= 0
                    sgtitle('Click on end of phase IV');
                    title('- Phase IV ends after the overshoot returns to baseline');
                    [t4c, ~ ] = ginput(1);
                    if isempty(t4c) == 0
                        [~,i_t4c] = min(abs(Tdata - t4c));
                    end
                    subplot(4,1,1); hold on
                    xline(Tdata(i_t4c),'r-','LineWidth',2)
                    subplot(4,1,2); hold on
                    xline(Tdata(i_t4c),'r-','LineWidth',2) 
                    subplot(4,1,3); hold on
                    xline(Tdata(i_t4c),'r-','LineWidth',2) 
                    subplot(4,1,4); hold on
                    xline(Tdata(i_t4c),'r-','LineWidth',2) 
                else
                    sgtitle('Click on end of phase IV');
                    title('- Phase IV ends after the overshoot returns to baseline');
                    [t4c, ~ ] = ginput(1);
                    if isempty(t4c) == 0
                        [~,i_t4c] = min(abs(Tdata - t4c));
                    end
                    subplot(3,1,1); hold on
                    xline(Tdata(i_t4c),'r-','LineWidth',2)
                    subplot(3,1,2); hold on
                    xline(Tdata(i_t4c),'r-','LineWidth',2) 
                    subplot(3,1,3); hold on
                    xline(Tdata(i_t4c),'r-','LineWidth',2) 
                end
                
                t4 = t4c;
                i_t4 = i_t4c;
        end
    end
    
    xt(1) = mean(Tdata(i_ts:i_t1));
    xt(2) = mean(Tdata(i_t1:i_t2l));
    xt(3) = Tdata(i_t2l+2);
    xt(4) = mean(Tdata(i_t3:i_t4));

    % add if statement for PthAvail and repeat 
    f = figure(2);clf
    set(gcf,'units','points','position',figSize)
    f.Units = 'pixels';
    sgtitle(patient,'Fontsize',fs,'FontWeight','bold','Interpreter', 'none')
    
    if PthAvail ~= 0 
        subplot(4,1,1); hold on
            p1 = plot(Traw,Praw,'Color',[0 0.4470 0.7410],'DisplayName','BP');
            p2 = plot(Tdata,SPdata,'b','LineWidth',2,'DisplayName','Systolic BP');
            xline(Tdata(i_ts),'-',' ','LineWidth',1,'LabelVerticalAlignment', 'top',...
                'LabelHorizontalAlignment','right','FontSize',fs)
            xline(Tdata(i_t2l),'-',' ','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_t1),'-',' ','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_t2e),'--',' ','LineWidth',1)
            xline(Tdata(i_t3),'-',' ','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_t4),'-','','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_PRT),'--','','LineWidth',1,'FontSize',fs)
            yline(SPbar,'--','Baseline before','linewidth',1,'FontSize',fs)
            text(xt,yt_BP,str,'FontSize',fs)
            axis([Tdata(i_r1) Tdata(i_r2) min(Praw)-5 max(SPdata)+10]);
            ylabel('BP (mmHg)');
            xlim([Tdata(i_r1) Tdata(i_r2)]);
            set(gca,'Fontsize',fs);
    
        subplot(4,1,2); hold on
            plot(Tdata,HRdata,'b','LineWidth',lw)
            xline(Tdata(i_ts),'-',' ','LineWidth',1,'LabelVerticalAlignment', 'top',...
                'LabelHorizontalAlignment','right','FontSize',fs)
            xline(Tdata(i_t2l),'-',' ','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_t1),'-',' ','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_t2e),'--',' ','LineWidth',1)
            xline(Tdata(i_t3),'-',' ','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_t4),'-','','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_PRT),'--','','LineWidth',1,'FontSize',fs)
            xlim([Tdata(i_r1) Tdata(i_r2)]);
            text(xt,yt_HR,str,'FontSize',fs)
            ylim([min(HRdata)-5 max(HRdata)+10]);
            ylabel('HR (bpm)');
            set(gca,'Fontsize',fs);

        subplot(4,1,3); hold on
            plot(Tdata,Pthdata,'b','LineWidth',lw)
            xline(Tdata(i_ts),'-',' ','LineWidth',1,'LabelVerticalAlignment', 'top',...
                'LabelHorizontalAlignment','right','FontSize',fs)
            xline(Tdata(i_t2l),'-',' ','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_t1),'-',' ','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_t2e),'--',' ','LineWidth',1)
            xline(Tdata(i_t3),'-',' ','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_t4),'-','','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_PRT),'--','','LineWidth',1,'FontSize',fs)
            xlim([Tdata(i_r1) Tdata(i_r2)]);
            text(xt,yt_Pth,str,'FontSize',fs)
            ylim([min(Pthdata)-5 max(Pthdata)+10]);
            ylabel('Pth (mmHg)');
            set(gca,'Fontsize',fs);
    
        subplot(4,1,4); hold on;
            set(gca,'Fontsize',fs); 
            xline(Tdata(i_ts),'-',' ','LineWidth',1,'LabelVerticalAlignment', 'top',...
                'LabelHorizontalAlignment','right','FontSize',fs)
            xline(Tdata(i_t2l),'-',' ','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_t1),'-',' ','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_t2e),'--',' ','LineWidth',1)
            xline(Tdata(i_t3),'-',' ','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_t4),'-','','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_PRT),'--','','LineWidth',1,'FontSize',fs)
            plot(Tdata,Rdata,'b','LineWidth',lw)
            xlim([Tdata(i_r1) Tdata(i_r2)]);
            xlabel('Time (s)')
            ylabel('Respiration')
    else
        subplot(3,1,1); hold on
            p1 = plot(Traw,Praw,'Color',[0 0.4470 0.7410],'DisplayName','BP');
            p2 = plot(Tdata,SPdata,'b','LineWidth',2,'DisplayName','Systolic BP');
            xline(Tdata(i_ts),'-',' ','LineWidth',1,'LabelVerticalAlignment', 'top',...
                'LabelHorizontalAlignment','right','FontSize',fs)
            xline(Tdata(i_t2l),'-',' ','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_t1),'-',' ','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_t2e),'--',' ','LineWidth',1)
            xline(Tdata(i_t3),'-',' ','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_t4),'-','','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_PRT),'--','','LineWidth',1,'FontSize',fs)
            yline(SPbar,'--','Baseline before','linewidth',1,'FontSize',fs)
            text(xt,yt_BP,str,'FontSize',fs)
            axis([Tdata(i_r1) Tdata(i_r2) min(Praw)-5 max(SPdata)+10]);
            ylabel('BP (mmHg)');
            xlim([Tdata(i_r1) Tdata(i_r2)]);
            set(gca,'Fontsize',fs);
    
        subplot(3,1,2); hold on
            plot(Tdata,HRdata,'b','LineWidth',lw)
            xline(Tdata(i_ts),'-',' ','LineWidth',1,'LabelVerticalAlignment', 'top',...
                'LabelHorizontalAlignment','right','FontSize',fs)
            xline(Tdata(i_t2l),'-',' ','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_t1),'-',' ','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_t2e),'--',' ','LineWidth',1)
            xline(Tdata(i_t3),'-',' ','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_t4),'-','','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_PRT),'--','','LineWidth',1,'FontSize',fs)
            xlim([Tdata(i_r1) Tdata(i_r2)]);
            text(xt,yt_HR,str,'FontSize',fs)
            ylim([min(HRdata)-5 max(HRdata)+10]);
            ylabel('HR (bpm)');
            set(gca,'Fontsize',fs);
    
        subplot(3,1,3); hold on;
            set(gca,'Fontsize',fs); 
            xline(Tdata(i_ts),'-',' ','LineWidth',1,'LabelVerticalAlignment', 'top',...
                'LabelHorizontalAlignment','right','FontSize',fs)
            xline(Tdata(i_t2l),'-',' ','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_t1),'-',' ','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_t2e),'--',' ','LineWidth',1)
            xline(Tdata(i_t3),'-',' ','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_t4),'-','','LineWidth',1,'FontSize',fs)
            xline(Tdata(i_PRT),'--','','LineWidth',1,'FontSize',fs)
            plot(Tdata,Rdata,'b','LineWidth',lw)
            xlim([Tdata(i_r1) Tdata(i_r2)]);
            xlabel('Time (s)')
            ylabel('Respiration')
    end

    answer = questdlg('Accept indices?', ...
        'Index Correction', ...
        'Yes','No ','Yes');
    msg_not_equal = 'Indices must be in ascending order';
    
    if i_ts < i_t1 == 0 
        uiwait(msgbox(msg_not_equal,'Error',"modal"));
        uiwait(msgbox('Phase I start > Phase II start','Fix',"modal"));
        answer = 'No ';
    elseif i_t1  < i_t2e == 0
        uiwait(msgbox(msg_not_equal,'Error',"modal"));
        uiwait(msgbox('Phase II start > Middle of Phase II','Fix',"modal"));
        answer = 'No ';
    elseif i_t2e < i_t2l == 0 
        uiwait(msgbox(msg_not_equal,'Error',"modal"));
        uiwait(msgbox('Middle of Phase II > Start of Phase III','Fix',"modal"));
        answer = 'No ';
    elseif i_t2l < i_t3  == 0
        uiwait(msgbox(msg_not_equal,'Error',"modal"));
        uiwait(msgbox('Start of Phase III > Start of Phase IV','Fix',"modal"));
        answer = 'No ';
    elseif i_t3  < i_t4  == 0
        msg_not_equal = 'Indices must be in ascending order';
        uiwait(msgbox(msg_not_equal,'Error',"modal"));
        uiwait(msgbox('Start of Phase IV > End of Phase IV','Fix',"modal"));
        answer = 'No ';
    end
end

%% Output
phases = [i_ts; i_t1linear; i_t1; i_t2elinear; i_t2emin; i_t2e; i_t2l; i_t3; i_t3linear; i_PRT; i_t4; i_t4linear];

end % function set_indices.m %