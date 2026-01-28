function clinicalratios(patient,plotmarkers,patids)
% function clinicalratios 
% Input: data and plotmarkers (structure with figure and screen settings)
% Output: table with clinical ratios, and amended data structure. This
% function also saves figures to the Figure directory and an excel file
% with the patient markers. Stored in the "Markers" folder.
% Uses: closeGenButton to close figures
% Description: Runs software

if nargin == 3
    load(strcat('../WS/',patids,patient,'_WS_temp.mat'), 'data','pat_char','raw_data');
else
    load(strcat('../WS/',patient,'_WS_temp.mat'), 'data','pat_char','raw_data');
end

% Values extracted from data
Tdata     = data.Tdata;  % Time (vector)
SPdata    = data.SPdata; % Systolic blood pressure (vector)
DPdata    = data.DPdata; % Diastolic blood pressure (vector)
PPdata    = SPdata - DPdata;  % Pulse pressure (vector)
HRdata    = data.HRdata; % Heart rate (vector)
RRdata    = data.RRdata; % RR intervals (vector)
SP_mean   = data.SPbar;  % Mean systolic blood pressure before VM (value)
HRb_mean  = data.Hbar;   % Mean heart rate before VM (value)
dt        = data.dt;     % Equidistant distance between Tdata points (value)

% Indices indicating VM phases
i_ts        = data.i_ts;    % Start  of phase 1
i_t1linear  = data.i_t1linear; % beginning of linear regression for phase 1
i_t1        = data.i_t1;    % End of phase 1
i_t2elinear = data.i_t2elinear; % point of maximum change in early Phase 2
i_t2emin    = data.i_t2emin; % Minimum point of Phase 2
i_t2e       = data.i_t2e;   % End of early phase 2
i_t2l       = data.i_t2l;   % End of late phase 2
i_t3        = data.i_t3;    % End of phase 3
i_t3linear  = data.i_t3linear; 
i_t4        = data.i_t4;    % End of phase 4 
i_t4linear  = data.i_t4linear;
i_PRT       = data.i_PRT;   % End of pressure recovery time

% Interval lengths (seconds)
T_1   = Tdata(i_t1)  - Tdata(i_ts);  % length of phase 1
T_2e  = Tdata(i_t2e) - Tdata(i_t1);  % length of early phase 2
T_2l  = Tdata(i_t2l) - Tdata(i_t2e); % length of late phase 2
T_3   = Tdata(i_t3)  - Tdata(i_t2l); % length of phase 3
T_PRT = Tdata(i_PRT) - Tdata(i_t3);  % pressure recovery time
T_4   = Tdata(i_t4)  - Tdata(i_t3);  % length of phase 4


%% Systolic blood pressures
% Phase I (maximum)
i_SP1_Max = i_t1;
i_SP1_Max_reg = i_t1linear;
SP1_Max   = SPdata(i_SP1_Max);
SP1_Max_reg = SPdata(i_SP1_Max_reg);

PP1_Max = PPdata(i_SP1_Max);

% Early phase II (minimum)
i_SP2e = i_t2e;
SP2e = SPdata(i_SP2e);
i_SP2e_min = i_t2emin;
SP2e_min   = SPdata(i_SP2e_min);
i_SP2e_linear = i_t2elinear;
SP2e_linear = SPdata(i_t2elinear);

% Late phase II (maximum and end)
SP2l_end  = SPdata(i_t2l);
[SP2l_Max, kt2l_Max] = max(SPdata(i_t2e:i_t2l));
i_t2l_Max = kt2l_Max + i_t2e - 1;

SP3_linear = SPdata(i_t3linear);

% Phase III (minimum)
SP3_end   = SPdata(i_t3);

% Early phase IV (within 5 seconds, minimum)
[SP4e_Max, i_SP4e_Max] = max(SPdata(i_t3:i_t3+round(6/dt)));
i_SP4e_Max = i_SP4e_Max + i_t3 - 1;
PP4_Max = PPdata(i_SP4e_Max);

i_SP4e_Max_reg = i_t4linear;
SP4e_Max_reg = SPdata(i_t4linear);


%% HR and RR intervals
% min HR (Max RR) early phase 2
[HR2e_min,k_HR2e_min] = min(HRdata(i_t1:i_t2e));
i_HR2e_min  = k_HR2e_min + i_t1 - 1;
RR2e_Max    = 60/HR2e_min;

% HR (RR) early phase 2
HR2e_end = HRdata(i_t2e);
RR2e_end = 60/HR2e_end;

% Max HR (min RR) early phase 4
[HR4e_Max,k_HR4e_Max] = max(HRdata(i_t3:i_t3+round(10/dt)));
i_HR4e_Max  = k_HR4e_Max + i_t3 - 1;
RR4e_min    = 60/HR4e_Max;
[HR_Max,i_HR_Max] = max(HRdata);
RR_min       = 60/HR_Max;

HRsmooth = sgolayfilt(HRdata,3,7);
HR1st = gradient(HRsmooth,dt);

% If the maximum heart rate is not between phase 2 early and phase 4, 
% use a relative maximum between phase 3 and 10 seconds later
[HR_Max,k_HR_Max] = max(HRdata(i_t3:i_t3+round(10/dt)));
i_HR_Max = k_HR_Max + i_t3 - 1;

% Find the changepoint where it is close enough to the cliff
HR4e_offset = round(0.7/dt);
i_HR4e_Max = findchangepts(HR1st(i_HR_Max:i_HR_Max+HR4e_offset),'Statistic','mean') + i_HR_Max;
while ~(HRdata(i_HR4e_Max) >= 1.1*HRdata(i_HR4e_Max + round(1/dt)))
    HR4e_offset = HR4e_offset + round(0.5/dt);
    if HR4e_offset > length(HRdata) - i_HR_Max
        i_HR4e_Max = findchangepts(HR1st(i_HR_Max:i_HR_Max+round(2/dt)),'Statistic','mean') + i_HR_Max;
        break;
    else
        i_HR4e_Max = findchangepts(HR1st(i_HR_Max:i_HR_Max+HR4e_offset),'Statistic','mean') + i_HR_Max;
    end
end
HR4e_Max = HRdata(i_HR4e_Max);

RR4e_linStart = 60/HR4e_Max;

% min HR (max RR) early phase 4 
[HR4e_min,k_HR4e_min] = min(HRdata(i_HR4e_Max:i_HR4e_Max+round(3/dt))); 
i_HR4e_min = k_HR4e_min + i_HR4e_Max - 1;
RR4e_Max   = 60/HR4e_min;

i_HR4e_min_reg = findchangepts(HR1st(i_HR4e_Max:i_HR4e_Max + round(3.5/dt)),'Statistic','mean') + i_HR4e_Max;
HR4e_min_reg = HRdata(i_HR4e_min_reg);
RR4e_min_reg = 60/HR4e_min_reg;

% min HR (Max RR) phase 4
if length(Tdata) > i_HR4e_Max+round(30/dt)
    [HR4_min,k_HR4_min] = min(HRdata(i_HR4e_Max:i_HR4e_Max+round(30/dt)));
else
    [HR4_min,k_HR4_min] = min(HRdata(i_HR4e_Max:end));
end
i_HR4_min = k_HR4_min + i_HR4e_Max-1;
RR4_Max   = 60/HR4_min;

% Mean HR and RR intervals 
HRa_mean= mean(HRdata(i_HR4e_min:end)); % After VM
RRb_mean = 60./HRb_mean*1000; % Before VM (ms)
RRa_mean = 60./HRa_mean*1000; % After VM (ms)

%% Pulse pressure

% Mean PP rest before VM
PPbar_b = mean(PPdata(1:i_ts));

% Phase I (maximum)
PP1_Max   = PPdata(i_SP1_Max);      % Pulse pressure at max pressure in phase 1 (at maximum pressure)
PP1_Max2  = max(PPdata(i_ts:i_t1)); % Max pulse pressure phase 1

% Early phase II (minimum)
PP2e_min  = PPdata(i_SP2e_min);       % Pulse pressure at min pressure early phase 2
PP2e_min2 = min(PPdata(i_t1:i_t2e));  % Min pulse pressure phase 2 early

dVPP  = 100*(PP1_Max-PP2e_min)/((PP1_Max+PP2e_min)/2);
dVPP1 = 100*(PP1_Max2-PP2e_min2)/((PP1_Max2+PP2e_min2)/2);
dVSP  = 100*(SPdata(i_SP1_Max)-SPdata(i_SP2e_min))/((SPdata(i_SP1_Max)+SPdata(i_SP2e_min))/2);

% Early phase IV (Max)
%PP4e_Max = max(PPdata(i_t3:i_t3+round(5/dt)));

% Mean PP rest after VM
PPbar_a = mean(PPdata(i_HR4e_min:end));


%% Figure properties
fs   = plotmarkers.fs;
fsL  = fs + 3;
ms   = plotmarkers.ms;
lw   = plotmarkers.lwt;
if lw > 2
    lw2 = lw -1;
else
    lw2 = 1;
end
figSize = plotmarkers.figSize;

% Time intervals used for regression lines
xt(1) = mean(Tdata(i_ts:i_t1));
xt(2) = mean(Tdata(i_t1:i_t2l));
xt(3) = Tdata(i_t2l+2);
xt(4) = mean(Tdata(i_t3:i_t4));
xtt   = mean(Tdata);

% Axis limits
HRmin = floor(min(HRdata)/5)*5-2;
HRmax = ceil(max(HRdata)/5)*5+2;
SPmin = floor(min(data.Praw)/10)*10-2;
SPmax = ceil(max(SPdata)/10)*10+10;
RRmin = floor(min(RRdata)/10)*10-2;
RRmax = ceil(max(RRdata)/10)*10+10;

yt_BP1(1:4) = SPmax + 6; 
yt_BP2      = SPmax + 16; 
str1 = {'I','II','III','IV'};
str2 = patient;

f = figure(1); hold on;
set(gcf,'units','points','position',figSize)
f.Units = 'pixels';
subplot(3,1,1); hold on
   p11 = plot(Tdata(i_SP1_Max),SP1_Max,'co','Markersize',ms,'Linewidth',lw2);
   p11a= plot(Tdata(i_SP1_Max_reg),SP1_Max_reg,'c*','MarkerSize',ms-2,'LineWidth',lw2);
   p21 = plot(Tdata(i_SP2e_min),SP2e_min,'o','color','#0000FF','Markersize',ms,'Linewidth',lw2);
   p21a= plot(Tdata(i_SP2e_linear),SP2e_linear,'*','color','#0000FF','MarkerSize',ms,'LineWidth',lw2);
   p31 = plot(Tdata(i_t2l_Max),SP2l_Max,'o','color','#7881C2','Markersize',ms,'Linewidth',lw2);
   p41 = plot(Tdata(i_t2l),SP2l_end,'ko','Markersize',ms,'Linewidth',lw2);
   p51 = plot(Tdata(i_t3),SP3_end,'o','color','#D95319','Markersize',ms,'Linewidth',lw2);
   p51a= plot(Tdata(i_t3linear),SP3_linear,'*','color','#D95319','MarkerSize',ms,'LineWidth',lw2);
   p61 = plot(Tdata(i_SP4e_Max),SP4e_Max,'o','color','#851803','Markersize',ms,'Linewidth',lw2);
   p61a= plot(Tdata(i_SP4e_Max_reg),SP4e_Max_reg,'*','color','#851803','Markersize',ms,'Linewidth',lw2);
   
   text(xt, yt_BP1,str1,'FontSize',fsL);
   text(xtt,yt_BP2,str2,'FontSize',fsL,'Interpreter', 'none');
   plot(data.Traw,data.Praw,'Color',[.7 .7 .7]);
   plot(Tdata,SPdata,'k','linewidth',lw2);
   xline(Tdata(i_ts),'-',' ','LineWidth',1);
   xline(Tdata(i_t2l),'-',' ','LineWidth',1);
   xline(Tdata(i_t1),'-',' ','LineWidth',1);
   xline(Tdata(i_t2e),'--',' ','LineWidth',1);
   xline(Tdata(i_t3),'-',' ','LineWidth',1);
   xline(Tdata(i_t4),'-','','LineWidth',1);
   xline(Tdata(i_PRT),'--','','LineWidth',1);
   l1 = plot(Tdata,SP_mean*ones(size(Tdata)),'k--','linewidth',1);
   text(Tdata(1)+1,SP_mean+10,'Baseline before','fontsize',fs);
   
   xlim([min(Tdata),max(Tdata)]);
   ylimmin = 5*floor((min(data.Praw)-5)/5);
   ylimmax = 5*ceil((max(SPdata)+5)/5);   
   ylim([ylimmin ylimmax]);
   
   set(gca,'Fontsize',fs);
   ylabel('BP (mmHg)');
   pleglines = [p11,p11a,p21,p21a,p31,p41,p51,p51a,p61a,p61];
   plegNames = {'1 Max','1 Max reg','2e min','2e reg','2l Max','2l end','3 end','3 reg','4e Max reg','4 Max'};
   pleg=legend(pleglines,plegNames,'Location','eastoutside');

   
   subplot(3,1,2); hold on
   plot(Tdata,HRdata,'k','linewidth',lw2);
   p12 = plot(Tdata(i_HR2e_min),HR2e_min,'co','Markersize',ms,'Linewidth',lw2);
   p22 = plot(Tdata(i_t2e),HR2e_end,'o','color','#0000FF','Markersize',ms,'Linewidth',lw2);
   p32 = plot(Tdata(i_HR4e_Max),HR4e_Max,'o','color','#D95319','Markersize',ms,'Linewidth',lw2);
   p42a= plot(Tdata(i_HR4e_min_reg),HR4e_min_reg,'*','color','#851803','Markersize',ms,'Linewidth',lw2);
   p42 = plot(Tdata(i_HR4e_min),HR4e_min,'o','color','#851803','Markersize',ms,'Linewidth',lw2);
   p52 = plot(Tdata(i_HR4_min),HR4_min,'o','Color','#77AC30','Markersize',ms,'Linewidth',lw2);
   p62 = plot(Tdata(i_HR_Max),HR_Max,'mo','Markersize',ms,'Linewidth',lw2);
   xline(Tdata(i_ts),'-',' ','LineWidth',1);
   xline(Tdata(i_t2l),'-',' ','LineWidth',1);
   xline(Tdata(i_t1),'-',' ','LineWidth',1);
   xline(Tdata(i_t2e),'--',' ','LineWidth',1);
   xline(Tdata(i_t3),'-',' ','LineWidth',1);
   xline(Tdata(i_t4),'-','','LineWidth',1);
   xline(Tdata(i_PRT),'--','','LineWidth',1);
   ylabel('HR (bpm)');
   set(gca,'Fontsize',fs);
   xlim([min(Tdata),max(Tdata)]);
   ylimmin = 5*floor((min(HRdata)-2)/5);
   ylimmax = 5*ceil((max(HRdata)+2)/5);   
   ylim([ylimmin ylimmax]);
   
   l2=plot(Tdata(1:i_ts),HRb_mean*ones(size(Tdata(1:i_ts))),'k--','linewidth',1);
   text(Tdata(1)+1,HRb_mean+10,'Baseline before','fontsize',fs);
   l3=plot(Tdata(i_HR4e_min:end),HRa_mean*ones(size(Tdata(i_HR4e_min:end))),'k--','linewidth',1);
   text(Tdata(end)-5,HRa_mean+10,'Baseline after','fontsize',fs);
  
   hleglines = [p12,p22,p32,p42a,p42,p52,p62];
   hlegNames = {'2e min','2e end','4e Max','4e min reg','4e min','4 min','Max'};
   hleg=legend(hleglines,hlegNames,'Location','eastoutside');
 
   subplot(3,1,3); hold on
   plot(Tdata,RRdata,'k','linewidth',lw2);
   p13=plot(Tdata(i_HR2e_min),RR2e_Max,'co','Markersize',ms,'Linewidth',lw2);
   p23=plot(Tdata(i_t2e),RR2e_end,'o','color','#0000FF','Markersize',ms,'Linewidth',lw2);
   p33=plot(Tdata(i_HR4e_Max),RR4e_min,'o','color','#D95319','Markersize',ms,'Linewidth',lw2);
   p43a=plot(Tdata(i_HR4e_min_reg),RR4e_min_reg,'*','color','#851803','Markersize',ms,'Linewidth',lw2);
   p43=plot(Tdata(i_HR4e_min),RR4e_Max,'o','color','#851803','Markersize',ms,'Linewidth',lw2);
   p53=plot(Tdata(i_HR4_min),RR4_Max,'o','Color','#77AC30','Markersize',ms,'Linewidth',lw2);
   p63=plot(Tdata(i_HR_Max),RR_min,'mo','Markersize',ms,'Linewidth',lw2);
   xline(Tdata(i_ts), '-',' ','LineWidth',1);
   xline(Tdata(i_t2l),'-',' ','LineWidth',1);
   xline(Tdata(i_t1), '-',' ','LineWidth',1);
   xline(Tdata(i_t2e),'--',' ','LineWidth',1);
   xline(Tdata(i_t3), '-',' ','LineWidth',1);
   xline(Tdata(i_t4), '-','','LineWidth',1);
   xline(Tdata(i_PRT),'--','','LineWidth',1);
   ylabel('RR (s)');
   xlabel('Time (s)');
   set(gca,'Fontsize',fs);
   xlim([min(Tdata),max(Tdata)]);
   ylimmin = 0.1*floor((min(RRdata)-0.05)/0.1);
   ylimmax = 0.1*ceil ((max(RRdata)+0.05)/0.1);   
   ylim([ylimmin ylimmax]);
   
   l4=plot(Tdata(1:i_ts),RRb_mean*ones(size(Tdata(1:i_ts))),'k--','linewidth',1);
   text(Tdata(1)+1,RRb_mean+0.15,'Baseline before','fontsize',fs);
   l5=plot(Tdata(i_HR4e_min:end),RRa_mean*ones(size(Tdata(i_HR4e_min:end))),'k--','linewidth',1);
   text(Tdata(end)-5,RRa_mean+0.15,'Baseline after','fontsize',fs);

   rleglines = [p13,p23,p33,p43a,p43,p53,p63];
   rlegNames = {'2e Max','2e end','4e min','4e Max reg','4e Max','4 Max','min'};
   rleg=legend(rleglines,rlegNames,'Location','eastoutside');   
   
%% Point Correction
answer_CR = questdlg('Do you want to accept markers?', ...
    'Points of Interest','Yes','No ','Yes');

if strcmp(answer_CR,'Yes') == 0
    liststr = {'SP1 Max','SP1 Max reg','SP2e min','SP2e reg','SP2l Max','SP2l end',...
        'SP3 end','SP3 reg','SP4 reg','SP4e Max','---------','HR2e min','HR2e end'...
        'HR4e Max','HR4e min reg','HR4e min','HR4 min','HR Max'};
    prompt  = {'Select marker to move',...
        'Hold control/command to select multiple markers',''};
    [indx,~] = listdlg('PromptString',prompt,...
        'SelectionMode','multiple','ListString',liststr,...
        'ListSize',[200,300],'Name','Point Selection');
    list = liststr(indx);
    
    for i = 1:length(list)
        if strcmp(list{i},'---------') == 1
        elseif strcmp(list{i},'SP1 Max') == 1
            uiwait(msgbox("Select phase I max SBP","modal"));
            [new_ind,SP1_Max] = ginput(1);
            [~,i_SP1_Max] = min(abs(Tdata-new_ind));
            PP1_Max = PPdata(i_SP1_Max);
            figure(1); hold on;
            subplot(3,1,1); hold on;
            delete(p11)
            p11 = plot(Tdata(i_SP1_Max),SP1_Max,'co','Markersize',ms,'Linewidth',lw2);
            delete(pleg)
            pleglines = [p11,p11a,p21,p21a,p31,p41,p51,p51a,p61a,p61];  
            pleg=legend(pleglines,plegNames,'Location','eastoutside');

        elseif strcmp(list{i},'SP1 Max reg') == 1
            uiwait(msgbox("Select phase I max SBP for regression","modal"));
            [new_ind,SP1_Max_reg] = ginput(1);
            disp('before')
            disp([i_SP1_Max i_SP1_Max_reg]);

            [~,i_SP1_Max_reg] = min(abs(Tdata-new_ind));
            disp('after')
            disp([i_SP1_Max i_SP1_Max_reg]);

            figure(1); hold on;
            subplot(3,1,1); hold on;
            delete(p11a)
            p11a = plot(Tdata(i_SP1_Max_reg),SP1_Max_reg,'c*','Markersize',ms-2,'Linewidth',lw2);
            delete(pleg)
            pleglines = [p11,p11a,p21,p21a,p31,p41,p51,p51a,p61a,p61];  
            pleg=legend(pleglines,plegNames,'Location','eastoutside'); 
        
        elseif strcmp(list{i},'SP2e reg') == 1
            uiwait(msgbox("Select early phase II min SBP for regression","modal"));
            [new_ind,SP2e_linear] = ginput(1);
            [~,i_SP2e_linear] = min(abs(Tdata-new_ind));
            PP2e_linear  = PPdata(i_SP2e_linear);
            figure(1); hold on;
            subplot(3,1,1); hold on;
            delete(p21a)
            p21a=plot(Tdata(i_SP2e_linear),SP2e_linear,'*','color','#0000FF','Markersize',ms,'Linewidth',lw2);
            delete(pleg)
            pleglines = [p11,p11a,p21,p21a,p31,p41,p51,p51a,p61a,p61];  
            pleg=legend(pleglines,plegNames,'Location','eastoutside');

        elseif strcmp(list{i},'SP2e min') == 1
            uiwait(msgbox("Select early phase II min SBP","modal"));
            [new_ind,SP2e_min] = ginput(1);
            [~,i_SP2e_min] = min(abs(Tdata-new_ind));
            PP2e_min  = PPdata(i_SP2e_min);
            figure(1); hold on;
            subplot(3,1,1); hold on;
            delete(p21)
            p21=plot(Tdata(i_SP2e_min),SP2e_min,'o','color','#0000FF','Markersize',ms,'Linewidth',lw2);
            delete(pleg)
            pleglines = [p11,p11a,p21,p21a,p31,p41,p51,p51a,p61a,p61];  
            pleg=legend(pleglines,plegNames,'Location','eastoutside');

        elseif strcmp(list{i},'SP2l Max') == 1
            uiwait(msgbox("Select late phase II max SBP","modal"));
            [new_ind,SP2l_Max] = ginput(1);
            [~,i_t2l_Max] = min(abs(Tdata-new_ind));
            figure(1); hold on;
            subplot(3,1,1); hold on;
            delete(p31)
            p31=plot(Tdata(i_t2l_Max),SP2l_Max,'o','color','#7881C2','Markersize',ms,'Linewidth',lw2);
            delete(pleg)
            pleglines = [p11,p11a,p21,p21a,p31,p41,p51,p51a,p61a,p61];  
            pleg=legend(pleglines,plegNames,'Location','eastoutside');
            
        elseif strcmp(list{i},'SP2l end') == 1
            uiwait(msgbox("Select late phase II end SBP","modal"));
            [new_ind,SP2l_end] = ginput(1);
            [~,i_SP2l_end] = min(abs(Tdata-new_ind));
            figure(1); hold on;
            subplot(3,1,1); hold on;
            delete(p41)
            p41=plot(Tdata(i_SP2l_end),SP2l_end,'ko','Markersize',ms,'Linewidth',lw2);
            delete(pleg)
            pleglines = [p11,p11a,p21,p21a,p31,p41,p51,p51a,p61a,p61];  
            pleg=legend(pleglines,plegNames,'Location','eastoutside');
            
        elseif strcmp(list{i},'SP3 end') == 1
            uiwait(msgbox("Select phase III end SBP","modal"));
            [new_ind,SP3_end] = ginput(1);
            [~,i_SP3_end] = min(abs(Tdata-new_ind));
            figure(1); hold on;
            subplot(3,1,1); hold on;
            delete(p51);
            p51=plot(Tdata(i_SP3_end),SP3_end,'o','color','#D95319','Markersize',ms,'Linewidth',lw2);
            delete(pleg);
            pleglines = [p11,p11a,p21,p21a,p31,p41,p51,p51a,p61a,p61];  
            pleg=legend(pleglines,plegNames,'Location','eastoutside');

        elseif strcmp(list{i},'SP3 reg') == 1
            uiwait(msgbox("Select phase III end SBP for regression","modal"));
            [new_ind,SP3_linear] = ginput(1);
            [~,i_SP3_linear] = min(abs(Tdata-new_ind));
            figure(1); hold on;
            subplot(3,1,1); hold on;
            delete(p51a);
            p51a=plot(Tdata(i_SP3_linear),SP3_linear,'*','color','#D95319','Markersize',ms,'Linewidth',lw2);
            delete(pleg);
            pleglines = [p11,p11a,p21,p21a,p31,p41,p51,p51a,p61a,p61];  
            pleg=legend(pleglines,plegNames,'Location','eastoutside');

        elseif strcmp(list{i},'SP4 reg') == 1
            uiwait(msgbox("Select early phase IV max SBP for regression","modal"));
            [new_ind,SP4e_Max_reg] = ginput(1);
            [~,i_SP4e_Max_reg] = min(abs(Tdata-new_ind));
            PP4e_Max_reg = PPdata(i_SP4e_Max_reg);
            figure(1); hold on;
            subplot(3,1,1); hold on;
            delete(p61a);
            p61a=plot(Tdata(i_SP4e_Max_reg),SP4e_Max_reg,'*','color','#851803','Markersize',ms,'Linewidth',lw2);
            delete(pleg);
            pleglines = [p11,p11a,p21,p21a,p31,p41,p51,p51a,p61a,p61];  
            pleg=legend(pleglines,plegNames,'Location','eastoutside');
            
        elseif strcmp(list{i},'SP4e Max') == 1
            uiwait(msgbox("Select early phase IV max SBP","modal"));
            [new_ind,SP4e_Max] = ginput(1);
            [~,i_SP4e_Max] = min(abs(Tdata-new_ind));
            PP4e_Max = PPdata(i_SP4e_Max);
            figure(1); hold on;
            subplot(3,1,1); hold on;
            delete(p61);
            p61=plot(Tdata(i_SP4e_Max),SP4e_Max,'o','color','#851803','Markersize',ms,'Linewidth',lw2);
            delete(pleg);
            pleglines = [p11,p11a,p21,p21a,p31,p41,p51,p51a,p61a,p61];  
            pleg=legend(pleglines,plegNames,'Location','eastoutside');
        
        elseif strcmp(list{i},'HR2e min') == 1
            uiwait(msgbox("Select early phase II min HR","modal"));
            [new_ind,HR2e_min] = ginput(1);
            [~,i_HR2e_min]     = min(abs(Tdata-new_ind));
            figure(1); hold on;
            subplot(3,1,2); hold on;
            delete(p12);
            p12=plot(Tdata(i_HR2e_min),HR2e_min,'co','Markersize',ms,'Linewidth',lw2);
            delete(hleg);
            hleglines = [p12,p22,p32,p42a,p42,p52,p62];
            hleg=legend(hleglines,hlegNames,'Location','eastoutside');

            subplot(3,1,3); hold on;
            delete(p13)
            RR2e_Max = 60./HR2e_min;
            p13=plot(Tdata(i_HR2e_min),RR2e_Max,'co','Markersize',ms,'Linewidth',lw2);
            delete(rleg);
            rleglines = [p13,p23,p33,p43a,p43,p53,p63];
            rleg=legend(rleglines,rlegNames,'Location','eastoutside');
             
        elseif strcmp(list{i},'HR2e end') == 1
            uiwait(msgbox("Select early phase II end HR","modal"));
            [new_ind,HR2e_end] = ginput(1);
            [~,i_HR2e_end] = min(abs(Tdata-new_ind));

            figure(1); hold on; 
           
            subplot(3,1,2); hold on;
            delete (p22);
            p22=plot(Tdata(i_HR2e_end),HR2e_end,'o','color','#0000FF','Markersize',ms,'Linewidth',lw2);
            delete(hleg);
            hleglines = [p12,p22,p32,p42a,p42,p52,p62];
            hleg=legend(hleglines,hlegNames,'Location','eastoutside');
            
            subplot(3,1,3); hold on;
            delete(p23)
            RR2e_end = 60./HR2e_end;
            p23=plot(Tdata(i_HR2e_end),RR2e_end,'o','color','#0000FF','Markersize',ms,'Linewidth',lw2);
            delete(rleg);
            rleglines = [p13,p23,p33,p43a,p43,p53,p63];
            rleg=legend(rleglines,rlegNames,'Location','eastoutside');
        
        elseif strcmp(list{i},'HR4e Max') == 1
            uiwait(msgbox("Select early phase IV max HR","modal"));
            [new_ind,HR4e_Max] = ginput(1);
            [~,i_HR4e_Max] = min(abs(Tdata-new_ind));
            figure(1); hold on;
            subplot(3,1,2); hold on;
            delete(p32);
            p32=plot(Tdata(i_HR4e_Max),HR4e_Max,'o','color','#D95319','Markersize',ms,'Linewidth',lw2);
            delete(hleg);
            hleglines = [p12,p22,p32,p42a,p42,p52,p62];
            hleg=legend(hleglines,hlegNames,'Location','eastoutside');
            
            subplot(3,1,3); hold on;
            delete(p33)
            RR4e_min = 60./HR4e_Max;
            p33=plot(Tdata(i_HR4e_Max),RR4e_min,'o','color','#D95319','Markersize',ms,'Linewidth',lw2);
            delete(rleg);
            rleglines = [p13,p23,p33,p43a,p43,p53,p63];
            rleg=legend(rleglines,rlegNames,'Location','eastoutside');

        elseif strcmp(list{i},'HR4e min reg') == 1
            uiwait(msgbox("Select early phase IV min HR for regression","modal"));
            [new_ind,HR4e_min_reg] = ginput(1);
            [~,i_HR4e_min_reg] = min(abs(Tdata-new_ind));
            figure(1); hold on;
            subplot(3,1,2); hold on;
            delete(p42a);
            p42a=plot(Tdata(i_HR4e_min_reg),HR4e_min_reg,'*','Color','#851803','Markersize',ms,'Linewidth',lw2);
            delete(hleg);
            hleglines = [p12,p22,p32,p42a,p42,p52,p62];
            hleg=legend(hleglines,hlegNames,'Location','eastoutside');
            
            subplot(3,1,3); hold on;
            delete(p43a)
            RR4e_Max = 60./HR4e_min_reg;
            p43a=plot(Tdata(i_HR4e_min_reg),RR4e_Max,'*','Color','#851803','Markersize',ms,'Linewidth',lw2);
            delete(rleg);
            rleglines = [p13,p23,p33,p43a,p43,p53,p63];
            rleg=legend(rleglines,rlegNames,'Location','eastoutside');
            
        elseif strcmp(list{i},'HR4e min') == 1
            uiwait(msgbox("Select early phase IV min HR","modal"));
            [new_ind,HR4e_min] = ginput(1);
            [~,i_HR4e_min] = min(abs(Tdata-new_ind));
            figure(1); hold on;
            subplot(3,1,2); hold on;
            delete(p42);
            p42=plot(Tdata(i_HR4e_min),HR4e_min,'o','Color','#851803','Markersize',ms,'Linewidth',lw2);
            delete(hleg);
            hleglines = [p12,p22,p32,p42a,p42,p52,p62];
            hleg=legend(hleglines,hlegNames,'Location','eastoutside');
            
            subplot(3,1,3); hold on;
            delete(p53)
            RR4e_Max = 60./HR4e_min;
            p53=plot(Tdata(i_HR4e_min),RR4e_Max,'o','Color','#851803','Markersize',ms,'Linewidth',lw2);
            delete(rleg);
            rleglines = [p13,p23,p33,p43a,p43,p53,p63];
            rleg=legend(rleglines,rlegNames,'Location','eastoutside');
            
        elseif strcmp(list{i},'HR4 min') == 1
            uiwait(msgbox("Select phase IV min HR","modal"));
            [new_ind,HR4_min] = ginput(1);
            [~,i_HR4_min] = min(abs(Tdata-new_ind));
            figure(1); hold on;
            subplot(3,1,2); hold on;
            delete(p62);
            p62=plot(Tdata(i_HR4_min),HR4_min,'o','Color','#77AC30','Markersize',ms,'Linewidth',lw2);
            delete(hleg);
            hleglines = [p12,p22,p32,p42a,p42,p52,p62];
            hleg=legend(hleglines,hlegNames,'Location','eastoutside');
            
            subplot(3,1,3); hold on;
            delete(p63);
            RR4_Max = 60./HR4_min;
            p63=plot(Tdata(i_HR4_min),RR4_Max,'o','Color','#77AC30','Markersize',ms,'Linewidth',lw2);
            delete(rleg);
            rleglines = [p13,p23,p33,p43a,p43,p53,p63];
            rleg=legend(rleglines,rlegNames,'Location','eastoutside');
        
        elseif strcmp(list{i},'HR Max') == 1
            uiwait(msgbox("Select Max HR","modal"));
            [new_ind,HR_Max] = ginput(1);
            [~,i_HR_Max] = min(abs(Tdata-new_ind));
            figure(1); hold on;
            subplot(3,1,2); hold on;
            p63=plot(Tdata(i_HR_Max),HR_Max,'mo','Markersize',ms,'Linewidth',lw2);
            delete(hleg);
            hleglines = [p12,p22,p32,p42a,p42,p52,p62,p63];
            hleg=legend(hleglines,hlegNames,'Location','eastoutside');
            
            subplot(3,1,3); hold on;
            delete(p63);
            RR_min = 60./HR_Max;
            p63=plot(Tdata(i_HR_Max),RR_min,'mo','Markersize',ms,'Linewidth',lw2);
            delete(rleg);
            rleglines = [p13,p23,p33,p43a,p43,p53,p63];
            rleg=legend(rleglines,rlegNames,'Location','eastoutside');
        end
    end
end

%% HR and RR intervals after correction
% Mean HR and RR intervals 
HRa_mean = mean(HRdata(i_HR4e_min:end)); % After VM
RRb_mean = 60./HRb_mean*1000; % Before VM
RRa_mean = 60./HRa_mean*1000; % After VM

%% Vagal baroreflex phase 2 early
% BRSv2 (HR to Tdata) overlapping segment
% x_BRSv2eHRTD  = Tdata (i_HR2e_min:i_reg2e_e);  % MO: Which one should we
% use?
% y_BRSv2eHRTD  = HRdata(i_HR2e_min:i_reg2e_e); 
x_BRSv2eHRTD    = Tdata (i_HR2e_min:i_t2e);
y_BRSv2eHRTD    = HRdata(i_HR2e_min:i_t2e); 
p_BRSv2eHRTD    = polyfit(x_BRSv2eHRTD,y_BRSv2eHRTD,1);
yfit_BRSv2eHRTD = polyval(p_BRSv2eHRTD,x_BRSv2eHRTD); 
BRSv2eHRTD      = p_BRSv2eHRTD(1);  % Slope

% R_2 calculations of BRSv2 (HR to Tdata)
resobs_BRSv2eHRTD = (y_BRSv2eHRTD - yfit_BRSv2eHRTD).^2;
SSRes_BRSv2eHRTD  = sum(resobs_BRSv2eHRTD);
ybar_BRSv2eHRTD   = mean(y_BRSv2eHRTD);
resm_BRSv2eHRTD   = (y_BRSv2eHRTD - ybar_BRSv2eHRTD).^2;
SSTot_BRSv2eHRTD  = sum(resm_BRSv2eHRTD);
R2_BRSv2eHRTD     = 1 - SSRes_BRSv2eHRTD / SSTot_BRSv2eHRTD;

% BRSv2 (RR to Tdata)
% x_BRSv2eRRTD  = Tdata (i_HR2e_min:i_reg2e_e); % MO: Which one should we
% use?
% y_BRSv2eRRTD  = RRdata(i_HR2e_min:i_reg2e_e); 
x_BRSv2eRRTD    = Tdata (i_HR2e_min:i_t2e); 
y_BRSv2eRRTD    = RRdata(i_HR2e_min:i_t2e); 
p_BRSv2eRRTD    = polyfit(x_BRSv2eRRTD,y_BRSv2eRRTD,1);
yfit_BRSv2eRRTD = polyval(p_BRSv2eRRTD,x_BRSv2eRRTD);  
BRSv2eRRTD      = p_BRSv2eRRTD(1)*1000;  % Slope

% R_2 calculations of BRSv2 (RR to Tdata)
resobs_BRSv2eRRTD = (y_BRSv2eRRTD - yfit_BRSv2eRRTD).^2;
SSRes_BRSv2eRRTD  = sum(resobs_BRSv2eRRTD);
ybar_BRSv2eRRTD   = mean(y_BRSv2eRRTD);
resm_BRSv2eRRTD   = (y_BRSv2eRRTD - ybar_BRSv2eRRTD).^2;
SSTot_BRSv2eRRTD  = sum(resm_BRSv2eRRTD);
R2_BRSv2eRRTD     = 1 - SSRes_BRSv2eRRTD / SSTot_BRSv2eRRTD;

% BRSv2 (SP to Tdata)
% x_BRSv2eSPTD  = Tdata (i_reg2e_b:i_reg2e_e);
% y_BRSv2eSPTD  = SPdata(i_reg2e_b:i_reg2e_e); 
x_BRSv2eSPTD    = Tdata (i_SP1_Max_reg:i_SP2e_linear); % CT: CHANGED
y_BRSv2eSPTD    = SPdata(i_SP1_Max_reg:i_SP2e_linear); % CT: CHANGED
p_BRSv2eSPTD    = polyfit(x_BRSv2eSPTD,y_BRSv2eSPTD,1);
yfit_BRSv2eSPTD = polyval(p_BRSv2eSPTD,x_BRSv2eSPTD);
BRSv2eSPTD      = p_BRSv2eSPTD(1);  % Slope

% R_2 calculations of BRSv2 (SP to Tdata)
resobs_BRSv2eSPTD = (y_BRSv2eSPTD - yfit_BRSv2eSPTD).^2;
SSRes_BRSv2eSPTD  = sum(resobs_BRSv2eSPTD);
ybar_BRSv2eSPTD   = mean(y_BRSv2eSPTD);
resm_BRSv2eSPTD   = (y_BRSv2eSPTD - ybar_BRSv2eSPTD).^2;
SSTot_BRSv2eSPTD  = sum(resm_BRSv2eSPTD);
R2_BRSv2eSPTD     = 1 - SSRes_BRSv2eSPTD / SSTot_BRSv2eSPTD;

% Overlapping segments early phase 2
i_reg2e_b = max(i_SP1_Max_reg,i_HR2e_min); % Overlapping segment start
i_reg2e_e = min(i_SP2e_linear,i_t2e);      % Overlapping segment end

% Time BP drop to HR increase
dtoff_2e = Tdata(i_HR2e_min)-Tdata(i_SP1_Max);

% BRSv (RR to SP)
x_BRSv2eRRSP    = SPdata(i_reg2e_b:i_reg2e_e);
y_BRSv2eRRSP    = RRdata(i_reg2e_b:i_reg2e_e); 
p_BRSv2eRRSP    = polyfit(x_BRSv2eRRSP,y_BRSv2eRRSP,1);
yfit_BRSv2eRRSP = polyval(p_BRSv2eRRSP,x_BRSv2eRRSP);
BRSv2eRRSP      = p_BRSv2eRRSP(1)*1000; % Slope note converted to ms

% R_2 calculations of BRSv (RR to SP)
resobs_BRSv2eRRSP = (y_BRSv2eRRSP - yfit_BRSv2eRRSP).^2;
SSRes_BRSv2eRRSP  = sum(resobs_BRSv2eRRSP);
ybar_BRSv2eRRSP   = mean(y_BRSv2eRRSP);
resm_BRSv2eRRSP   = (y_BRSv2eRRSP - ybar_BRSv2eRRSP).^2;
SSTot_BRSv2eRRSP  = sum(resm_BRSv2eRRSP);
R2_BRSv2eRRSP     = 1 - SSRes_BRSv2eRRSP / SSTot_BRSv2eRRSP;

% BRSv (HR to SP)
x_BRSv2eHRSP    = SPdata(i_reg2e_b:i_reg2e_e); 
y_BRSv2eHRSP    = HRdata(i_reg2e_b:i_reg2e_e); 
p_BRSv2eHRSP    = polyfit(x_BRSv2eHRSP,y_BRSv2eHRSP,1);
yfit_BRSv2eHRSP = polyval(p_BRSv2eHRSP,x_BRSv2eHRSP);
BRSv2eHRSP      = p_BRSv2eHRSP(1);  % Slope

% R_2 calculations of BRSv (HR to SP)
resobs_BRSv2eHRSP = (y_BRSv2eHRSP - yfit_BRSv2eHRSP).^2;
SSRes_BRSv2eHRSP  = sum(resobs_BRSv2eHRSP);
ybar_BRSv2eHRSP   = mean(y_BRSv2eHRSP);
resm_BRSv2eHRSP   = (y_BRSv2eHRSP - ybar_BRSv2eHRSP).^2;
SSTot_BRSv2eHRSP  = sum(resm_BRSv2eHRSP);
R2_BRSv2eHRSP     = 1 - SSRes_BRSv2eHRSP / SSTot_BRSv2eHRSP;

%% Vagal baroreflex phase 4 early
% HR vs Tdata
x_BRSv4eHRTD    = Tdata (i_HR4e_Max:i_HR4e_min_reg); 
y_BRSv4eHRTD    = HRdata(i_HR4e_Max:i_HR4e_min_reg);
p_BRSv4eHRTD    = polyfit(x_BRSv4eHRTD,y_BRSv4eHRTD,1);
yfit_BRSv4eHRTD = polyval(p_BRSv4eHRTD,x_BRSv4eHRTD);
BRSv4eHRTD      = p_BRSv4eHRTD(1); % Slope 

% R_2 calculations of HR vs Tdata
resobs_BRSv4eHRTD = (y_BRSv4eHRTD - yfit_BRSv4eHRTD).^2;
SSRes_BRSv4eHRTD  = sum(resobs_BRSv4eHRTD);
ybar_BRSv4eHRTD   = mean(y_BRSv4eHRTD);
resm_BRSv4eHRTD   = (y_BRSv4eHRTD - ybar_BRSv4eHRTD).^2;
SSTot_BRSv4eHRTD  = sum(resm_BRSv4eHRTD);
R2_BRSv4eHRTD     = 1 - SSRes_BRSv4eHRTD / SSTot_BRSv4eHRTD;

% RR vs Tdata
x_BRSv4eRRTD    = Tdata (i_HR4e_Max:i_HR4e_min_reg); 
y_BRSv4eRRTD    = RRdata(i_HR4e_Max:i_HR4e_min_reg); 
p_BRSv4eRRTD    = polyfit(x_BRSv4eRRTD,y_BRSv4eRRTD,1);
yfit_BRSv4eRRTD = polyval(p_BRSv4eRRTD,x_BRSv4eRRTD);
BRSv4eRRTD      = p_BRSv4eRRTD(1)*1000; % Slope 

% R_2 calculations of RR vs Tdata
resobs_BRSv4eRRTD = (y_BRSv4eRRTD - yfit_BRSv4eRRTD).^2;
SSRes_BRSv4eRRTD  = sum(resobs_BRSv4eRRTD);
ybar_BRSv4eRRTD   = mean(y_BRSv4eRRTD);
resm_BRSv4eRRTD   = (y_BRSv4eRRTD - ybar_BRSv4eRRTD).^2;
SSTot_BRSv4eRRTD  = sum(resm_BRSv4eRRTD);
R2_BRSv4eRRTD     = 1 - SSRes_BRSv4eRRTD / SSTot_BRSv4eRRTD;

% BP vs Tdata
x_BRSv4eSPTD    = Tdata(i_t3linear:i_SP4e_Max_reg); 
y_BRSv4eSPTD    = SPdata(i_t3linear:i_SP4e_Max_reg);
p_BRSv4eSPTD    = polyfit(x_BRSv4eSPTD,y_BRSv4eSPTD,1);
yfit_BRSv4eSPTD = polyval(p_BRSv4eSPTD,x_BRSv4eSPTD);
BRSv4eSPTD      = p_BRSv4eSPTD(1); % Slope 

% R_2 calculations of BP vs Tdata
resobs_BRSv4eSPTD = (y_BRSv4eSPTD - yfit_BRSv4eSPTD).^2;
SSRes_BRSv4eSPTD  = sum(resobs_BRSv4eSPTD);
ybar_BRSv4eSPTD   = mean(y_BRSv4eSPTD);
resm_BRSv4eSPTD   = (y_BRSv4eSPTD - ybar_BRSv4eSPTD).^2;
SSTot_BRSv4eSPTD  = sum(resm_BRSv4eSPTD);
R2_BRSv4eSPTD     = 1 - SSRes_BRSv4eSPTD / SSTot_BRSv4eSPTD;

% BP vs RR
dtoff_4e = Tdata(i_HR4e_Max)-Tdata(i_t3); % Time from phase 3 end to max heart rate early phase 4

ix4 = max(i_t3linear,i_HR4e_Max);
iy4 = min(i_SP4e_Max_reg,i_HR4e_min_reg);
x_BRSv4eRRSP    = SPdata(ix4:iy4);
y_BRSv4eRRSP    = RRdata(ix4:iy4); 
p_BRSv4eRRSP    = polyfit(x_BRSv4eRRSP,y_BRSv4eRRSP,1);
yfit_BRSv4eRRSP = polyval(p_BRSv4eRRSP,x_BRSv4eRRSP);
BRSv4eRRSP      = p_BRSv4eRRSP(1)*1000; % Slope

% R_2 calculations of BP vs RR
resobs_BRSv4eRRSP = (y_BRSv4eRRSP - yfit_BRSv4eRRSP).^2;
SSRes_BRSv4eRRSP  = sum(resobs_BRSv4eRRSP);
ybar_BRSv4eRRSP   = mean(y_BRSv4eRRSP);
resm_BRSv4eRRSP   = (y_BRSv4eRRSP - ybar_BRSv4eRRSP).^2;
SSTot_BRSv4eRRSP  = sum(resm_BRSv4eRRSP);
R2_BRSv4eRRSP     = 1 - SSRes_BRSv4eRRSP / SSTot_BRSv4eRRSP;

% BP vs HR
x_BRSv4eHRSP    = SPdata(ix4:iy4); 
y_BRSv4eHRSP    = HRdata(ix4:iy4); 
p_BRSv4eHRSP    = polyfit(x_BRSv4eHRSP,y_BRSv4eHRSP,1);
yfit_BRSv4eHRSP = polyval(p_BRSv4eHRSP,x_BRSv4eHRSP);
BRSv4eHRSP      = p_BRSv4eHRSP(1); % Slope
 
% R_2 calculations of BP vs HR
resobs_BRSv4eHRSP = (y_BRSv4eHRSP - yfit_BRSv4eHRSP).^2;
SSRes_BRSv4eHRSP  = sum(resobs_BRSv4eHRSP);
ybar_BRSv4eHRSP   = mean(y_BRSv4eHRSP);
resm_BRSv4eHRSP   = (y_BRSv4eHRSP - ybar_BRSv4eHRSP).^2;
SSTot_BRSv4eHRSP  = sum(resm_BRSv4eHRSP);
R2_BRSv4eHRSP     = 1 - SSRes_BRSv4eHRSP / SSTot_BRSv4eHRSP;

%% Adrenergic baroreflex late phase 2
x_alpha_BRSa    = Tdata(i_t2e:i_t2l_Max); 
y_alpha_BRSa    = SPdata(i_t2e:i_t2l_Max);
p_alpha_BRSa    = polyfit(x_alpha_BRSa,y_alpha_BRSa,1);
yfit_alpha_BRSa = polyval(p_alpha_BRSa,x_alpha_BRSa); 
BRSa2lSPTD      = p_alpha_BRSa(1); % Slope

% R_2 calculations
resobs_alpha_BRSa = (y_alpha_BRSa - yfit_alpha_BRSa).^2;
SSRes_alpha_BRSa  = sum(resobs_alpha_BRSa);
ybar_alpha_BRSa   = mean(y_alpha_BRSa);
resm_alpha_BRSa   = (y_alpha_BRSa - ybar_alpha_BRSa).^2;
SSTot_alpha_BRSa  = sum(resm_alpha_BRSa);
R2_BRSa2lSPTD     = 1 - SSRes_alpha_BRSa / SSTot_alpha_BRSa;

% Time BP drop to HR increase
dtoff_2l = Tdata(i_t2l)-Tdata(i_t2l_Max);

%% Palamarchuk et al. 
% Blood pressure differences
A = SP2e_min - SP_mean;      %  Drop in phase II from baseline (neg)
B = SP3_end  - SP2l_end;     %  Drop over phase III (neg)
C = SP2l_end - SP2e_min;     %  Change in late phase II
D = SP_mean  - SP3_end;      %  Change in phase III from baseline (pos)
E = SP4e_Max - SP_mean;      %  Overshoot in phase IV (pos)

% Triangle late phase II
adj_a  = x_alpha_BRSa(end)-x_alpha_BRSa(1);
opp_a  = yfit_alpha_BRSa(end)-yfit_alpha_BRSa(1);
hyp_a  = sqrt(opp_a^2 + adj_a^2);

% Adrenergic baroreflex sensitivity
BRSa  = abs(B) / T_PRT; 
BRSa1 = abs(A+0.75*B)/T_PRT; 

% Alpha adrenergic 
alpha_BRSa     = C/T_2l; % equivalent to regression line

% Angle alpha
alpha          = atand(opp_a/adj_a);
alpha_Area     = adj_a*opp_a/2;
alpha_Area2    = T_2l*C/2;  

% Traingle early phse IV
x_beta_BRSa    = Tdata(i_t3:i_PRT);
y_beta_BRSa    = SPdata(i_t3:i_PRT);
p_beta_BRSa    = polyfit(x_beta_BRSa,y_beta_BRSa,1);
yfit_beta_BRSa = polyval(p_beta_BRSa,x_beta_BRSa); 

adj_b          = x_beta_BRSa(end) - x_beta_BRSa(1);
opp_b          = yfit_beta_BRSa(end)-yfit_beta_BRSa(1);
hyp_b          = sqrt(opp_b^2 + adj_b^2);

% Beta adrenergic 
cosbeta_BRSb    = (T_PRT * D)/hyp_b; 
beta_BRSa       = D/T_PRT; % equivalent to regression line

beta          = atand(opp_b/adj_b); 
beta_Area     = adj_b*opp_b/2;  
 beta_Area2   = D*T_PRT/2;

% Alpha and Beta
BRSa_Area  = alpha_Area*beta_Area;
BRSa_Area2 = alpha_Area2*beta_Area2;

% Global BRS 
BRSg     = BRSa *(BRSv2eRRSP/1000); % BRSa (mmHg/s)  * (s/mmHg)
BRSg1    = BRSa1*(BRSv2eRRSP/1000); % BTSa1 (mmHg/s) * (s/mmHg)

% Valsalva ratio
VR1 = HR4e_Max/HR4_min;
VR2 = HR_Max/HR4_min;
 
%% Figure properties
fs   = plotmarkers.fs;
fsL  = fs + 3;
ms   = plotmarkers.ms;
lwt  = plotmarkers.lwt;

% Time intervals
xt(1) = mean(Tdata(i_ts:i_t1));
xt(2) = mean(Tdata(i_t1:i_t2l));
xt(3) = Tdata(i_t2l+2);
xt(4) = mean(Tdata(i_t3:i_t4));

% Axis limits
HRmin = floor(min(HRdata)/5)*5-2;
HRmax = ceil(max(HRdata)/5)*5+2;
SPmin = floor(min(data.Praw)/10)*10-2;
SPmax = ceil(max(SPdata)/10)*10+10;
RRmin = floor(min(RRdata)/10)*10-2;
RRmax = ceil(max(RRdata)/10)*10+10;

f = figure(1); hold on;
subplot(3,1,1); hold on
   p71 = plot(x_BRSv2eSPTD,yfit_BRSv2eSPTD,'linewidth',lw,'color','#00FFFF') ;
   p81 = plot(x_alpha_BRSa,yfit_alpha_BRSa,'linewidth',lw,'color','#0000FF');
   p101= plot(x_BRSv4eSPTD,yfit_BRSv4eSPTD,'linewidth',lw,'color','#D95319');
   delete(pleg);
   pleglines = [p11,p11a,p21,p21a,p31,p41,p51,p51a,p61a,p61,p71,p81,p101];  
   pleg=legend(pleglines,plegNames,'Location','eastoutside');
subplot(3,1,2); hold on
   p72 = plot(x_BRSv2eHRTD, yfit_BRSv2eHRTD, 'linewidth',lw,'color','#00FFFF');
   p82 = plot(x_BRSv4eHRTD, yfit_BRSv4eHRTD, 'linewidth',lw,'color','#D95319');
   delete(hleg);
   hleglines = [p12,p22,p32,p42a,p42,p52,p62,p72,p82];
   hleg=legend(hleglines,[hlegNames,'reg 2e','reg 4e'],'Location','eastoutside');
            
subplot(3,1,3); hold on
   p73 = plot(x_BRSv2eRRTD, yfit_BRSv2eRRTD,'linewidth',lw,'color','#00FFFF');
   p83 = plot(x_BRSv4eRRTD, yfit_BRSv4eRRTD,'linewidth',lw,'color','#D95319');
   delete(rleg);
   rleglines = [p13,p23,p33,p43a,p43,p53,p63,p73,p83];
   rleg=legend(rleglines,[rlegNames,'reg 2e','reg 4e'],'Location','eastoutside');
   
  
   resobs = (y_BRSv2eRRSP - yfit_BRSv2eRRSP).^2;
   SSRes  = sum(resobs);
   ybar   = mean(y_BRSv2eRRSP);
   resm   = (y_BRSv2eRRSP-ybar).^2;
   SSTot  = sum(resm);
   Rsqr1  = 1 - SSRes/SSTot;
   
   figFolder = plotmarkers.figFolder;
   fileName = strcat(patient,'_ratios.png');
   if nargin == 3
       fullFileName = fullfile(figFolder,'Data/',patids,fileName);
   else
       fullFileName = fullfile(figFolder,'Data',fileName);
   end
   exportgraphics(f,fullFileName);

   bbutton = uicontrol('Parent',f,'Style','pushbutton',...
      'Position',[25,10,figSize(3)*0.25,figSize(4)*0.05],...
      'String','Save and exit','fontsize',fs,'Callback',@closeGenButton);
  
   uiwait(f)
   
   screenSize     = get(0,'screensize');
   figSizeF2      = screenSize;
   figSizeF2(1:2) = 0.12*screenSize(3:4);
   figSizeF2(3)   = 0.4*screenSize(3);
   figSizeF2(4)   = 0.6*screenSize(4);
   formatSpec = '%.2f'; 

   f2 = figure(2);
   set(gcf,'units','points','position',figSizeF2)
   f2.Units = 'pixels'; 
   sgtitle(str2,'Fontsize',fs,'FontWeight','bold','interpreter','none')
       
   subplot(1,2,1);
       h = plot(x_BRSv2eRRSP,y_BRSv2eRRSP,'co',x_BRSv2eRRSP,yfit_BRSv2eRRSP,'c','linewidth',lw);
       set(gca,'fontsize',fs);
       title(strcat('Early Phase II, R^2 = ', num2str(Rsqr1,3)));
       xlabel('SBP (mmHg)');
       ylabel('RR (s)');

       resobs = (y_BRSv4eRRSP - yfit_BRSv4eRRSP).^2;
       SSRes  = sum(resobs);
       ybar   = mean(y_BRSv4eRRSP);
       resm   = (y_BRSv4eRRSP-ybar).^2;
       SSTot  = sum(resm);
       Rsqr2  = 1 - SSRes/SSTot;
       
    subplot(1,2,2);
       h = plot(x_BRSv4eRRSP,y_BRSv4eRRSP,'o',x_BRSv4eRRSP,yfit_BRSv4eRRSP,'color','#D95319','linewidth',lw);
       set(gca,'fontsize',fs);
       title(strcat('Early Phase IV, R^2 = ', num2str(Rsqr2,3)));
       xlabel('SBP (mmHg)');
       
   fileName = strcat(patient,'_ratios_regresion.png');
   if nargin == 3
       fullFileName = fullfile(figFolder,'Data/',patids,fileName);
   else
       fullFileName = fullfile(figFolder,'Data',fileName);
   end
   exportgraphics(f2,fullFileName);
  
   bbutton = uicontrol('Parent',f2,'Style','pushbutton',...
      'Position',[25,2,figSize(3)*0.2,figSize(4)*0.05],...
      'String','Save and exit','fontsize',fs,'Callback',@closeGenButton);
  
   uiwait(f2)

%% Outputs 
% Markers added to data
data.i_HR4e_min = i_HR4e_min;
data.Hbara      = HRa_mean;
data.i_HR4e_Max = i_HR4e_Max;        

% Change units to milliseconds
RR2e_Max = RR2e_Max*1000;
RR2e_end = RR2e_end*1000;
RR4e_min = RR4e_min*1000;
RR4e_Max = RR4e_Max*1000;
RR4_Max  = RR4_Max*1000;

name = patient;

ClinicalRatios = table({name},T_1, T_2e, T_2l, T_3, T_PRT, ... 
    SP_mean, SP1_Max, SP2e_min, SP2l_Max, SP2l_end, SP3_end, SP4e_Max, ...
    PPbar_b, PP1_Max, PP2e_min, PP4_Max, PPbar_a, ....
    HRb_mean, HRa_mean, HR2e_min, HR2e_end, HR4e_Max, HR4e_min, HR4_min, HR_Max, VR1, VR2, ...
    RRb_mean, RRa_mean, RR2e_Max, RR2e_end, RR4e_min, RR4e_Max, RR4_Max, ...
    BRSv2eHRTD, R2_BRSv2eHRTD, BRSv2eRRTD, R2_BRSv2eRRTD, ...
    BRSv2eSPTD, R2_BRSv2eSPTD, BRSv2eHRSP, R2_BRSv2eHRSP, dtoff_2e, ...
    BRSv2eRRSP, R2_BRSv2eRRSP, BRSa2lSPTD, R2_BRSa2lSPTD, dtoff_2l, ...
    BRSv4eHRTD, R2_BRSv4eHRTD, BRSv4eRRTD, R2_BRSv4eRRTD, ...
    BRSv4eSPTD, R2_BRSv4eSPTD, BRSv4eHRSP, R2_BRSv4eHRSP, ...
    BRSv4eRRSP, R2_BRSv4eRRSP, dtoff_4e, ...
    A, B, C, D, E, ...
    BRSa, BRSa1, ...
    alpha_BRSa, beta_BRSa, ...
    alpha, beta,...
    alpha_Area, alpha_Area2, beta_Area, beta_Area2, BRSa_Area, BRSa_Area2, ...
    BRSg, BRSg1, dVPP, dVPP1, dVSP); 

 % Save results for each patient in markers
 if nargin == 3
     writetable(ClinicalRatios,strcat('../Markers/',patids,patient,'_markers.xlsx'))
 else
     writetable(ClinicalRatios,strcat('../Markers/',patient,'_markers.xlsx'))
 end
 
 % Appending workspace
 
 if nargin == 3
     save(strcat('../WS/',patids,patient,'_WS.mat'),'data','pat_char','raw_data');
 else
     save(strcat('../WS/',patient,'_WS.mat'),'data','pat_char','raw_data');
 end
 
end % function clinicalratios %