function create_WS_indices(patient,plotmarkers,patids)
% function create_WS_indices
% Input: patient name and figure settings
% Output: augmented data strucure and temporary workspace (saved in WS)
% Uses: set_indices to determine VM phases
% Description: generates valsalva maneuver phases 
% create_WS_BP loads BP and raw data file and creates workspace for patient

% global Height Sex Weight SPinds DPraw DPinds dt ECGraw ...
%       Hraw HRdata PPraw Praw Rdata RRdata SPraw Traw 

if nargin == 3
    load(strcat('../WS/',patids,patient,'_WS_temp.mat'));
else
    load(strcat('../WS/',patient,'_WS_temp.mat'));
end

PatientID   = pat_char.PatientID;
Sex         = pat_char.sex;
Height      = pat_char.height;
Weight      = pat_char.weight;

Traw        = raw_data.Traw-raw_data.Traw(1);
ECGraw      = raw_data.ECGraw;
Praw        = raw_data.Praw;
if isfield(raw_data,'dt')
    dt      = raw_data.dt;
else
    dt      = data.dt;
end

PthAvail    = raw_data.PthAvail;
if PthAvail ~= 0
    Pthraw  = raw_data.Pthraw;
end

HRraw      = raw_data.HRraw;
RRraw      = data.RRraw;
Rraw       = data.Rraw;

SPinds      = data.SPinds;
DPinds      = data.SPinds;

SPraw       = data.SPraw;
DPraw       = data.DPraw;
PPraw       = SPraw - DPraw;

%% Figure properties
fs = plotmarkers.fs;
lw = plotmarkers.lwt;
if lw > 2
    lw2 = lw -1;
else
    lw2 = 1;
end
figSize = plotmarkers.figSize;

SPbar = trapz(SPraw(1:1000-1))/(1000 - 1); % was 100 very short to take the average over

%% If There Is No Pth data Available, Find Start and Stop
if PthAvail == 0
    
    f = figure(1); clf; hold on;
    set(gcf,'units','points','position',figSize)
    f.Units = 'pixels';
    sgtitle(patient,'Fontsize',fs,'FontWeight','bold','Interpreter', 'none')

    subplot(4,1,1); hold on
      plot(Traw,Praw,'Color',[0 0.4470 0.7410]);
      plot(Traw,SPraw,'b','LineWidth',lw2);
      set(gca,'Fontsize',fs);
      title('VM start (min before first peak)','VM end (max before drop & overshoot)');
      ylabel('BP (mmHg)')
      legend('BP','SBP')
      xlim([Traw(1) Traw(end)]);
      ylim([min(Praw)-5 max(Praw)+5])

    subplot(4,1,2); hold on
      plot(Traw,HRraw,'b','LineWidth',lw2,'Displayname','HR (bpm)');
      set(gca,'Fontsize',fs);
      ylabel('HR (bpm)');
      xlim([Traw(1) Traw(end)]);
      ylim([min(HRraw)-2 max(HRraw)+2])

    subplot(4,1,3); hold on
      plot(Traw,Rraw,'b','LineWidth',lw2);
      set(gca,'Fontsize',fs);
      ylabel('Resp (mV)');
      xlim([Traw(1) Traw(end)]);

    subplot(4,1,4); hold on
      plot(Traw,ECGraw,'b');
      set(gca,'Fontsize',fs);
      ylabel('ECG (mV)');
      xlabel('Time (s)')
      xlim([Traw(1) Traw(end)]);
    
    [val_times, ~ ] = ginput(2);
    val_start       = val_times(1);
    val_end         = val_times(2);
    [~,i_ts_raw]    = min(abs(Traw - val_start));
    [~,i_te_raw]    = min(abs(Traw - val_end));
    
    if i_ts_raw >= i_te_raw || isempty(i_te_raw) || isempty(i_ts_raw)
        figure(1); hold on
        title('Error -- click start/stop again')
        [val_times, ~ ] = ginput(2);
        val_start = val_times(1);
        val_end   = val_times(2);
        [~,i_ts_raw] = min(abs(Traw - val_start)); % Find the indices of
        % the time data points that are closest to the VM start
        [~,i_te_raw] = min(abs(Traw - val_end)); % Find the indices of
        % the time data points that are closest to the VM end
    end
else
    ID1 = find(Pthraw>10,1);
    ID2 = find(Pthraw>30,1,'last');
    val_start = Traw(ID1);
    val_end   = Traw(ID2);
    
    disp('Pth available')
    i_ts_raw = ID1; 
    i_te_raw = ID2;
end
% if Pth is avail make sure to generate i_ts_raw, i_te_raw, val_start and val_end
% end of the if statement
% if there is pth data then take out bottom graph and insert pth graph at
% the top and sense the points and give an indicator line of where those
% points are (val_start and val_end) then have a pop up box for corrections

%% Calculate times to cut data before and after VM 

% Calculate available time at start and end of maneuver
timeAvailableS = Traw(i_ts_raw) - Traw(1);   % Breath hold start
timeAvailableE = Traw(end) - Traw(i_te_raw); % Breath hold end

% Rest time before and after (25 secs)
restTimeS = 25; 
restTimeE = 25;
 
% Test if the rest time asked for exceeds available time
 dtS = 0;
 dtE = 0;
 if timeAvailableS > restTimeS
     dtS = timeAvailableS-restTimeS;
     disp(strcat(num2str(PatientID),{' '},num2str(dtS),{' '},'rest time before too long'));
end
if timeAvailableE > restTimeE
     dtE = timeAvailableE-restTimeE;
     disp(strcat(num2str(PatientID),{' '},num2str(dtE),{' '},'rest time after too long'));
end

% Cut raw data
icutS  = find(Traw >= Traw(1)   + dtS,1);
icutE  = find(Traw >= Traw(end) - dtE,1); 

Traw   = Traw(icutS:icutE)-Traw(icutS);
Praw   = Praw(icutS:icutE);
ECGraw = ECGraw(icutS:icutE);
SPraw  = SPraw(icutS:icutE);
DPraw  = DPraw(icutS:icutE);
PPraw  = PPraw(icutS:icutE);
Rraw   = Rraw(icutS:icutE);
HRraw  = HRraw(icutS:icutE);
RRraw  = RRraw(icutS:icutE);
if PthAvail ~= 0
     Pthraw  = Pthraw(icutS:icutE);
end

% Subsampling at 100 Hz
v = (1:100:length(Traw))'; 
Tdata       = Traw(v);
HRdata      = HRraw(v);
Rdata       = Rraw(v);
RRdata      = RRraw(v);
SPdata      = SPraw(v);
DPdata      = DPraw(v);
PPdata      = SPdata-DPdata;
PPdata      = PPraw(v);
if PthAvail ~= 0
     Pthdata     = Pthraw(v);
end

% set indices structure to pass in
set_ind.Traw        = Traw;
set_ind.Praw        = Praw;
set_ind.Rdata       = Rdata;
set_ind.Tdata       = Tdata;
set_ind.SPdata      = SPdata;
set_ind.DPdata      = DPdata;
set_ind.PPdata      = PPdata;
set_ind.HRdata      = HRdata;
set_ind.RRdata      = RRdata;
set_ind.val_start   = val_start-dtS;
set_ind.val_end     = val_end-dtS;
if PthAvail ~= 0
    set_ind.Pthdata = Pthdata; % Pthdata (subsampled)
end
set_ind.PthAvail = PthAvail;

% Set Indices should not be subsampled(cut Pth above)
[phases]   = set_indices(set_ind, patient, plotmarkers);
i_ts       = phases(1); % breath hold start
i_t1linear = phases(2);
i_t1       = phases(3); % end of phase I
i_t2elinear= phases(4); % end of early phase II
i_t2emin   = phases(5); % minimum point of phase II
i_t2e      = phases(6);
i_t2l      = phases(7); % end of late phase II
i_t3       = phases(8); % end of phase III
i_t3linear = phases(9);
i_PRT      = phases(10); % end of pressure recovery time 
i_t4       = phases(11); % end of phase IV
i_t4linear = phases(12);
if i_t4linear > length(Tdata)
    disp('i_t4linear too large');
end
% Baseline Values
i_ts2 = round(i_ts/2);
SPbar = trapz(SPdata(1:i_ts - 1))/(i_ts - 1);
Hbar  = trapz(HRdata(1:i_ts2))/(i_ts2);
Rbar  = trapz(Rdata(1:i_ts - 1))/(i_ts - 1);
PPbar = trapz(PPdata(1:i_ts - 1))/(i_ts -1);

% Calculate mean distance between timepoints
dt = mean(diff(Tdata));

% Find max and min HR during respiration
HRdatarespmax = max(HRdata(1:i_ts));
HRdatarespmin = min(HRdata(1:i_ts));
HRdatarespamp = HRdatarespmax - HRdatarespmin;

[HRdatapeaks,~] = findpeaks(HRdata(1:i_ts),'MinPeakProminence',0.25 * HRdatarespamp);
[~,v_locs_data] = findpeaks(-HRdata(1:i_ts),'MinPeakProminence',0.25 * HRdatarespamp);
HRdatavalleys   = HRdata(v_locs_data);

if isempty(HRdatapeaks)
    [HRdatapeaks,~] = max(HRdata(1:i_ts-round(5/dt)));
end
if isempty(HRdatavalleys)
    [HRdatavalleys,~] = min(HRdata(1:i_ts-round(5/dt)));
end

HmaxR = mean(HRdatapeaks);
HminR = mean(HRdatavalleys);

% Make thoracic pressure signal that is strictly positive
Rnew  = 4*(Rdata - mean(Rdata(1:i_ts-1)))/(max(Rdata(1:i_ts-1))-min(Rdata(1:i_ts-1))); % 0 mean value
rhogh = 1060*0.2*9.81/133.322 - 6;
PR    = Rnew + rhogh;

 if PthAvail == 0 %%%% if there is no Pth data
    %Thoracic pressure
    Pth = zeros(size(Tdata));
    for i = 1:length(Tdata)
        t = Tdata(i);
        if (t > Tdata(i_ts)) && (t < Tdata(i_t2l))
            p = 40*(1 - exp(-2*(t - Tdata(i_ts))));
        else
            p = PR(i);
        end
        Pth(i) = p;
    end
 else
    Pth = zeros(size(Tdata));
    for i = 1:length(Tdata)
        t = Tdata(i);
        if (t > Tdata(i_ts)) && (t < Tdata(i_t2l))
            p = Pthdata(i);
        else
            p = PR(i);
        end
        Pth(i) = p;
    end
 end
 Pthdata = movmean(Pth,10); 
 Pthbar  = trapz(Pthdata(1:i_ts-1))/(i_ts-1); 
 
%% Find derivatives of pressure data    
% Derivative of SBP on subsampled data
dSBPdt = diff(SPdata)/dt; % Forward difference
dSBPdt = [dSBPdt;dSBPdt(end)];

% Derivative of PP on subsampled data
dPPdt = diff(PPdata)/dt;
dPPdt = [dPPdt;dPPdt(end)];

% Mean values before VM
dSBPdt_bar = trapz(dSBPdt(1:i_ts - 1))/(i_ts - 1);
dPPdt_bar  = trapz(dPPdt(1:i_ts - 1))/(i_ts - 1);

% Derivative of Pth
dPthdt = diff(Pthdata)/dt;
dPthdt = [dPthdt;dPthdt(end)];
dPthdt_bar = trapz(dPthdt(1:i_ts - 1))/(i_ts - 1);

% Make Data Structure 
data.Traw       = Traw;    
data.HRraw      = HRdata;  
data.PPraw      = PPraw;   
data.Praw       = Praw;
data.ECGraw     = ECGraw;
data.RRdata     = RRdata;
data.Tdata      = Tdata;
data.PPdata     = PPdata;
data.SPdata     = SPdata;
data.HRdata     = HRdata;
data.Pthdata    = Pthdata;
data.Rdata      = Rdata;
data.DPdata     = DPdata;
data.RRdata     = RRdata;

data.val_start  = val_start;
data.val_end    = val_end;

data.i_ts        = i_ts;
data.i_t1linear  = i_t1linear;
data.i_t1        = i_t1;
data.i_t2elinear = i_t2elinear;
data.i_t2emin    = i_t2emin;
data.i_t2e       = i_t2e;
data.i_t2l       = i_t2l;
data.i_t3        = i_t3;
data.i_t3linear  = i_t3linear;
data.i_PRT       = i_PRT;
data.i_t4        = i_t4;
data.i_t4linear  = i_t4linear;

data.HminR      = HminR;
data.HmaxR      = HmaxR;
data.Pthbar     = Pthbar;
data.SPbar      = SPbar;
data.PPbar      = PPbar;
data.Hbar       = Hbar;
data.Rbar       = Rbar;

data.dt         = dt;

data.dSBPdt     = dSBPdt;
data.dPPdt      = dPPdt;
data.dPthdt     = dPthdt;
data.dSBPdt_bar = dSBPdt_bar;
data.dPPdt_bar  = dPPdt_bar;
data.dPthdt_bar = dPthdt_bar;

% Four Panel Figure
% X- and Y- axis limits for plots
Tlims   = [Tdata(1) Tdata(end)];
Plims   = [min(Praw)-5 max(Praw)+2];
Hlims   = [min(HRdata)-2 max(HRdata)+2];
Pthlims = [min(Pthdata)-2 max(Pthdata)+2];
ECGlims = 1000*[min(data.ECGraw)*1.1, max(data.ECGraw)*1.1];
dlims   = [-60 60]; 

% Times for VM phases
tsVM  = Tdata(data.i_ts);
t1VM  = Tdata(data.i_t1);
t2eVM = Tdata(data.i_t2e);
t2lVM = Tdata(data.i_t2l);
t3VM  = Tdata(data.i_t3);
t4VM  = Tdata(data.i_t4);

if data.i_t4 <= length(Tdata)
    t4VM = Tdata(data.i_t4);
end

% Identity
I = ones(2,1);

% X values for shaded regions for each VM phase
x1 = [tsVM t1VM t1VM tsVM];
x2 = [t1VM t2lVM t2lVM t1VM];
x3 = [t2lVM t3VM t3VM t2lVM];
x4 = [t3VM t4VM t4VM t3VM];

% Y values for shaded regions for each VM phase
yH   = [Hlims(1) Hlims(1) Hlims(2) Hlims(2)];
yP   = [Plims(1) Plims(1) Plims(2) Plims(2)];
yPth = [Pthlims(1) Pthlims(1) Pthlims(2) Pthlims(2)];
yd   = [dlims(1) dlims(1) dlims(2) dlims(2)];
yECG = [ECGlims(1) ECGlims(1) ECGlims(2) ECGlims(2)];

% Colors
gray  = [.9 .9 .9];
lgray = [.95 .95 .95];

f = figure(3); clf;
set(gcf,'units','points','position',figSize)
f.Units = 'pixels';
sgtitle(patient,'Fontsize',fs,'FontWeight','bold','Interpreter', 'none')

subplot(3,1,1); hold on;
  patch(x1,yP,gray)
  patch(x2,yP,lgray)
  patch(x3,yP,gray)
  patch(x4,yP,lgray)
  plot(t2eVM * I,Plims,'k:','LineWidth',1)
  plot(Tlims, SPbar * I,'k:','LineWidth',1)
  plot(Tlims, PPbar * I,'k:','linewidth',1)
  h1 = plot(data.Traw,data.Praw,'Color',[0 0.4470 0.7410]);
  h2 = plot(Tdata,SPdata,'b','linewidth',lw2);
  set(gca,'Fontsize',fs); 
  legend([h1 h2],'BP','SBP','location','southwest')
  ylabel('BP (mmHg)')
  xlim(Tlims);
  ylim(Plims);

subplot(3,1,2); hold on;
  patch(x1,yH,gray)
  patch(x2,yH,lgray)
  patch(x3,yH,gray)
  patch(x4,yH,lgray)
  plot(t2eVM * I,Hlims,'k:','LineWidth',1)
  plot(Tlims, Hbar * I,'k:','LineWidth',1)
  h1 = plot(Tdata, HRdata,'b','LineWidth',lw);
  set(gca,'Fontsize',fs)
  ylabel('HR (bpm)')
  xlim(Tlims)
  ylim(Hlims)

subplot(3,1,3); hold on; 
  patch(x1,yPth,gray)
  patch(x2,yPth,lgray)
  patch(x3,yPth,gray)
  patch(x4,yPth,lgray)
  plot(t2eVM * I,Pthlims,'k:','LineWidth',1)
  plot(Tlims, Pthbar*ones(2,1),'k:','linewidth',1)
  plot(Tdata, Pthdata,'Color','b','LineWidth',lw)
  set(gca,'Fontsize',fs); 
  xlabel('Time (s)');
  ylabel('Pth (mmHg)');
  xlim(Tlims);
  ylim(Pthlims);
  figFolder = plotmarkers.figFolder;
  fileName = strcat(patient,'_VMphases.png');
  fullFileName = fullfile(figFolder,'Data',fileName);
  exportgraphics(f,fullFileName);
   
  bbutton = uicontrol('Parent',f,'Style','pushbutton',...
      'Position',[25,10,figSize(3)*0.25,figSize(4)*0.05],...
      'String','Save and exit','fontsize',fs,'Callback',@closeGenButton);
      
  uiwait(f)

% Save workspace
if nargin == 3
    s = strcat('../WS/',patids,patient,'_WS_temp.mat');
else
    s = strcat('../WS/',patient,'_WS_temp.mat');
end
save(s,'data','-append');

end