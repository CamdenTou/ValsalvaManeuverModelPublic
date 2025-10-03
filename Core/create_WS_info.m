function [] = create_WS_info(patient,plotmarkers)
% function create_WS_info
% Input: patient name and figure settings
% Output: saved temporary workspace (saved in WS)
% Uses:
% Description: Selects data be analyzed and saves temporary workspace with
% patient characteristics and raw ECG, HR, and BP data.

%% Input patient ID, age, height, and weight
prompt1 = {'Patient ID Number','Age','Sex (male (m)/female(f))','Height (cm)','Weight (kg)'};
definput = {'123','25','f','165','77'};
dlgtitle = patient; 
dims = [1 75];

[answer1] = inputdlg(prompt1,dlgtitle,dims,definput);

% Cancel - return to previous screen
if sum(size(answer1) == [0 0]) == 2
    return
end

% Patient Info
PatientID   = str2num(answer1{1});
Age         = str2num(answer1{2});
Sex         = answer1{3};
Height      = str2num(answer1{4});
Weight      = str2num(answer1{5});

pat_char.PatientID  = PatientID;
pat_char.age        = Age;
pat_char.sex        = Sex;
pat_char.height     = Height;
pat_char.weight     = Weight;

%% Input channel indicies % SC added %
prompt2 = {'ECG Channel','HR Channel','BP Channel', ... 
   'Thoracic Pressure (if using enter channel, if not enter 0)'};
% side note about this prompt: if thoracic data is ommitted and it is
% before the other channels, just substitute the thoracic channel number
% with 0 (ex: Thoracic is first but want to omit: {'2','4','3','0'})
definput = {'2','4','3','1'};
dlgtitle = patient; 
dims = [1 75];

[answer2] = inputdlg(prompt2,dlgtitle,dims,definput);

% Cancel - return to previous screen
if sum(size(answer2) == [0 0]) == 2
    return
end

% Patient Info
ECG     = str2num(answer2{1});
HR      = str2num(answer2{2});
BP      = str2num(answer2{3});
PTh     = str2num(answer2{4});


%% Load LabChart ECG and blood pressure (BP) data
filename = strcat(patient,'.mat');
load(strcat('../Labchart/',filename)); % must be .mat file

start_time_of_file = 0;

% Indecies for channels for ECG, HR, BP (and PTh depending)
if PTh == 0
    channel_inds = [ECG, HR, BP];
else
    channel_inds = [ECG, HR, BP, PTh];
end

% Assemble data
tstart = 0;
sz = size(dataend);
endtime = dataend(1,1);
cols = sz(2);
if cols>1
    for i=2:cols
        endtime=endtime+(dataend(1,i)-dataend(end,i-1));
    end
end

dat = zeros(endtime,4);
if length(tickrate) ~= 1 
    if mean(tickrate) == tickrate(1)
        tickrate = tickrate(1);
    else
        error(strcat('Check tickrate - something strange. Index ',num2str(pt))) 
    end
end

dummy_time_vec = tstart:1/tickrate:endtime/tickrate-1/tickrate;
val_dat(:,1) = dummy_time_vec+start_time_of_file;
for j = 1:length(channel_inds)
    alldata = data(datastart(channel_inds(j),1):dataend(channel_inds(j),1));
    for i = 2:cols
        alldata=[alldata data(datastart(channel_inds(j),i):dataend(channel_inds(j),i))];
    end
    val_dat(:,j+1) = data(datastart(channel_inds(j)):dataend(channel_inds(j)));
end

% Create vector with time (Traw), ECG (ECGraw), heart rate (Hraw), BP
% (Praw), and thoracic pressure (Thraw)
Traw   = val_dat(:,1)-val_dat(1,1);
ECGraw = val_dat(:,2);
Hraw   = val_dat(:,3);
Praw   = val_dat(:,4);
if PTh ~= 0
    Pthraw = val_dat(:,5);
end

%% Plot data and select time series data for analysis

% figure properties
fs = plotmarkers.fs;
lw = plotmarkers.lwt;
figSize = plotmarkers.figSize;

f = figure(1); clf;
set(gcf,'units','points','position',figSize)
f.Units = 'pixels';

sgtitle(patient,'Fontsize',fs,'FontWeight','bold','Interpreter','none');
subplot(3,1,1)
    plot(Traw,ECGraw,'b')
    set(gca,'fontsize',fs)
    title('Select Data (~20s before and after VM)') 
    ylabel('ECG (mV)')
    xlim([0 max(Traw)])
subplot(3,1,2)
    plot(Traw,Hraw,'b','linewidth',lw)
    set(gca,'fontsize',fs)
    ylabel('HR (bpm)')
    xlim([0 max(Traw)])
    ylim([min(Hraw)-2 max(Hraw)+2])
subplot(3,1,3)
    plot(Traw,Praw,'b')
    set(gca,'fontsize',fs)
    xlabel('Time (s)')
    ylabel('BP (mmHg)')
    xlim([0 max(Traw)])
    ylim([min(Praw)-5 max(Praw)+5])

% Establish data start and end points using user input
[x,~] = ginput(2);
if x(1)<0 == 1
    istart = 1;     % Special case if x(1) is negative
else
    istart = max(find(Traw<x(1)));
end
iend   = max(find(Traw<x(2)));

% Save selected data
TrawO  = Traw;
Traw   = Traw(istart:iend);
ECGraw = ECGraw(istart:iend);
Hraw   = Hraw(istart:iend);
Praw   = Praw(istart:iend);
if PTh ~= 0
    Pthraw = Pthraw(istart:iend);
end


figure(1); hold on
sgtitle(patient,'Fontsize',fs+2,'FontWeight','bold','Interpreter','none');
subplot(3,1,1); hold on
    plot(Traw,ECGraw,'r')
    set(gca,'fontsize',fs)
    title('LabChart data')
    ylabel('ECG (mV)')
    xlim([0 max(TrawO)])
subplot(3,1,2); hold on
    plot(Traw,Hraw,'r','linewidth',lw)
    ylabel('HR (bpm)')
    set(gca,'fontsize',fs)
    xlim([0 max(TrawO)])
    ylim([min(Hraw)-10 max(Hraw)+2])
subplot(3,1,3); hold on
    plot(Traw,Praw,'r')
    xlabel('Time (s)')
    ylabel('BP (mmHg)')
    xlim([0,max(TrawO)])
    ylim([min(Praw)-5 max(Praw)+5])

figFolder = plotmarkers.figFolder;
fileName = strcat(patient,'_dataAnalyzed.png');
fullFileName = fullfile(figFolder,'Data',fileName);
exportgraphics(f,fullFileName);

bbutton = uicontrol('Parent',f,'Style','pushbutton',...
      'Position',[25,10,figSize(3)*0.25,figSize(4)*0.05],...
      'String','Save and exit','fontsize',fs,'Callback',@closeGenButton);
 
uiwait(f)

%% Save selected data to structure

raw_data.TrawOrg    = Traw;
raw_data.Traw       = Traw-Traw(1);
raw_data.ECGraw     = ECGraw;
raw_data.Hraw       = Hraw;
raw_data.Praw       = Praw;
raw_data.PthAvail   = PTh;
if PTh ~= 0
    raw_data.Pthraw     = Pthraw;
end


%% Save temmporary workspace 
s = strcat('../WS/',patient,'_WS_temp.mat');
save(s,'raw_data','pat_char');

end  % function create_WS_info.m %

