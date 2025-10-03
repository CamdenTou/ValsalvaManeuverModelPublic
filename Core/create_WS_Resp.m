function [] = create_WS_Resp(patient,plotmarkers)
% function create_WS_Resp
% Input: patient name and figure settings
% Output: saved temporary workspace (saved in WS)
% Description: Generates respiratory signal from ECG data using algorithm
% from Randall et al. (ADD REFERENCE)

%% Load temporary workspace and assign vectors
filename = strcat(patient,'_WS_temp.mat');
load(strcat('../WS/',filename));

Traw = raw_data.Traw0;
R    = data.R;
TR   = data.TR;
S    = data.S;
TS   = data.TS;

%% Figure properties
fs = plotmarkers.fs;
lw = plotmarkers.lwt;
figSize = plotmarkers.figSize;

%% Calculate respiratory signal
dt = mean(diff(Traw)); 
Fs = round(1/dt); 

% Find amplitude as difference of R- and S-waves
RS  = R - S;
TRs = sort(TR);
TQs = sort(TS);

% Take average time between Q-and R-waves 
TRS = (TRs + TQs)/2;

% Check for positive amplitudes 
xRS = find(RS > 0); 
RS  = RS(xRS); 
TRS = TRS(xRS); 

% Check to ensure amplitudes are <= 25% or >= 200% mean amplitude
mRS = mean(RS); 
xRS = find(RS >= .25*mRS & RS <= 2*mRS);
RS  = RS(xRS);
TRS = TRS(xRS);

%% Interpolate through amplitude points 

% Include first and last time points for piecewise cubic Hermite splines 
RS  = [RS(1); RS; RS(end)];
TRS = [Traw(1); TRS; Traw(end)];

% Use PCHIP algorithm to interpolate through amplitudes (this algorithm
% preserves derivative behavior)
S = griddedInterpolant(TRS,RS,'pchip'); 
QRamplitudes = S(Traw); % Evaluate at the time points 

%% Find peaks and valleys minimum peak distance of 1.5 seconds
[p,ploc] = findpeaks(QRamplitudes,Fs,'MinPeakDistance',1.5); 
[v,vloc] = findpeaks(-QRamplitudes,Fs,'MinPeakDistance',1.5);
v = -v;

% Reorganize time vector and sort filtered peaks and valleys of QR amp signal
Tloc   = [ploc; vloc]; 
[Tloc,ii] = sort(Tloc,'ascend'); 
Points = [p; v]; 
Points = Points(ii);
Tloc   = [Traw(1); Tloc; Traw(end)];
Points = [QRamplitudes(1); Points; QRamplitudes(end)];
ss     = griddedInterpolant(Tloc,Points,'pchip'); 

[Ex,T_Ex] = findpeaks(ss(Traw),Fs); 
[In,T_In] = findpeaks(-ss(Traw),Fs); 
In = -In; 

% Check to ensure first time point is included
if T_Ex(1) < T_In(1) 
    T_In = [Traw(1); T_In]; 
    In = [ss(Traw(1)); In]; 
end 
% Check to ensure last time point is included 
if T_In(end) > T_Ex(end) 
    T_Ex = [T_Ex; Traw(end)];
    Ex = [Ex; ss(Traw(end))];
end 

%% Impose restrictions from Widjaja paper 

% Amplitude of each breath 
Ramp = Ex - In;
% Difference in previous and subsequent amplitudes 
surroundingampdiff = abs(Ramp(3:end) - Ramp(1:end-2));
% Find 15% amplitude difference 
Ramp15diff = .15*surroundingampdiff; 
F = [Ramp(1); Ramp15diff; Ramp(end)]; 
% Remove amplitudes less than 15% of difference of the surrounding
% amplitudes
ii = [];
for i = 1:length(Ramp)
    if Ramp(i) >= F(i)
        ii = [ii i];
    end 
end 

%Combine filtered expiration and inspiration signals and sort
Ex   = Ex(ii);
T_Ex = T_Ex(ii);
In   = In(ii);
T_In = T_In(ii);
[newT,i] = sort([T_Ex; T_In],'ascend'); 
newR = [Ex; In];
newR = newR(i); 

% Tack on end time points for PCHIP interpolation
if isempty(find(newT == Traw(1),1)) == 1
    newT = [Traw(1); newT];
    newR = [ss(Traw(1)); newR];
end 
if isempty(find(newT == Traw(end),1)) == 1
    newT = [newT; Traw(end)];
    newR = [newR; ss(Traw(end))]; 
end
if isempty(find(newT == Traw(1) | newT == Traw(end),1)) == 1 
    newT = [Traw(1); newT; Traw(end)];
    newR = [ss(Traw(1)); newR; ss(Traw(end))]; 
end
Rspline = griddedInterpolant(newT,newR,'pchip');
Rraw   = Rspline(Traw); 

%% Plot results
f = figure(1); clf; hold on;
set(gcf,'units','points','position',figSize)
f.Units = 'pixels';
title(patient,'Interpreter','none'); hold on;
plot(Traw,Rraw,'b-','Linewidth',lw);
set(gca,'fontsize',fs)
xlim([Traw(1),Traw(end)])
ylabel('Respiration')
xlabel('Time (s)')

figFolder = plotmarkers.figFolder;
fileName = strcat(patient,'_Respiration.png');
fullFileName = fullfile(figFolder,'Data',fileName);
exportgraphics(f,fullFileName);

bbutton = uicontrol('Parent',f,'Style','pushbutton',...
    'Position',[25,10,figSize(3)*0.25,figSize(4)*0.05],...
    'String','Save and exit','fontsize',fs,'Callback',@closeGenButton);

uiwait(f)

%% Save temporary workspace
data.Rraw = Rraw;

s = strcat('../WS/',patient,'_WS_temp.mat');
save(s,'data','-append');

end % function create_WS_Resp.m %
