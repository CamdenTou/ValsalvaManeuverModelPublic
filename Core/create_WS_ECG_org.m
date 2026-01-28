function [] = create_WS_ECG(patient,plotmarkers)
% function create_WS_ECG
% Input: patient name and figure settings
% Output: saved temporary workspace (saved in WS)
% Uses: correct_ECG to correct the heart rate.
% Description: loads temporary workspace with raw data, corrects the ECG
% removing spurious peaks, calculates RR interval and saves corrected ECG
% in the temporary workspacce. This function also saves R and S peaks.

%% Load temporary workspace
filename = strcat(patient,'_WS_temp.mat');
load(strcat('../WS/',filename)); % load patient info

Traw   = raw_data.Traw;
ECGraw = raw_data.ECGraw;
Hraw   = raw_data.Hraw;
Praw   = raw_data.Praw;

Traw = Traw -Traw(1); % ensure that the time starts at zero
raw_data.Traw0 = Traw;

dt   = mean(diff(Traw)); 
raw_data.dt = dt;

%% Determines baseline of ECG signal with medfilt1

%Filter out P waves and QRS complex with a window of 200 ms
smoothECG = medfilt1(ECGraw,round(.2/dt)); 
 
%Filter out T waves with a window of 600 ms 
smoothECG2 = medfilt1(smoothECG,round(.6/dt)); 

%Baseline corrected ECG signal
BaselineCorrectedECG = ECGraw - smoothECG2; 

%% Savitsky-Golay Filter 

%Savitsky-Golay Filter filters out very low frequency signals. The window
%must be odd and within .15 seconds
SVGwindow = round(.15/dt); 
if mod(SVGwindow,2) == 0
    SVGwindow = SVGwindow + 1;
end 
%Check to ensure order is less than window size 
if SVGwindow > 5
    SVGorder = 3; 
else
    SVGorder = 3; 
end 
smoothECG3 = sgolayfilt(BaselineCorrectedECG,SVGorder,SVGwindow); 

%% Accentuate peaks to easily find them 

% Extract Q and S peaks 
accentuatedpeaks = BaselineCorrectedECG - smoothECG3; 

% Finds Q and S points with minimum peak distance of 200 ms 
[~,z]  = findpeaks(accentuatedpeaks,'MinPeakDistance',round(.2/dt)); 
zz     = mean(accentuatedpeaks(z)); 
[~,iR] = findpeaks(accentuatedpeaks,'MinPeakHeight',zz,...
    'MinPeakDistance',round(.2/dt)); 
[~,y]  = findpeaks(-accentuatedpeaks,'MinPeakDistance',round(.2/dt)); 
yy     = -mean(accentuatedpeaks(y));
[~,iQ] = findpeaks(-accentuatedpeaks,'MinPeakHeight',yy,...
    'MinPeakDistance',round(.2/dt)); 

% R peaks and S peaks of ECG waveforms needed to extract respiration
R  = ECGraw(iR);
TR = Traw(iR); 
S  = ECGraw(iQ);
TS = Traw(iQ); 

% Check to ensure all peaks are >= 25% of the mean or <= 200% of the mean
mR = mean(R);
mQ = mean(S);
xR = find(R <= 2*mR & R >= .25*mR);
iR = iR(xR);
xQ = find(S >= 2*mQ & S <= .25*mQ);
iQ = iQ(xQ);
R  = R(xR);
TR = TR(xR);
S  = S(xQ);
TS = TS(xQ);

nR = length(R);
nS = length(S);

% Corrects the ECG removing spurious markers and adding missing markers
[TR,TS,R,S] = correct_ECG(patient,Traw,ECGraw,TR,TS,R,S,plotmarkers);

T_RRint = sort(TR(1:end-1));
RRint   = diff(sort(TR));

data.R  = R;
data.TR = TR;
data.S  = S;
data.TS = TS;
data.RRint = RRint;
data.T_RRint = T_RRint;

% Save temporary workspace
s = strcat('../WS/',patient,'_WS_temp.mat');
save(s,'raw_data','data','-append');
end % function create_WS_ECG.m %