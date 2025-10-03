function [HR_new,answer_HR] = spline_HR(patient,HR,T_RRint,Traw,x,plotmarkers)
% function spline_HR 
% Input: patient name, HR data, time vector, RR intervals, list of points to be connected, and plot markers
% Output: corrected heart rates and response to dialogue box.
% Uses: N/A
% Description: inserts a spline merging orginal and corrected heart rate
% values. Returns a figure with corrected heart rates saved to
% Figures/WS_data

% Generate new heart rate vector
HR_new = HR;

% Identify indicies in the time-vector where the corrected heart rate
% should be inserted
for i = 1:length(x)
    [~,i_xR(i)] =  min(abs(T_RRint-x(i))); %find index for RR
end

lw = plotmarkers.lwt;
if lw > 2
    lw2 = lw -1;
else
    lw2 = 1;
end

% Linear spline
tA = T_RRint(i_xR(1)); 
tB = T_RRint(i_xR(2)); 
yA = HR(i_xR(1));
yB = HR(i_xR(2));
a  = (yB-yA)/(tB-tA);
b  = yB - a*tB;
t  = T_RRint(i_xR(1):i_xR(2));
Ht = a*t+b;

% Update heart rate vector
HR_new(i_xR(1):i_xR(2)) = Ht;
    
% Plot corrected heart rate and save the figure
f = figure(1); hold on
subplot(2,1,1);hold on
g = plot(t,Ht,'Og-','LineWidth',lw2);
ylim([min(HR)-5 max(HR)+5])
xlim([min(T_RRint) max(T_RRint)])

subplot(2,1,2);hold on
xlim([min(Traw) max(Traw)])

answer_HR = questdlg('Accept changes?','Heart Rate','Yes','Undo','Add change','Yes');

if strcmp(answer_HR,'Yes') == 1
    return
elseif strcmp(answer_HR,'Undo') == 1
    HR_new = HR;
    figure(1); hold on
    delete(g)

    answer_HR = questdlg('Accept changes?','Heart Rate','Yes','Add change','Yes');

elseif strcmp(answer_HR,'Add change') == 1
    return
end

end

