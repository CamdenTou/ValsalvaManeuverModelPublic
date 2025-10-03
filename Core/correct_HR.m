function [HR_new,answer_HR] = correct_HR(HR,Traw,T_RRint,patient,plotmarkers)
% function correct_HR
% Input: HR data, time vector, RR intervals, patient name and plot markers
% Output: corrected heart rates and response to dialogue box.
% Uses: spline_HR.m to concatenate corrected HR points
% Description: Reads ECG data and corrects the HR if needed.

%% Figure properties
cfs = plotmarkers.fs;
figSize = plotmarkers.figSize;

answer_HR = 'Add change';
HR_new = HR;

%% Correcting heart rate 
while strcmp(answer_HR,'Add change') == 1

    f = figure(1); hold on
    subplot(2,1,1); hold on;
    title('scrolling (enter when done)')
    pan on
    xlim([Traw(1), 40])
    pause;

    if ismember(findall(0,'type','figure'),f) == 0
        return
    end
    
    uiwait(msgbox("Click at points to connect",'modal'));
    
    if ismember(findall(0,'type','figure'),f) == 0
        return
    end
    
    [x,~] = ginput(2);
    if ~isempty(x) == 1
        [HR_new, answer_HR] = spline_HR(patient,HR_new,T_RRint,Traw,x,plotmarkers);
    else
        return
    end

end

end % function correct_HR %