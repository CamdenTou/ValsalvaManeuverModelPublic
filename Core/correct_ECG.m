function [TR,TS,R,S] = correct_ECG(patient,Traw,ECGraw,TR,TS,R,S,plotmarkers)
% function correct_ECG
% Input: patient name, vectors with time (Traw), ECG signal (ECGraw) times
% and values for R and S peaks as well as figure setting (plotmarkers)
% Output: corrected vectors and times for R and S peaks
% Interactive function removing spurious and adding missing R and S peaks.
% Uses: ECG_add_point and ECG_delete_point to add and remove wrong ECG R
% and S points.
% Description: Interactive function removing spurious and adding missing R and S peaks.

%% Figure properties
fs = plotmarkers.fs;
ms = plotmarkers.ms; 
lw = plotmarkers.lwt;
figSize = plotmarkers.figSize;
if lw > 2
    lw2 = lw-1;
else
    lw2 =1;
end
if lw > 2
    lw3 = 1;
else
    lw3 = 0.75;
end

%% Loop correcting ECG
notdone = true;
while notdone
    
    nR = length(R);
    nS = length(S);
    
    f = figure(1);clf; hold on
    set(gcf,'units','points','position',figSize)
    f.Units = 'pixels';

    h1 = plot(Traw,ECGraw,'k-','linewidth',lw3);
    h2 = plot(TR,R,'ro','markersize',ms,'linewidth',lw2);
    h3 = plot(TS,S,'bo','markersize',ms,'linewidth',lw2);
    hold off
    set(gca,'FontSize',fs)
    xlim([Traw(1), 15])
    xlabel('Time (s)')
    ylabel('Raw ECG (mV)')
    title(patient,'interpreter','none')
    str1 = {'Number of S points: ',num2str(nS)};
    str2 = {'Number of R points: ',num2str(nR)};
    a1 = annotation('textbox', [0.91, 0.4, 0, 0], 'string', str1,'fontsize',fs);
    a2 = annotation('textbox', [0.91, 0.9, 0, 0], 'string', str2,'fontsize',fs);
    xL=xlim;
    yL=ylim;
    pan on
    bgcolor = f.Color;
    %pause;
    
    % Clicking correction code
    answer_ECG = questdlg('Do you want to correct RS points?', ...
        'ECG wave correction','Yes','No','Yes');

    answer_repeat = 'Yes';
    if strcmp(answer_ECG,'Yes')==1
        while strcmp(answer_repeat,'Yes') == 1
            f = figure(1); hold on
            title(patient,'scrolling (enter when done)')
            pan on
            pause;  

            title(patient,'Place cross-hair on point to change')
            
            if ismember(findall(0,'type','figure'),f) == 0
                return
            end

            [qr,y] = ginput(1);
            answer_QS = questdlg('Add or remove point?', ...
                'ECG wave correction','Add','Remove','Cancel','Cancel');
            
            if strcmp(answer_QS,'Remove') == 1
                
                [TR,R,TS,S] = ECG_delete_point(TR,R,TS,S,qr,y);
                
                figure(1); hold on
                delete(h2);
                delete(h3);
                h2 = plot(TR,R,'ro','markersize',ms,'linewidth',lw2);
                h3 = plot(TS,S,'bo','markersize',ms,'linewidth',lw2);
                
            elseif strcmp(answer_QS,'Add') == 1
                [TR,R,TS,S,new] = ECG_add_point(TR,R,TS,S,qr,y);
                
                figure(1);hold on
                a_new = plot(TR(new+1),R(new+1),'o','color','#77AC30',...
                    'markersize',ms,'linewidth',lw);
                pause(1);
            end
            nR = length(R);
            nS = length(S);
            
            figure(1); pan on
            delete(a1);
            delete(a2);
            
            str1 = {'Number of S points: ',num2str(nS)};
            str2 = {'Number of R points: ',num2str(nR)};
            a1 = annotation('textbox', [0.91, 0.4, 0, 0], 'string', str1,'fontsize',fs);
            a2 = annotation('textbox', [0.91, 0.9, 0, 0], 'string', str2,'fontsize',fs);
            
            title(patient,'Continue inspecting points - press enter when done');
            
            figure(1); hold on
            if exist('a_new','var') == 1
                delete(a_new)
            end
            delete(h2);
            delete(h3);
            h2 = plot(TR,R,'ro','markersize',ms,'linewidth',lw2);
            h3 = plot(TS,S,'bo','markersize',ms,'linewidth',lw2);            
            pause
            
            answer_repeat = questdlg('Do you want to correct RS points?', ...
                'ECG wave correction','Yes','No','Yes');
        end
    end
    
    % Check to ensure lengths of vectors are the same
    if nR == nS
        notdone = false;
    else
        msg_not_equal = {'Number of R and S peaks must be equal. ',...
            'Inspect data.'};
        uiwait(msgbox(msg_not_equal,'Error',"modal"));
    end
    
end % while loop %

end % function correct_ECG % 
