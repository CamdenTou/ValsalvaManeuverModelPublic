function [] = create_WS_test(patient,plotmarkers)

    disp(patient)
    disp(plotmarkers)

    fs = plotmarkers.fs;
    lw = plotmarkers.lwt;
    figSize = plotmarkers.figSize;

if lw > 2
    lw2 = lw-1;
else
    lw2 =1;
end

filename = strcat(patient,'_WS_temp.mat');
load(strcat('../WS/',filename)); 

Traw   = raw_data.Traw0;
ECGraw = raw_data.ECGraw;

RRint   = data.RRint;
T_RRint = data.T_RRint;

HR      = 60./RRint;

f = figure(1); clf; hold on;
set(gcf,'units','points','position',figSize)
f.Units = 'pixels';
sgtitle(patient,'Fontsize',fs+2,'FontWeight','bold','Interpreter','none')
sub1 = subplot(2,1,1); hold on;
    plot(T_RRint,HR,'bo-','linewidth',lw2);
    set(gca,'fontsize',fs)
    ylabel('HR (bpm)')
    xlim([Traw(1),Traw(end)])
    ylim([min(HR)-5 max(HR)+5]);

sub2 = subplot(2,1,2); hold on;
    plot(Traw,ECGraw,'b-')
    set(gca,'fontsize',fs)
    xlabel('Time (s)')
    ylabel('ECG (mV)')
    ylim([-1.25e-3 1.25e-3]);
    xlim([Traw(1),Traw(end)])
    ylim([-1.25e-3 1.25e-3]);
linkaxes([sub1, sub2],'x')

figFolder = plotmarkers.figFolder;
fileName = strcat(patient,'_HeartRateECG');
fullFileName = fullfile(figFolder,'Data',fileName);
saveas(f,fullFileName,'png');

bbutton = uicontrol('Parent',f,'Style','pushbutton',...
    'Position',[25,10,figSize(3)*0.25,figSize(4)*0.05],...
    'String','Save and exit','fontsize',fs,'Callback',@closeButton);

cbutton = uicontrol('Parent',f,'Style','pushbutton',...
    'Position',plotmarkers.posCB,...
    'String','Correct Heart Rate','fontsize',fs,'Callback',@correctButton);

uiwait(f)


function [] = closeButton(a,~)
    f = a.Parent;
    HR_new = HR;
    close(f)
end

function [] = correctButton(a,~) 
    
    [HR_new,answer_HR] = correct_HR(HR,Traw,T_RRint,patient,plotmarkers);
    
    f = figure(1); clf; hold on;
    set(gcf,'units','points','position',figSize)
    f.Units = 'pixels';
    sgtitle(patient,'Fontsize',fs+2,'FontWeight','bold','Interpreter','none')
    sub1 = subplot(2,1,1); hold on;
        plot(T_RRint,HR_new,'bo-','linewidth',lw2);
        set(gca,'fontsize',fs)
        ylabel('HR (bpm)')
        xlim([Traw(1),Traw(end)])
        ylim([min(HR_new)-5 max(HR_new)+5]);


end