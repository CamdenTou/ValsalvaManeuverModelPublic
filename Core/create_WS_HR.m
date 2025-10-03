function [] = create_WS_HR(patient,plotmarkers)
% function create_WS_HR
% Input: patient name and figure settings
% Output: saved temporary workspace (saved in WS)
% Uses: correct_HR to correct heart rate values
% Description: loads temporary workspace with corrected ECG data calculates
% heart rate and saves temporary workspace with Heart rate output

% nested functions bc using user info (blue color is not global just bc of
% nested)

global HR_new % has to be global because of nested functions

%% Figure properties
fs = plotmarkers.fs;
lw = plotmarkers.lwt;
figSize = plotmarkers.figSize;

if lw > 2
    lw2 = lw-1;
else
    lw2 =1;
end

%% Load temporary workspace and assign vectors
filename = strcat(patient,'_WS_temp.mat');
load(strcat('../WS/',filename),'data','raw_data'); 
addpath('./fdaM/');

Traw   = raw_data.Traw0;    % time starts at zero
ECGraw = raw_data.ECGraw;

RRint   = data.RRint;
T_RRint = data.T_RRint;

HR      = 60./RRint;
%[HR, ~, ~, ~] = auto_smooth_hr(T_RRint, HR);

%% ECG to HR Calculations
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
    xlim([Traw(1),Traw(end)])
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
end % function close botton %

function [] = correctButton(a,~) 
    
    [HR_new] = correct_HR(HR,Traw,T_RRint,patient,plotmarkers);
    
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
        
    sub2 = subplot(2,1,2); hold on;
        plot(Traw,ECGraw,'b-')
        set(gca,'fontsize',fs)
        xlabel('Time (s)')
        ylabel('ECG (mV)')
        xlim([Traw(1),Traw(end)])
        %ylim([-1.25e-3 1.25e-3]);
        
    linkaxes([sub1, sub2],'x')

    figFolder = plotmarkers.figFolder;
    fileName = strcat(patient,'_HeartRateECG');
    fullFileName = fullfile(figFolder,'Data',fileName);
    saveas(f,fullFileName,'png');

    bbutton = uicontrol('Parent',f,'Style','pushbutton',...
    'Position',[25,10,figSize(3)*0.25,figSize(4)*0.05],...
    'String','Save and exit','fontsize',fs,'Callback',@closeGenButton);
   
end % function correct botton %

% Interpolate over step function and evaluate at Tdata
HRraw = interp1(T_RRint,HR_new,Traw,'pchip');

% Use HR to find RR intervals
RRraw = 60./HRraw;

%% Save temporary workspace
% data.HR     = HR;
% data.HR_new = HR_new;   
raw_data.HRraw = HRraw;   % user corrected HR data
data.RRraw     = RRraw;

plotmarkers.fs = fs;
plotmarkers.figSize = figSize;
plotmarkers.sub1 = sub1;
plotmarkers.sub2 = sub2;
plotmarkers.lw2 = lw2;

% Save temporary workspace
s = strcat('../WS/',patient,'_WS_temp.mat');
save(s,'data','raw_data','plotmarkers','-append');
end % function create_WS_HR %


function [smoothed_hr, best_nbasis, best_lambda, best_method] = auto_smooth_hr(t, hr_inst)
    % Candidate parameters
    candidate_nbasis = 3:1:min(15, floor(length(t)/2)); % number of basis functions in expansion
    candidate_lambdas = logspace(-8, 2, 50); % penalizing 2nd derivative
    best_gcv = inf;
    range = [min(t), max(t)];
    
    % Try different smoothing methods
    methods = {'bspline', 'fourier'};
    best_method = 'bspline';
    best_nbasis = [];
    best_lambda = [];
    
    for method_idx = 1:length(methods)
        method = methods{method_idx};
        
        switch method
            case {'bspline', 'fourier'}
                % Basis function approaches
                if strcmp(method, 'fourier')
                    % Fourier requires odd number of basis functions
                    nbasis_candidates = candidate_nbasis(mod(candidate_nbasis,2)==1);
                else
                    nbasis_candidates = candidate_nbasis;
                end
                
                for nbasis = nbasis_candidates
                    % Create basis
                    switch method
                        case 'bspline'
                            basis = create_bspline_basis(range, nbasis);
                        case 'fourier'
                            period = range(2) - range(1);
                            basis = create_fourier_basis(range, nbasis, period);
                        % case 'monom'
                        %     basis = create_monom_basis(range, nbasis);
                    end
                    
                    phi = eval_basis(t, basis);
                    R = eval_penalty(basis, 2);  % Second-derivative penalty
                    
                    for lambda = candidate_lambdas
                        A = phi' * phi + lambda * R;
                        c = A \ (phi' * hr_inst);
                        hr_hat = phi * c;
                        
                        % Compute GCV
                        S = phi / A * phi';
                        df = trace(S);
                        rss = sum((hr_inst - hr_hat).^2);
                        n = length(hr_inst);
                        gcv = (n * rss) / (n - df)^2;
                        
                        if gcv < best_gcv
                            best_gcv = gcv;
                            best_nbasis = nbasis;
                            best_lambda = lambda;
                            best_method = method;
                        end
                    end
                end
                
            % case 'pchip'
            %     % Piecewise Cubic Hermite Interpolating Polynomial
            %     hr_hat = pchip(t, hr_inst, t);
            % 
            %     % Compute GCV (approximate df as number of knots)
            %     df = length(t)/10;  % Empirical approximation
            %     rss = sum((hr_inst - hr_hat).^2);
            %     n = length(hr_inst);
            %     gcv = (n * rss) / (n - df)^2;
            % 
            %     if gcv < best_gcv
            %         best_gcv = gcv;
            %         best_nbasis = NaN;  % Not applicable
            %         best_lambda = NaN;  % Not applicable
            %         best_method = 'pchip';
            %     end
        end
    end
    
    % Apply best method
    switch best_method
        case 'bspline'
            basis = create_bspline_basis(range, best_nbasis);
            phi = eval_basis(t, basis);
            R = eval_penalty(basis, 2);
            c = (phi' * phi + best_lambda * R) \ (phi' * hr_inst);
            smoothed_hr = phi * c;
            
        case 'fourier'
            period = range(2) - range(1);
            basis = create_fourier_basis(range, best_nbasis, period);
            phi = eval_basis(t, basis);
            R = eval_penalty(basis, 2);
            c = (phi' * phi + best_lambda * R) \ (phi' * hr_inst);
            smoothed_hr = phi * c;
            
        case 'monom'
            basis = create_monom_basis(range, best_nbasis);
            phi = eval_basis(t, basis);
            R = eval_penalty(basis, 2);
            c = (phi' * phi + best_lambda * R) \ (phi' * hr_inst);
            smoothed_hr = phi * c;
            
        case 'pchip'
            smoothed_hr = pchip(t, hr_inst, t);
    end
end