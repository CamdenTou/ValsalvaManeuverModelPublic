function per_p=get_periods(t, ch1, no_reg, p_peak_thresh,patient,plotmarkers)
%returns the indices of the period boundaries of the data in fnm.dat
%also requires number of regions into which the signal is to be divided (no_reg)
%no_reg is usually 1 or 3
%p_peak_thresh is whatever fraction of the peaks clicked on in a given region
%another peak must be for detection
%per_p is for pressure and per_v is for flow

thresh=-150;
t_bound=max(t);

fs = plotmarkers.fs;
lw = plotmarkers.lwt;
if lw > 2
    lw2 = lw - 1;
else
    lw2 = 1;
end
figSize = plotmarkers.figSize;

%partition data into regions
ch1_cell=cell(3,no_reg); %cells for storing partitioned data
%time in row 1, p in row 2, starting index in row 3, periods vector in row 4
ch1_mat{1,1}=t(find(t<=t_bound(1)));
ch1_mat{2,1}=ch1(find(t<=t_bound(1))); 
ch1_mat{3,1}=0;                        
for ind=2:no_reg
    ch1_mat{1,ind}=t(find((t>t_bound(ind-1)) & (t<=t_bound(ind))));
    ch1_mat{2,ind}=ch1(find((t>t_bound(ind-1)) & (t<=t_bound(ind))));
    ch1_mat{3,ind}=max(find(t<=t_bound(ind-1)));    
end

% Selecting the distance between peaks
f = figure(1);clf;
set(gcf,'units','points','position',figSize)
f.Units = 'pixels';

for ind=1:no_reg
    xID = min(ch1_mat{1,ind})+.2*(max(ch1_mat{1,ind}));
    ID = max(find(ch1_mat{1,ind}<xID));
    subplot(no_reg,1,ind);
    plot(ch1_mat{1,ind}(1:ID),ch1_mat{2,ind}(1:ID),'b','Linewidth',lw2);
    axis([min(ch1_mat{1,ind}(1:ID)) max(ch1_mat{1,ind}(1:ID)) min(ch1_mat{2,ind}(1:ID))*0.95 max(ch1_mat{2,ind}(1:ID))*1.05]);
end
title(patient,'Click on two consecutive peaks (from left to right)','interpreter','none');
xlabel('Time (s)')
ylabel('BP (mmHg)');
set(gca,'fontsize',fs)

[t_peaks p_peaks]= ginput(2*no_reg);

% Determine pts for each region and compute period vectors
for ind=1:no_reg
    pts(ind)=ceil(.4*length(find((ch1_mat{1,ind}>=t_peaks(2*ind-1)) & (ch1_mat{1,ind}<=t_peaks(2*ind)))));
    ch1_mat{4,ind}=peaks2(ch1_mat{2,ind},pts(ind),p_peak_thresh*p_peaks(2*ind))+ch1_mat{3,ind};
    ch1_mat{5,ind}=peaks2(-ch1_mat{2,ind},pts(ind),thresh)+ch1_mat{3,ind};
end
        
% Create period boundary vectors
per_p=[]; peak_p=[]; 
for ind=1:no_reg
    peak_p=[peak_p; ch1_mat{4,ind}];
    per_p=[per_p; ch1_mat{5,ind}];
end

% Deliminates period boundaries that are very close
bad_per_p=0; 
while(bad_per_p)
    diff_p=diff(per_p); 
    min_diff_p=.01*mean(diff(ch1_mat{5,1})); 
    bad_per_p=min(find(diff_p<=min_diff_p)); 
    if(bad_per_p)
      per_p=per_p([1:bad_per_p (bad_per_p+2):length(per_p)]);
    end
end

per_pc=per_p; % Copies of period vectors
clear per_p 
for ind=1:length(peak_p)
    if(length(per_pc(max(find(per_pc<peak_p(ind)))))>0)
        per_p(ind)=per_pc(max(find(per_pc<peak_p(ind))));
    end
end

bad_per_p=1;
while(bad_per_p) 
    diff_p=diff(per_p);
    min_diff_p=.01*mean(diff(ch1_mat{5,1}));
    bad_per_p=min(find(diff_p<=min_diff_p));
    if(bad_per_p)
      per_p=per_p([1:bad_per_p (bad_per_p+2):length(per_p)]);
    end
end

%gets rid of repeats in p data, which happen for some reason...
diff_p=diff(per_p);
diff_p=diff_p(find(diff_p));
per_p=[cumsum([per_p(1) diff_p])];

per_p=per_p(find(per_p)); %assures only nonzero periods are returned

end