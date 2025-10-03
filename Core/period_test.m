function [SPinds, DPinds] = period_test(Traw,Praw,patient,plotmarkers)
% Uses get_periods2 to create files that contain period boundaries
% may or may not need to comment out -eps difference in get_periods2
% just see what works best

pN = Praw;
pd  = csaps(Traw,pN,0.999999,Traw);

per_p = get_periods(Traw, pd, 1, 0, patient,plotmarkers);

fs = plotmarkers.fs;
lw = plotmarkers.lwt;
figSize = plotmarkers.figSize;

id1=1;
for i = 1:length(per_p)-1
    id2 = max(find(Traw<=Traw(per_p(i+1))));
    pAomean(i)  = 1/(Traw(per_p(i+1))-Traw(per_p(i)))*trapz(Traw(id1:id2),pd(id1:id2));
    id1 = id2;

    RAo  = [per_p(i) per_p(i+1)];
    
    [maxpAo idxAo] = max(pd(RAo(1):RAo(2)));
    mINDAo(i)      = RAo(1)+idxAo-1;    
end


HR  = 60./diff(Traw(per_p));
tHR = Traw(per_p(1:end-1));
t_per = Traw(per_p);
p_per = Praw(per_p);

SPinds = mINDAo;
DPinds = per_p;

