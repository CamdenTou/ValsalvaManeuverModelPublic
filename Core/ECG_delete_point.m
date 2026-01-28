function [TR,R,TS,S]=ECG_delete_point(TR,R,TS,S,qr,y) 
% function ECG_delete_point
% Input: location and value of R (TR, R) and S (TS,S) points within an ECG.
%        If y > 0 an R point is removed else an S point is removed.
% Output: augmented list of R and S points along with their location.
% Description: ECG_delete_point removes an R or S point from the list.

% TR - R time points (R>0)
% TS - S time points (S<0)

if y > 0 % the point is an erroneous R
    [~,i] = min(abs(TR-qr));
    TR = TR([1:(i-1) (i+1):end]);
    R  =  R([1:(i-1) (i+1):end]);
elseif y < 0 % point is an erroneous Q
    [~,i] = min(abs(TS-qr));
    TS = TS([1:(i-1) (i+1):end]);
    S  =  S([1:(i-1) (i+1):end]);
end

end