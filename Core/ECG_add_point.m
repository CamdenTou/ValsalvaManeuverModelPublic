function [TR,R,TS,S,i]=ECG_add_point(TR,R,TS,S,tc,y) 
% function ECG_add_point
% Input: location and value of R (TR, R) and S (TS,S) points within an ECG.
%        If y > 0 an R point is added else an S point is added.
% Output: augmented list of R and S points along with their location, and
% an index i = 0 if the point added is the first point.
% Description: ECG_add_point adds a R or S point to the list.

% TR - R time points (R>0)
% TS - S time points (S<0)
if y > 0 % the point is R
    if tc < TR(1)
        TR = [tc TR(1:end)']'; 
        R  = [y  R(1:end)']';
        i = 0;
    else
        [~,i] = max(find(TR<tc));
        TR = [TR(1:i)' tc TR(i+1:end)']';
        R  = [R(1:i)' y R(i+1:end)']';
    end
elseif y < 0 % point is S
    if tc < TS(1)
       TS  = [tc TS(1:end)']'; 
        S  = [y S(1:end)']';
        i = 0;
    elseif tc < TS(end)
        [~,i] = max(find(TS<tc));
        TS = [TS(1:i)' tc TS(i+1:end)']';
        S  = [S(1:i)'  y S(i+1:end)']';
    else
        TS = [TS(1:end)' tc]';
        S  = [S(1:end)'  y]';
        i  = length(TS);
    end
end

end