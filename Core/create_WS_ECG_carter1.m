cdfunction [] = create_WS_ECG6(patient, plotmarkers)
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

Traw = Traw - Traw(1); % ensure that the time starts at zero
raw_data.Traw0 = Traw;

dt   = mean(diff(Traw)); 
fs = 1/dt;
raw_data.dt = dt;

%% Fixed Parameter Selection
baseline_window = 0.2;     % Baseline removal window (seconds)
sg_window = 0.1;           % Savitzky-Golay window (seconds)
sg_order = 3;              % Savitzky-Golay polynomial order
bandpass_range = [1 40]; % Bandpass frequencies (Hz)
wavelet = 'db4';           % Wavelet type
level = 6;                 % Wavelet decomposition level
energy_threshold = 0.99;   % Energy threshold for reconstruction

%% 1. Baseline Removal
% win_samples1 = round(baseline_window / dt);
% smooth1 = medfilt1(ECGraw, win_samples1, 'omitnan', 'truncate');
% win_samples2 = round(3*baseline_window / dt);
% smooth2 = medfilt1(smooth1, win_samples2, 'omitnan', 'truncate');
% ECG_baseline = ECGraw - smooth2;

%% 2. Savitzky-Golay Smoothing
% win_sg = round(sg_window / dt);
% if mod(win_sg,2) == 0
%     win_sg = win_sg + 1;  % Ensure odd window length
% end
% ECG_sg = sgolayfilt(ECG_baseline, sg_order, win_sg);

%% 3. Bandpass Filtering
[b,a] = butter(4, bandpass_range/(fs/2), 'bandpass');
ECG_band = filtfilt(b, a, ECGraw);

%% 4. Outlier Removal
Q1 = quantile(ECG_band, 0.25);
Q3 = quantile(ECG_band, 0.75);
IQRv = Q3 - Q1;
LB = Q1 - 4*IQRv;
UB = Q3 + 4*IQRv;
ECG_clean = ECG_band;

for i = 1:length(ECG_band)
    if ECG_band(i) < LB || ECG_band(i) > UB
        idxs = max(1,i-2):min(length(ECG_band), i+2);
        neighbors = ECG_band(idxs);
        neighbors(neighbors == ECG_band(i)) = [];
        ECG_clean(i) = mean(neighbors, 'omitnan');
    end
end

%% 5. Wavelet Reconstruction 
[c, l] = wavedec(ECG_clean, level, wavelet);
energy = abs(c).^2;
total_energy = sum(energy);
cumE = cumsum(sort(energy, 'descend'));
idx_keep = find(cumE >= energy_threshold*total_energy, 1);
[~, sorted_idx] = sort(energy, 'descend');
c_thresh = zeros(size(c));
c_thresh(sorted_idx(1:idx_keep)) = c(sorted_idx(1:idx_keep));
ECG_recon = waverec(c_thresh, l, wavelet);

%% 6. Peak Detection
% Estimate refractory period
if isfield(raw_data,'HRraw')
    idx10 = find(Traw > 10, 1);
    if ~isempty(idx10) && idx10 > 1
        HRavg_first10 = mean(raw_data.HRraw(1:idx10));
        num_beats_first10 = HRavg_first10 / 6;
        average_refractory_first10 = 10 / num_beats_first10;
    else
        average_refractory_first10 = 0.8;
    end
else
    average_refractory_first10 = 0.8;
end

min_peak_dist_samples = round(0.5 * average_refractory_first10 * fs);
win_s_sec = 0.2; win_q_sec = 0.3;
win_s = round(win_s_sec * fs); win_q = round(win_q_sec * fs);
th = 0.30 * max(ECG_recon);
[~, locs] = findpeaks(ECG_recon, 'MinPeakDistance', min_peak_dist_samples, 'MinPeakHeight', th);

% Gap filling
TR_initial = Traw(locs); 
RR = diff(TR_initial);
RR_thresh = 1.05 * average_refractory_first10;
baseline_RR = 0.8;

% Prepare storage for new interpolated peaks
new_times = [];
new_vals = [];

% (Optional) store for debugging
data.test = TR_initial;
data.gap_thresh = RR_thresh;

for k = 1:length(RR)
    if RR(k) > RR_thresh
        t_left = TR_initial(k); 
        t_right = TR_initial(k+1);
        sig_left = ECG_recon(locs(k)); 
        sig_right = ECG_recon(locs(k+1));
        gap_len = t_right - t_left;

        % --- Local RR before gap
        % Find last peak time <= t_left - 2.5 sec
        idx_pre = find(TR_initial <= t_left - 2.5, 1, 'last');
        if isempty(idx_pre) || idx_pre >= k
            RR_before = baseline_RR;
        else
            pre_peaks = locs(idx_pre:k);
            if length(pre_peaks) > 1
                RR_before = mean(diff(Traw(pre_peaks)));
            else
                RR_before = baseline_RR;
            end
        end

        % --- Local RR after gap
        % Find first peak time >= t_right + 2.5 sec
        idx_post = find(TR_initial >= t_right + 2.5, 1, 'first');
        if isempty(idx_post) || idx_post <= k+1
            RR_after = baseline_RR;
        else
            % idx_post is an index into TR_initial; but we want peaks from k+1 up to idx_post
            % Convert idx_post in TR_initial to index in locs (they align: TR_initial(i) = Traw(locs(i)))
            post_peaks = locs((k+1):idx_post);
            if length(post_peaks) > 1
                RR_after = mean(diff(Traw(post_peaks)));
            else
                RR_after = baseline_RR;
            end
        end

        % Validate RR estimates against a max valid RR
        max_valid_RR = 1.05 * average_refractory_first10;
        if RR_before > max_valid_RR, RR_before = baseline_RR; end
        if RR_after > max_valid_RR, RR_after = baseline_RR; end

        % Determine final RR for interpolation spacing
        final_RR = mean([RR_before, RR_after], 'omitnan');
        if isnan(final_RR) || final_RR <= 0
            final_RR = baseline_RR;
        end

        % Compute how many peaks to insert in this gap
        num_insert = max(0, round(gap_len / final_RR) - 1);
        % Also limit by minimum peak distance in samples
        max_possible = floor(gap_len / (min_peak_dist_samples * dt)) - 1;
        if max_possible < 0
            num_insert = 0;
        else
            num_insert = min(num_insert, max_possible);
        end

        % Insert linearly spaced peaks in time, with linear interpolation of amplitude
        for j = 1:num_insert
            t_interp = t_left + j * (gap_len / (num_insert + 1));  % exact time spacing
            % linear interpolation of signal:
            val_interp = sig_left + (sig_right - sig_left) * (t_interp - t_left) / gap_len;

            new_times(end+1,1) = t_interp;
            new_vals(end+1,1) = val_interp;
        end
    end
end

% --- Combine original peaks and inserted peaks in time domain ---
% Original:
orig_times = TR_initial;                 % times of detected peaks
orig_vals  = ECG_recon(locs);            % values of detected peaks

if isempty(new_times)
    % No gaps filled
    TR = orig_times;
    R  = orig_vals;
else
    all_times = [orig_times; new_times];
    all_vals  = [orig_vals;  new_vals];
    % Sort by time:
    [TR_sorted, sort_idx] = sort(all_times);
    TR = TR_sorted;
    R  = all_vals(sort_idx);
end

data.times_after_gap = TR;
data.peaks_after_gap = R;

%% 7 - Q/S Peak Detection
% Parameters:
win_q_time = 0.2;  % 50 ms before R-peak
win_s_time = 0.3;  % 50 ms after R-peak

qs_plot_locs = [];  % to store final indices of Q or S

for i = 1:length(TR)
    t_r = TR(i);

    % --- Q-search: Find trough in [t_r - win_q_time, t_r)
    t_start_q = t_r - win_q_time;
    t_end_q = t_r;
    idx_q = find(Traw >= t_start_q & Traw < t_end_q);
    q_idx = NaN;
    if ~isempty(idx_q)
        [~, min_idx_rel] = min(ECG_recon(idx_q));
        q_idx = idx_q(min_idx_rel);
    end

    % --- S-search: Find trough in (t_r, t_r + win_s_time]
    t_start_s = t_r;
    t_end_s = t_r + win_s_time;
    idx_s = find(Traw > t_start_s & Traw <= t_end_s);
    s_idx = NaN;
    if ~isempty(idx_s)
        [~, min_idx_rel] = min(ECG_recon(idx_s));
        s_idx = idx_s(min_idx_rel);
    end

    % --- Choose Q or S peak (whichever is deeper trough)
    if ~isnan(q_idx) && ~isnan(s_idx)
        if ECG_recon(q_idx) < ECG_recon(s_idx)
            qs_plot_locs(end+1) = q_idx;
        else
            qs_plot_locs(end+1) = s_idx;
        end
    elseif ~isnan(q_idx)
        qs_plot_locs(end+1) = q_idx;
    elseif ~isnan(s_idx)
        qs_plot_locs(end+1) = s_idx;
    end
end

qs_plot_locs = unique(qs_plot_locs, 'stable');
TS = Traw(qs_plot_locs);
S = ECG_recon(qs_plot_locs);


%% 8. Quality Filtering
% % Filter R-peaks based on amplitude
% mR = mean(R);
% valid_R = (R >= 0.25*mR) & (R <= 2*mR);
% TR = TR(valid_R);
% R = R(valid_R);
% TS = TS(valid_R);  % Keep corresponding S-peaks
% S = S(valid_R);
% 
% % Filter S-peaks based on amplitude
% mS = mean(S);
% valid_S = (S <= 2*mS) & (S >= 0.25*mS);  % S-waves are negative
% TR = TR(valid_S);
% R = R(valid_S);
% TS = TS(valid_S);
% S = S(valid_S);

%% 9. Pair R and S peaks
% Sort and ensure chronological order
[TR, sortR] = sort(TR);
R = R(sortR);
[TS, sortS] = sort(TS);
S = S(sortS);


%% 10. Peak Correction
[TR, TS, R, S] = correct_ECG(patient, Traw, ECG_recon, TR, TS, R, S, plotmarkers);

%% 11. Calculate RR Intervals
T_RRint = TR(1:end-1);
RRint = diff(TR);

%% 12. Save Results
data.R = R;
data.TR = TR;
data.S = S;
data.TS = TS;
data.RRint = RRint;
data.T_RRint = T_RRint;

% Save temporary workspace
save(strcat('../WS/', filename), 'raw_data', 'data', '-append');

end