function shift_samples_per_trial=detect_wavelet_and_pick_index(mic_signal,trial_timing,sampling_freq)
%% define the wavelet 
Fs = 48000; % Sampling rate in Hz
t_interval = 0.2; % Time interval in seconds (200 ms)

% Calculate the number of samples
num_samples = round(Fs * t_interval);

% Create the time vector
x = (0:num_samples-1) / Fs;

y_wavelet=.5*sin(2*pi*40*x).*exp(-.5*(x-.1).^2/(0.025^2));
y_wavelet_rs=resample(y_wavelet,sampling_freq,Fs);
%% define a low pass filter 
cutoff_freq = 60;      % Low-pass cutoff frequency in Hz (adjust as needed)
nyquist = sampling_freq / 2;       % Nyquist frequency
filter_order = 4;       % 4th-order Butterworth
[b, a] = butter(filter_order, cutoff_freq / nyquist, 'low');

%% clean up the signal and wavelet; 
filt_orig_signal = filtfilt(b, a, mic_signal);
filt_wavelet     = filtfilt(b, a, y_wavelet_rs);
% make their hilbert version
filt_orig_signal_line=abs(hilbert(zscore(filt_orig_signal)));
filt_wavelet_rs_line=abs(hilbert(zscore(filt_wavelet)));
shift_samples_per_trial=[];
for trial_id=1:length(trial_timing)
    trial_key_index=trial_timing{trial_id,2};
    audio_start=trial_key_index{1};
    audio_end=trial_key_index{2};
    filt_orig_trial=filt_orig_signal(:,floor([(audio_start+sampling_freq*-.0):(audio_end+sampling_freq*+1.5)]));
    %filt_orig_trial_line=filt_orig_signal_line(:,floor([(trial_timing{trial_id}(1)+sampling_freq*-.5):(trial_timing{trial_id}(1)+sampling_freq*+1.5)]));
    
    threshold_value = 0.075;  % For instance; you may adjust after experimentation
    [foundWavelet, lag, peakCorr,shiftSamples] = checkWaveletPresence(filt_orig_trial, filt_wavelet, threshold_value);
    if not(foundWavelet)
        %figure;
        %plot(filt_orig_trial)
        %keyboard;
        shift_samples_per_trial=[shift_samples_per_trial;0];
    else
        shift_samples_per_trial=[shift_samples_per_trial;shiftSamples];
    end 
    if foundWavelet
        fprintf('Wavelet found! Max corr = %.3f at lag %d.\n', peakCorr, lag);
    else
        fprintf('Wavelet not detected (max corr = %.3f < threshold).\n', peakCorr);
    end

end 
end 


function [isPresent, bestLag, corrPeak,shiftSamples] = checkWaveletPresence(signal, wavelet, threshold)
% checkWaveletPresence
% Performs cross-correlation between 'signal' and 'wavelet', then checks 
% if the peak correlation exceeds 'threshold'.
%
% Inputs:
%   signal   : 1D array (the data in which we search for the wavelet)
%   wavelet  : 1D array (the wavelet "template" we want to detect)
%   threshold: scalar (minimum normalized correlation to declare "present")
%
% Outputs:
%   isPresent : logical (true if wavelet is detected, false otherwise)
%   bestLag   : sample lag where cross-correlation is highest
%   corrPeak  : peak normalized cross-correlation value

    % ------------------------
    % 1. Cross-correlation
    % ------------------------
    [r, lags] = xcorr(signal, wavelet);
    
    % We often want to look at absolute correlation in case 
    % the wavelet might be inverted in the signal
    [maxVal, idx] = max(abs(r));
    bestLag = lags(idx);

    % ------------------------
    % 2. Normalization (optional but recommended)
    % ------------------------
    % Normalized cross-correlation: correlation / (||signal|| * ||wavelet||)
    norm_factor = norm(signal) * norm(wavelet);
    corrPeak = maxVal / norm_factor;

    % ------------------------
    % 3. Compare Against Threshold
    % ------------------------
    isPresent = (corrPeak >= threshold);

       %bestIndex = bestLag + length(wavelet);
% Correct the bestIndex calculation
    bestIndex = bestLag + 1;  % Add 1 to convert lag into 1-based indexing

% Compute the number of samples needed to shift
    shiftSamples = 1 - bestIndex;  % If bestIndex == 1, shiftSamples == 0 (no shift needed)


    % --- 3) Compute shiftSamples so wavelet begins at index 1 ---
    % If wavelet is found at 'bestIndex', to move wavelet's start to index 1, 
    % the shift = 1 - bestIndex. (Negative = shift left, positive = shift right.)
    %shiftSamples = 1 - bestIndex;
end