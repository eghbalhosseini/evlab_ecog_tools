function [good_channels,bad_channels,ops_out] = ...
    channel_selection_from_60Hz_noise(signal, sr, varargin)

% Uses 60 Hz noise to detect good/bad channels
% 
% Derived from Schalk lab code
% 
% signal: [time x electrode] singal matrix
% 
% sr: sampling rate of the signal
% 
% figure_fname: name of file to save plots to
% 
% 2016-08-12 - Created, Sam NH
% 
% 2016-09-20 - Modified to use the mode instead of the median to detect outliers
% 
% 2016-09-23 - Switched back to median, but use alternative z-scoring procedure
% to detect bad electrodes, Sam NH
% 
% 2019-06-12 - Added functionality to compute 60 Hz noise within arrays

% 2021-11-18: changed parameter setting and how they get processed, ehoseini@mit.edu


p=inputParser();
addParameter(p, 'bw', 0.6);
addParameter(p, 'frac', 0.2); % fraction of electrodes used to calculate standard deviation
addParameter(p, 'thresh', 5);% threshold needed to calculate
addParameter(p, 'array_inds', ones(1,size(signal,2)));% indices for different arrays
addParameter(p, 'min_nchannels', 10); % minimum number of channels needed to calculate distribution
addParameter(p, 'chnames', {});
parse(p, varargin{:});
ops = p.Results;

n_channels = size(signal,2);

noNaN_channels = find(all(~isnan(signal)));
signal_noNaN = signal(:, noNaN_channels);

% measure 60 Hz power in rms units
[b,a] = iirpeak(60/(sr/2), ops.bw/(sr/2));
rms60Hz = mean(sqrt(filter(b, a, signal_noNaN).^2),1);

% infer mode
% [N,bin_centers] = hist(noise60Hz_rms, 1000);
% [~,xi] = max(N);
% mode_60Hz = bin_centers(xi);

unique_arrays = unique(ops.array_inds);
rms60Hz_zscore = nan(1, n_channels);
for i = 1:length(unique_arrays)
    
    array_chan = find(ops.array_inds(noNaN_channels)==unique_arrays(i));
    
    % fraction of channels to use to calculate z-score
    frac_corresponding_to_min = ops.min_nchannels / length(array_chan);
    frac = min(max(frac_corresponding_to_min, ops.frac), 1);
    
    % zsore using samples nearest the median
    if frac < 1
        rms60Hz_zcore_within_array = zscore_using_central_samples(rms60Hz(array_chan), frac);
        rms60Hz_zscore(noNaN_channels(array_chan)) = rms60Hz_zcore_within_array;
    end
end

% select good channels
bad_channels = (abs(rms60Hz_zscore) > ops.thresh) | ~noNaN_channels;
good_channels = find(~bad_channels);
bad_channels=setdiff(1:n_channels, good_channels);

ops_out=ops;
ops_out.rms60Hz=rms60Hz;

end 