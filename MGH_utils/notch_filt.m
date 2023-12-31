function [signal,ops_out] = notch_filt(signal, sr, varargin)

% Uses 60 Hz noise to detect good/bad channels
% 
% Derived from Schalk lab code
% 
% signal: [time x electrode] singal matrix
% 
% sr: sampling rate of the signal
% 
% 2016-08-12: Created, Sam NH
% 
% 2019-01-21: Altered parameter handling, Sam NH
%
% 2021-11-19: updated the code with evlab setting, ehoseini@mit.edu

p=inputParser();
addParameter(p, 'bw', 1);
addParameter(p, 'freqs', [60, 120, 180, 240]);
parse(p, varargin{:});
ops = p.Results;

n_channels = size(signal,2);

% notch filter parameters
b = cell(1,length(ops.freqs)); 
a = cell(1,length(ops.freqs));
for i = 1:length(ops.freqs)
    [b{i},a{i}] = iirnotch(ops.freqs(i) / (sr/2), ops.bw / (sr/2));
end

% fvtool(b{2}, a{2}, 'Fs', ecog_sr);
% figure;
% [h,t] = impz(b{1},a{1},[],ecog_sr);
% plot(t,h);
% xlim([0 1])

% apply notch filter
pbar=ProgressBar(n_channels);
for i = 1:n_channels
    for j = 1:length(ops.freqs)
        signal(:,i) = filtfilt(...
            b{j},a{j},signal(:,i));
    end
    pbar.step([],[],[]);
end 

ops_out=ops;
ops_out.notch_filter_b=[b];
ops_out.notch_filter_a=[a];
end 

