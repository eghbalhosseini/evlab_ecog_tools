function env = bandpass_envelopes(signal, signal_sr, ...
    band_in_Hz, filter_order, varargin)

% Measure envelopes from bandpass filters.
%
% -- Inputs --
%
% signal: [time x electrode] signal matrix
%
% signal_sr: sampling rate of the signal
%
% env_sr: sampling rate to use for the envelopes (typically much lower than the
% signal sampling rate, e.g. 100 Hz vs. 1200 Hz)
%
% band_in_Hz: 2-dimensional vector of frequency cutoffs
%
% filter_orders: order of the butterworth filters
%
% figure_directory: directory to save figures to
%
% -- Outputs --
%
% envelopes: [time x bands x electrode] matrix with envelopes
%
% 2016-1-26: Created by Sam NH
%
% 2016-09-23: Changes to handling of optional inputs and plotting, Sam NH
% 
% 2021-11-18: edits by Eghbal H. to fit Evlab pipeline

%% Setup
n_channels = size(signal,2);
p=inputParser();
addParameter(p, 'good_channels', 1:n_channels);
addParameter(p, 'filter_order', 6);
parse(p, varargin{:});
ops = p.Results;
%% Construct an FDESIGN object and call its BUTTER method.
h = fdesign.bandpass('N,F3dB1,F3dB2', filter_order, ...
    band_in_Hz(1), band_in_Hz(2), signal_sr);
Hd = design(h,'butter');
[B, A] = sos2tf(Hd.sosMatrix,Hd.scaleValues);
%% apply filter
fprintf('Applying filter...\n');
subb = filtfilt(B,A,signal);
%% measure envelope
fprintf('Calculating envelopes...\n');
env = abs(hilbert(subb));
%% truncate
env(env < 0) = 0;
end 

