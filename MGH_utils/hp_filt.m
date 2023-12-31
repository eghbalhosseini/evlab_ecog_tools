function [signal,ops_out] = hp_filt(signal, sr, varargin)

% High-pass filter ECoG signal to remove very low frequencies (e.g. < 0.5 Hz)
% 
% Derived from Schalk lab code
% 
% signal: [time x electrode] singal matrix
% 
% sr: sampling rate of the signal
% 
% 2016-08-12: Created, Sam NH
% 
% 2018-01-21: Changed parameter handling
%
% 2021-11-18: changed parameter setting and how they get processed, ehoseini@mit.edu
% def
p=inputParser();
addParameter(p, 'order', 4);
addParameter(p, 'cutoff', 0.5);
parse(p, varargin{:});
ops = p.Results;
% design filter
Hd = design(fdesign.highpass('N,F3dB', ops.order, ops.cutoff, sr), 'butter');
% convert to pole form
[b, a] = sos2tf(Hd.sosMatrix,Hd.scaleValues);
% apply filter
signal = filtfilt(b,a,signal);

ops_out=ops;
ops_out.hp_filter_b=b;
ops_out.hp_filter_a=a;
end 