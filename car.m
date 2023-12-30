function [out_signal,signal_mean] = car(signal,chans_to_exclude,channelDim)

out_signal = zeros(size(signal));

all_chans = 1:size(signal,channelDim);
chans_selected = setxor(chans_to_exclude,all_chans);

fprintf(1, '> Common average filtering signal \n');

fprintf(1,'[');

if channelDim == 1
    signal_mean = mean(signal(chans_selected,:),channelDim);
elseif channelDim == 2
    signal_mean = mean(signal(:,chans_selected),channelDim);
end

for channelInd = 1:size(signal,channelDim)
    if channelDim == 1
        out_signal(channelInd,:) = signal(channelInd,:) - signal_mean;
    elseif channelDim == 2
        out_signal(:,channelInd) = signal(:,channelInd) - signal_mean;
    end
    fprintf(1,'.');
end


fprintf(1,'] done\n');
end