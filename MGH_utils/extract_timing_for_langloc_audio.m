function [all_trial_timing,trial_based_frame] = extract_timing_for_langloc_audio(varargin)
%EXTRACT_TIMING_FOR_LANGLOC_AUDIO Summary of this function goes here
%   Detailed explanation goes here

% Define default values
    TrigMat1 = [];
    pre_start_ = 0;
    exclusion_ = {};

    % Parse name-value pairs
    for i = 1:2:length(varargin)
        switch varargin{i}
            case 'trigger'
                TrigMat1 = varargin{i+1};
            case 'sampling'
                sampling_freq = varargin{i+1};
            case 'pre_start'
                pre_start_ = varargin{i+1};
            case 'exclusion'
                exclusion_ = varargin{i+1};
            otherwise
                error('Unknown parameter name: %s', varargin{i});
        end
    end
bit1=TrigMat1(:,1);
bit2=TrigMat1(:,2);
bit3=TrigMat1(:,3);
bit4=TrigMat1(:,4);
bit5=TrigMat1(:,5);
bit6=TrigMat1(:,6);
bit7=TrigMat1(:,7);
bit8=TrigMat1(:,8);

start_expt_bit=bit1;
end_expt_bit=bit2;
start_aud_bit=bit3;
cond_bit=bit4;
end_aud_bit=bit5;
probe_bit=bit6;
fixation_bit=bit7;
expr_frame=(start_expt_bit & ~end_expt_bit);
expr_frame(1:pre_start_)=0;
for kk=1:length(exclusion_)
    exc=exclusion_{kk};
    exc_st=exc(1);
    exc_end=exc(2);
    expr_frame(exc_st:exc_end)=0;
end 
fixation_on=find(diff(fixation_bit.*expr_frame)==+1)+1;
fixation_off=find(diff(fixation_bit.*expr_frame)==-1)+1;
trial_start_temp=find(diff(start_aud_bit.*expr_frame)==+1)+1;
trial_start=trial_start_temp-.2*sampling_freq;
trial_end_temp=find(diff(probe_bit.*expr_frame)==-1)+1;
trial_end=trial_end_temp+.2*sampling_freq;

audio_end=find(diff(start_aud_bit.*expr_frame)==-1)+1;
audio_start=find(diff(start_aud_bit.*expr_frame)==+1)+1;
pre_probe_start=find(diff(end_aud_bit.*expr_frame)==-1)+1;
probe_start=find(diff(probe_bit.*expr_frame)==1)+1;
probe_end=find(diff(probe_bit.*expr_frame)==-1)+1;
post_probe_start=find(diff(probe_bit.*expr_frame)==-1)+1;


trial_keys={audio_start,audio_end,pre_probe_start,probe_start,post_probe_start};
key_tags={'audio_start','audio_end','pre_probe','probe','post_probe'};
buffer_idx=50*sampling_freq;
trial_based_frame=[1:max(trial_end)+buffer_idx];

all_trial_timing={};
for tr=1:length(trial_start)
    start_idx=trial_start(tr);
    end_idx=trial_end(find(trial_end>start_idx,1,'first'));

    tr_key_idx=cellfun(@(x) x((start_idx<x ) & (x<end_idx))',trial_keys,'uni',false);
    % do it with substraction
    all_trial_timing{tr,1}=[start_idx,end_idx];
    all_trial_timing{tr,2}=tr_key_idx;
end 


end

