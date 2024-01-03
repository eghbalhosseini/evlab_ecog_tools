function [all_trial_timing,trial_based_frame] = extract_timing_for_ANNSet1(varargin)
%EXTRACT_TIMING_FOR_LANGLOC_AUDIO Summary of this function goes here
%   Detailed explanation goes here
% Define default values
TrigMat1 = [];
pre_start_ = 0;
exclusion_ = {};
all_trial_timing={};
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

start_exp_bit=bit1;
end_exp_bit=bit2;
start_audio_bit=bit3;
cond_bit=bit4;
repeat_bit=bit5;
end_audio_bit=bit6;
record_frame=find(start_exp_bit);

TrialTrig=[start_exp_bit,end_exp_bit,start_audio_bit,cond_bit,repeat_bit,end_audio_bit];
trial_keys={cond_bit,repeat_bit};
key_tags={'trial_cond','trial_repeat'};
trial_keys_idx=cellfun(@find,trial_keys,'uni',false);

expr_frame=(start_exp_bit & ~end_exp_bit);
expr_frame(1:pre_start_)=0;
for kk=1:length(exclusion_)
    exc=exclusion_{kk};
    exc_st=exc(1);
    exc_end=exc(2);
    expr_frame(exc_st:exc_end)=0;
end

trial_audio_st=find(diff(start_audio_bit.*expr_frame)==1)+1;
trial_audio_end=find(diff(end_audio_bit.*expr_frame)==1)+1;



max_trial_time=[0];
buffer_idx=50*sampling_freq;


for tr=1:length(trial_audio_st)
    start_idx=trial_audio_st(tr);
    end_idx=trial_audio_end(find(trial_audio_end>start_idx,1,'first'));
    if ~isempty(end_idx)
        tr_key_idx=cellfun(@(x) x((start_idx<=x ) & (x<=end_idx))',trial_keys_idx,'uni',false);
        if any(bit8(floor(start_idx):floor(end_idx))==1)
        else
            tr_key_idx=cellfun(@(x) x((start_idx<x ) & (x<end_idx))',trial_keys,'uni',false);
            % do it with substraction
            all_trial_timing{tr,1}=[start_idx,end_idx];
            all_trial_timing{tr,2}=tr_key_idx;
            max_trial_time=max(max_trial_time,end_idx);
        end

    end
end
trial_based_frame=[1:min(max(max_trial_time)+buffer_idx,length(TrigMat1))];

end

