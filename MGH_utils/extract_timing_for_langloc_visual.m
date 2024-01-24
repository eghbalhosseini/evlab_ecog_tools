function [all_trial_timing,trial_based_frame] = extract_timing_for_langloc_visual(varargin)
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
bit1=TrigMat1(:,8);
bit2=TrigMat1(:,7);
bit3=TrigMat1(:,6);
bit4=TrigMat1(:,5);
bit5=TrigMat1(:,4);
bit6=TrigMat1(:,3);
bit7=TrigMat1(:,2);
bit8=TrigMat1(:,1);
% fixation 01101
fix_time=bit1 & ~bit2 & bit3 & bit4 & ~bit5;
% word1 00001
word1_time=bit1 & ~bit2 & ~bit3 & ~bit4 & ~bit5;
% word2 00010
word2_time=~bit1 & bit2 & ~bit3 & ~bit4 & ~bit5;
% word3 00011
word3_time=bit1 & bit2 & ~bit3 & ~bit4 & ~bit5;
% word4 00100
word4_time=~bit1 & ~bit2 & bit3 & ~bit4 & ~bit5;
% word5 00101
word5_time=bit1 & ~bit2 & bit3 & ~bit4 & ~bit5;
% word6 00110
word6_time=~bit1 & bit2 & bit3 & ~bit4 & ~bit5;
% word7 00111
word7_time=bit1 & bit2 & bit3 & ~bit4 & ~bit5;
% word8 01000
word8_time=~bit1 & ~bit2 & ~bit3 & bit4 & ~bit5;
% word9 01001
word9_time=bit1 & ~bit2 & ~bit3 & bit4 & ~bit5;
% word10 01010
word10_time=~bit1 & bit2 & ~bit3 & bit4 & ~bit5;
% word11 01011
word11_time=bit1 & bit2 & ~bit3 & bit4 & ~bit5;
% word12 01100
word12_time=~bit1 & ~bit2 & bit3 & bit4 & ~bit5;
% preprobe 01110
preprobe_time=~bit1 & bit2 & bit3 & bit4 & ~bit5;
% probe 10000
probe_time=~bit1 & ~bit2 & ~bit3 & ~bit4 & bit5;
% extra time 01111
extra_time=bit1 & bit2 & bit3 & bit4 & ~bit5;
SN_cond=bit6;
record_frame=bit8;
expr_frame=(record_frame & ~bit7);

expr_frame(1:pre_start_)=0;
for kk=1:length(exclusion_)
    exc=exclusion_{kk};
    exc_st=exc(1);
    exc_end=exc(2);
    expr_frame(exc_st:exc_end)=0;
end 

trial_keys={fix_time,...
    word1_time,...
    word2_time,...
    word3_time,...
    word4_time,...
    word5_time,...
    word6_time,...
    word7_time,...
    word8_time,...
    word9_time,...
    word10_time,...
    word11_time,...
    word12_time,...
    preprobe_time,...
    probe_time,...
    extra_time};
key_tags={'fix','word_1','word_2','word_3','word_4','word_5','word_6',...
    'word_7','word_8','word_9','word_10','word_11','word_12',...
    'preprobe','probe','extra_probe'};

trial_keys_idx=cellfun(@find,trial_keys,'uni',false);
% get start and end times for each trial 
trial_start=find(diff(fix_time)==1)+1;
trial_end=find(diff(extra_time)==-1)+1;
buffer_idx=10*sampling_freq;
trial_based_frame=[1:max(trial_end)+buffer_idx];
buffer_idx=50*sampling_freq;

max_trial_time=[0];
all_trial_timing={};
for tr=1:length(trial_start)
    start_idx=trial_start(tr);
    end_idx=trial_end(find(trial_end>start_idx,1,'first'));
    x=find(fix_time);
    tr_key_idx=cellfun(@(x) x((start_idx<x ) & (x<end_idx))',trial_keys_idx,'uni',false);
    % do it with substraction
    all_trial_timing{tr,1}=[start_idx,end_idx];
    all_trial_timing{tr,2}=tr_key_idx;
        max_trial_time=max(max_trial_time,end_idx);
end 

trial_based_frame=[1:min(max(max_trial_time)+buffer_idx,length(TrigMat1))];

end

