function [all_trial_timing,trial_based_frame,exp_seq] = extract_timing_for_LangMathMusic(varargin)
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
%% predefined seqeuences 
conditions={'sent','nonw','math','music','sent_repeat'};
exp_seq=struct;
exp_seq.awake_ids=[1,2,1,3,1,4,1,3,1,2,1];
exp_seq.awake_cond=arrayfun(@(x) conditions{x}, exp_seq.awake_ids,'uni',false);

exp_seq.under_ids=[1 1 2 1 1 4 1 1 2 5 5];
exp_seq.under_cond=arrayfun(@(x) conditions{x}, exp_seq.under_ids,'uni',false);

exp_seq.awake_id_blocks=repmat(exp_seq.awake_ids,1,4);
exp_seq.awake_cond_blocks=repmat(exp_seq.awake_cond,1,4);


exp_seq.under_id_blocks=repmat(exp_seq.under_ids,1,4);
exp_seq.under_cond_blocks=repmat(exp_seq.under_cond,1,4);
repmat({'awake'},1, length(exp_seq.awake_cond_blocks))

exp_seq.trial_id=horzcat(exp_seq.awake_id_blocks,exp_seq.under_id_blocks);
exp_seq.trial_cond=horzcat(exp_seq.awake_cond_blocks,exp_seq.under_cond_blocks);
exp_seq.block_cond=horzcat(repmat({'awake'},1, length(exp_seq.awake_cond_blocks)), ...
    repmat({'under'},1, length(exp_seq.under_cond_blocks)));

aud_start_bit=bit1;
question_str_bit=bit2;
trial_end_bit=bit3;
expr_frame=aud_start_bit*0+1;
expr_frame(1:pre_start_)=0;
for kk=1:length(exclusion_)
    exc=exclusion_{kk};
    exc_st=exc(1);
    exc_end=exc(2);
    expr_frame(exc_st:exc_end)=0;
end 
audio_on=find(diff(aud_start_bit.*expr_frame)==+1)+1;
audio_end=find(diff(aud_start_bit.*expr_frame)==-1)+1;
trial_start_temp=find(diff(aud_start_bit.*expr_frame)==+1)+1;
trial_start=trial_start_temp-.5*sampling_freq;
trial_end=find(diff(trial_end_bit.*expr_frame)==+1)+1;
question_on=find(diff(question_str_bit.*expr_frame)==+1)+1;
question_end=find(diff(question_str_bit.*expr_frame)==-1)+1;




trial_keys={audio_on,audio_end,question_on,question_end};
key_tags={'audio_start','audio_end','question_start','question_end'};
buffer_idx=50*sampling_freq;

max_trial_time=[0];
all_trial_timing={};
for tr=1:length(trial_start)
    start_idx=trial_start(tr);
    %block_cond=exp_seq.block_cond{tr};
    %trial_cond=exp_seq.trial_cond{tr};
    end_idx=trial_end(find(trial_end>start_idx,1,'first'));

    tr_key_idx=cellfun(@(x) x((start_idx<x ) & (x<end_idx))',trial_keys,'uni',false);
    % do it with substraction
    all_trial_timing{tr,1}=[start_idx,end_idx];
    all_trial_timing{tr,2}=tr_key_idx;
    max_trial_time=max(max_trial_time,end_idx);
end 

trial_based_frame=[1:min(max(max_trial_time)+buffer_idx,length(TrigMat1))];

end

