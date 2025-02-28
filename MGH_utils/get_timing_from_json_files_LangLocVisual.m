function [trial_timing_dec] = get_timing_from_json_files_LangLocVisual(dataout,all_trial_timing,events_table,audio_align_path,with_wavelet)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
stim_table=events_table(find(~ismember(events_table.condition,'F')),:);
%% save data as an object for ease of further processing 
trial_timing_dec={};
trial_timing={};
key_tags={'fix','word_1','word_2','word_3','word_4','word_5','word_6',...
    'word_7','word_8','word_9','word_10','word_11','word_12',...
    'preprobe','probe','extra_probe'};
words=arrayfun(@(x) ['word',num2str(x)],1:12,'uni',false);
words=[words,'probe'];

all_condition=stim_table.condition;
for k=1:size(all_trial_timing,1)        
        trial_trig_start=all_trial_timing{k,1}(1);
        trial_trig_end=all_trial_timing{k,1}(2);
        [P,Q] = rat(dataout.ops.fsDownsample/dataout.ops.sr);
        trial_type=all_condition{k};
        
        stim_row=stim_table(k,:);
        assert(stim_row.condition{:}==trial_type)
        
        table_words=cell(size(key_tags));
        table_words(:)={''};
        corres=cellfun(@(x) find(ismember(erase(key_tags,'_'),x)), words,'uni',false);
        for idx=1:length(words)
            table_words(corres{idx})=stim_row.(words{idx});
        end 
        table_words{1}='';% fixation string 
        trial_trig_start_dec=floor(trial_trig_start*P/Q);
        trial_trig_end_dec=floor(trial_trig_end*P/Q);
        trial_idx_dec=trial_trig_start_dec:trial_trig_end_dec;
        trial_timing_dec{k,1}=key_tags';
        trial_timing{k,1}=key_tags';
        % 
        trial_trig_idx=all_trial_timing{k,2};
        trial_trig_idx_dec={};
        trial_trig={};
        for key=1:length(key_tags)
            state_idx=trial_trig_idx{key};
            state_trig_start_dec=floor(min(state_idx)*P/Q);
            state_trig_end_dec=floor(max(state_idx)*P/Q);
            trial_trig_idx_dec{key}=[state_trig_start_dec,state_trig_end_dec];
            trial_trig{key}=[min(state_idx),max(state_idx)];
            table_words{key}=table_words{key};
        end                                                              
        %trial_timing_dec{k,2}=trial_trig_idx_dec';
        trial_timing{k,2}=trial_trig;
        pairs_times=cell2mat(trial_trig_idx_dec');
        
        trial_timing_dec{k,1}=table(key_tags',table_words',pairs_times(:,1),pairs_times(:,2),'VariableNames',{'key','string','start','end'});
        trial_timing_dec{k,2}=trial_type;
end 

% figure;
% histogram(all_timing_diff)
% xlabel('ms')
% assert(max(abs(all_timing_diff))<50)

end