function [events_table] = extract_behavioral_events_for_langloc_audio(varargin)
%EXTRACT_TIMING_FOR_LANGLOC_AUDIO Summary of this function goes here
%   Detailed explanation goes here

% Define default values
    d_events = [];


    % Parse name-value pairs
    for i = 1:2:length(varargin)
        switch varargin{i}
            case 'behavior_files'
                d_events = varargin{i+1};
            case 'sampling'
                sampling_freq = varargin{i+1};
            otherwise
                error('Unknown parameter name: %s', varargin{i});
        end
    end
all_tags=[];
all_condition=[];
all_strings=[];
all_audio_filenames=[];
all_probes=[];
all_probe_ans=[];
all_resp_ans=[];
all_react_time=[];
all_trial_comp=[];
all_block_id=[];
all_events_Table={};
for nn=1:size(d_events)
    event_name=strrep(d_events(nn).name,'._','');
    event_table=readtable(strcat(d_events(nn).folder,'/',event_name));
    tags=event_table.trial_onset*sampling_freq;
    all_events_Table{nn,1}=[event_table];
    condition=event_table.final_condition;
    non_fixation=find(~strcmp(condition,'F'));
    string=event_table.final_audio_transcript;
    audio_filename=event_table.final_audio_filename;
    probe_ans=event_table.final_probe;
    resp_ans=event_table.response;
    %reaction_time=event_table.RT;
    trial_complete=event_table.trial_completed;
    
    all_tags=[all_tags;tags(non_fixation)];
    all_condition=[all_condition;condition(non_fixation)];
    all_strings=[all_strings;string(non_fixation)];
    all_probe_ans=[all_probe_ans;probe_ans(non_fixation)];
    all_audio_filenames=[all_audio_filenames;audio_filename(non_fixation)];
    all_resp_ans=[all_resp_ans;resp_ans(non_fixation)];
    %all_react_time=[all_react_time;reaction_time(non_fixation)];
    all_trial_comp=[all_trial_comp;trial_complete(non_fixation)];
    all_block_id=[all_block_id;trial_complete(non_fixation)*0+nn];
    
end

events_table=all_events_Table{1};
for tt=1:size(all_events_Table,1)
    event_table=all_events_Table{tt,1};
try 
    TJNew= removevars(event_table,{'date_time'});
    all_events_Table{tt}=TJNew;
catch
    all_events_Table{tt}=event_table;
end 
end 
%events_table.list=cell2mat(cellfun(@(x) 0, events_table.list,'uni',false))
for nn=2:size(all_events_Table,1)
    events_table=[events_table;all_events_Table{nn}];
end 
