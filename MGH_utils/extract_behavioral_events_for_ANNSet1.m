function [events_table] = extract_behavioral_events_for_ANNSet1(varargin)
%extract_behavioral_events_for_ANNSet1 Summary of this function goes here
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
final_list=[];
trial_num=[];
trial_onset=[];
audio_file_name={};
condition={};
sentence_txt={};
audio_end_time=[]; 
pressed_space_to_next=[];
completed_trial=[];
resume_number=[];
block_id=[];
all_events_Table={};

for nn=1:length(d_events)
    event_table=readtable(strcat(d_events(nn).folder,'/',d_events(nn).name));
    tags=event_table.trial_onset*sampling_freq;
    all_events_Table{nn,1}=[event_table];
    for nnn=1:height(event_table)
        block_id=[block_id;nn];
        final_list=[final_list;event_table.trial_onset(nnn)];
        trial_num=[trial_num;event_table.trial(nnn)];
        trial_onset=[trial_onset;event_table.trial_onset(nnn)];
        audio_file_name=[audio_file_name;event_table.final_audio_filename(nnn)];
        condition=[condition;event_table.final_condition(nnn)];
        sentence_txt=[sentence_txt;event_table.final_audio_transcript(nnn)];
        audio_end_time=[audio_end_time;event_table.audio_ended(nnn)];
        pressed_space_to_next=[pressed_space_to_next;event_table.pressed_space_to_continue(nnn)];
        completed_trial=[completed_trial;event_table.trial_completed(nnn)];
        resume_number=[resume_number;event_table.resume_number(nnn)];
    end 
end

events_table=all_events_Table{1};
for nn=2:size(all_events_Table,1)
    events_table=[events_table;all_events_Table{nn}];
end 
