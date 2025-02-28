function [events_table] = extract_behavioral_events_for_langloc_visual(varargin)
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
all_word_1=[];
all_probe_ans=[];
all_resp_ans=[];
all_react_time=[];
all_trial_comp=[];
all_block_id=[];
all_events_Table={};
for nn=1:size(d_events)
    event_table=readtable(strcat(d_events(nn).folder,'/',d_events(nn).name));
    tags=event_table.actual_onset*sampling_freq;
    
    condition=event_table.condition;
    non_fixation=find(~strcmp(condition,'F'));
    all_events_Table{nn,1}=[event_table(non_fixation,:)];
    word_1=event_table.word1;
    probe_ans=event_table.probe;
    resp_ans=event_table.response;
    reaction_time=event_table.RT;
    trial_complete=event_table.trial_completed;
    
    all_tags=[all_tags;tags(non_fixation)];
    all_condition=[all_condition;condition(non_fixation)];
    all_word_1=[all_word_1;word_1(non_fixation)];
    all_probe_ans=[all_probe_ans;probe_ans(non_fixation)];
    all_resp_ans=[all_resp_ans;resp_ans(non_fixation)];
    all_react_time=[all_react_time;reaction_time(non_fixation)];
    all_trial_comp=[all_trial_comp;trial_complete(non_fixation)];
    all_block_id=[all_block_id;trial_complete(non_fixation)*0+nn];
    
end


%events_table.list=cell2mat(cellfun(@(x) 0, events_table.list,'uni',false))

for tt=1:size(all_events_Table,1)
    event_table=all_events_Table{tt,1};
try 
    TJNew= removevars(event_table,{'date_time'});
    all_events_Table{tt}=TJNew;
catch
    all_events_Table{tt}=event_table;
end 
end 
events_table=[];
for nn=1:size(all_events_Table,1)
    events_table=[events_table;all_events_Table{nn}];
end 
%events_table.list=cell2mat(cellfun(@(x) 0, events_table.list,'uni',false))

mit_langloc_dir="/Users/eghbalhosseini/MyCodes/MIT_U01_experiments/U01_langloc_vJan2021/";
%d_file=dir(strcat(mit_langloc_dir,'/*.csv'));
%d_list=d_file(cellfun(@(x) ~isempty(x),regexpi({d_file(:).name},'test_list[1-9]+.csv')))
%d_list=arrayfun(@(x) fullfile(d_list(x).folder,d_list(x).name),1:size(d_list),'uni',false)

load(fullfile(mit_langloc_dir,'materials.mat'))
mat_r1r2r3=[materials.run1;materials.run2;materials.run3];

mat_r1r2r3=mat_r1r2r3(find(~ismember(mat_r1r2r3.condition,'F')),:);

event_keys = strcat(string(events_table.list), '_', string(events_table.trial));
mat_keys = strcat(string(mat_r1r2r3.list), '_', string(mat_r1r2r3.trial));

matching_rows = ismember(mat_keys, event_keys);
% Filter mat_r1r2r3 to keep only the matching rows
filtered_mat_r1r2r3 = mat_r1r2r3(matching_rows, :);

assert(all(ismember(events_table.word1,filtered_mat_r1r2r3.word1)));
assert(all(events_table.list==cell2mat(filtered_mat_r1r2r3.list)));
assert(all(events_table.trial==cell2mat(filtered_mat_r1r2r3.trial)));
words=arrayfun(@(x) ['word',num2str(x)],1:12,'uni',false);
words=[words,'probe'];
    
% replace them 
for word=words
        events_table.(word{1})=filtered_mat_r1r2r3.(word{1});
end 

events_table.final_list=events_table.list;
events_table.trial_onset=events_table.actual_onset;

end