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

all_events_Table={};
for nn=1:size(d_events)
    filename=strcat(d_events(nn).folder,'/',d_events(nn).name);
    opts = detectImportOptions(filename);
    opts = setvartype(opts,'char');
    event_table=readtable(filename,opts);
    condition=event_table.condition;
    non_fixation=find(~strcmp(condition,'F'));
    event_table=event_table(non_fixation,:);
    TJNew= removevars(event_table,{'date_time'});
    all_events_Table{nn,1}=[TJNew];
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
%vents_table.list=cell2mat(cellfun(@(x) 0, events_table.list,'uni',false))
for nn=2:size(all_events_Table,1)
    events_table=[events_table;all_events_Table{nn}];
end

events_table.list = cellfun(@str2num,events_table.list);
events_table.planned_onset = cellfun(@str2num,events_table.planned_onset);
events_table.actual_onset = cellfun(@str2num,events_table.actual_onset);
events_table.probe_answer = cellfun(@str2num,events_table.probe_answer);
events_table.response = cellfun(@str2num,events_table.response);
events_table.RT = cellfun(@str2num,events_table.RT);
events_table.trial = cellfun(@str2num,events_table.trial);
events_table.trial_completed = cellfun(@str2num,events_table.trial_completed);
events_table.final_list=events_table.list;

end