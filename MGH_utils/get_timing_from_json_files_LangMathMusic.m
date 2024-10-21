function [trial_timing_dec] = get_timing_from_json_files_LangLocAudio(dataout,all_trial_timing,events_table,audio_align_path,with_wavelet)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
all_timing_diff=[];
initial_table_keys={'trial_start_end',};
key_tags={'audio_start','audio_end','question_start','question_end'};
all_condition={};
all_audio_length=[];
%
trial_timing_dec={};
% 
% Specify the file name
filename = '/Users/eghbalhosseini/MyCodes/Task_MIT_LangMathMusic/IntraOp_stim_20200731.xlsx';

% Get the sheet names
Stim_table=struct;
[~, sheetNames] = xlsfinfo(filename);
sheetNames=sheetNames(1:4);

% Loop through each sheet and read the data
for i = 1:numel(sheetNames)
    sheetName = sheetNames{i};
    data = readtable(filename, 'Sheet', sheetName);
    Stim_table.(sheetName)=data;
    % Display or process the data

end


if with_wavelet
    awake_audio_path='/Users/eghbalhosseini/MyCodes/Task_MIT_LangMathMusic/stimuli/awake/';
    under_audio_path='/Users/eghbalhosseini/MyCodes/Task_MIT_LangMathMusic/stimuli/going_under/';
    delay_val=0.0;
else
    awake_audio_path='/Users/eghbalhosseini/MyCodes/Task_MIT_LangMathMusic/stimuli/awake/';
    under_audio_path='/Users/eghbalhosseini/MyCodes/Task_MIT_LangMathMusic/stimuli/going_under/';
    delay_val=0.0;
end 
all_audio_filenames=events_table.fname;
all_phases=events_table.phase;
    
for k=1:size(all_trial_timing,1)
    table_keys=initial_table_keys;
    trial_trig_idx_dec=cell(size(initial_table_keys));
    trial_trig_start=all_trial_timing{k,1}(1);
    trial_trig_end=all_trial_timing{k,1}(2);
    trial_trig_keys=all_trial_timing{k,2};
    [P,Q] = rat(dataout.ops.fsDownsample/dataout.ops.sr);
    trial_trig_start_dec=floor(trial_trig_start*P/Q);
    trial_trig_end_dec=floor(trial_trig_end*P/Q);
    trial_trig_keys_dec=cellfun(@(x) x*P/Q,trial_trig_keys,'UniformOutput',false);
    trial_trig_keys_dec=cell2mat(trial_trig_keys_dec);
    trial_sent_file=all_audio_filenames{k};
    trial_sent_file=erase(trial_sent_file,'.wav');
    trial_id=extract(trial_sent_file,digitsPattern);
    trial_type=extract(trial_sent_file,lettersPattern);
    phase=all_phases{k};

    table_keys=horzcat(table_keys,{'fix','audio_onset','pause_onset','question_onset','response_onset'});
    table_words={};
    if strcmp(trial_type{1},'sent')
        lookup_table=Stim_table.sentences;
        lookup_table=lookup_table(ismember(lookup_table.Component,phase),:);
        lookup_row=lookup_table(ismember(erase(lookup_table.fname,"'"),sprintf('sentence_%s.wav',trial_id{1})),:);
        assert(size(lookup_row,1)==1);
        transcript =lookup_row.Sentence{1};
        question = lookup_row.Question;
        qfname = erase(lookup_row.qfname,"'");
        % add timings 
        % {'audio_start','audio_end','question_start','question_end'};
        trial_trig_idx_dec{ismember(table_keys,'fix')}=[trial_trig_start_dec,trial_trig_keys_dec(1)];
        table_words{ismember(table_keys,'fix')}='';
        trial_trig_idx_dec{ismember(table_keys,'audio_onset')}=[trial_trig_keys_dec(1),trial_trig_keys_dec(2)];
        table_words{ismember(table_keys,'audio_onset')}=transcript;
        trial_trig_idx_dec{ismember(table_keys,'pause_onset')}=[trial_trig_keys_dec(2),trial_trig_keys_dec(3)];
        table_words{ismember(table_keys,'pause_onset')}='';
        trial_trig_idx_dec{ismember(table_keys,'question_onset')}=[trial_trig_keys_dec(3),trial_trig_keys_dec(4)];
        table_words{ismember(table_keys,'question_onset')}=question{1};
        trial_trig_idx_dec{ismember(table_keys,'response_onset')}=[trial_trig_keys_dec(4),trial_trig_end_dec];
        table_words{ismember(table_keys,'response_onset')}='';
    elseif strcmp(trial_type{1},'nonword')
        lookup_table=Stim_table.nonword_sequences;
        lookup_table=lookup_table(ismember(lookup_table.Component,phase),:);
        lookup_row=lookup_table(ismember(erase(lookup_table.fname,"'"),sprintf('nonword_%s.wav',trial_id{1})),:);
        assert(size(lookup_row,1)==1);
        nonword_keys=lookup_row.Properties.VariableNames(contains(lookup_row.Properties.VariableNames,'nonword_'));
        transcript=strjoin(cellfun(@(x)lookup_row.(x), nonword_keys),' ');        
        question = lookup_row.Question;
        % add timings 
        % {'audio_start','audio_end','question_start','question_end'};
        trial_trig_idx_dec{ismember(table_keys,'fix')}=[trial_trig_start_dec,trial_trig_keys_dec(1)];
        table_words{ismember(table_keys,'fix')}='';
        trial_trig_idx_dec{ismember(table_keys,'audio_onset')}=[trial_trig_keys_dec(1),trial_trig_keys_dec(2)];
        table_words{ismember(table_keys,'audio_onset')}=transcript;
        trial_trig_idx_dec{ismember(table_keys,'pause_onset')}=[trial_trig_keys_dec(2),trial_trig_end_dec];
        table_words{ismember(table_keys,'pause_onset')}='';
        trial_trig_idx_dec{ismember(table_keys,'question_onset')}=[trial_trig_keys_dec(2),trial_trig_end_dec];
        table_words{ismember(table_keys,'question_onset')}=question{1};
        trial_trig_idx_dec{ismember(table_keys,'response_onset')}=[trial_trig_keys_dec(2),trial_trig_end_dec];
        table_words{ismember(table_keys,'response_onset')}='';


            
    elseif strcmp(trial_type{1},'math')
        lookup_table=Stim_table.math;
        lookup_table=lookup_table(ismember(lookup_table.Var1,phase),:);
        lookup_row=lookup_table(ismember(erase(lookup_table.fname,"'"),sprintf('math_%s.wav',trial_id{1})),:);
        assert(size(lookup_row,1)==1);
        transcript =lookup_row.Var3{1};
        question = '';
        % add timings 
        % {'audio_start','audio_end','question_start','question_end'};
        trial_trig_idx_dec{ismember(table_keys,'fix')}=[trial_trig_start_dec,trial_trig_keys_dec(1)];
        table_words{ismember(table_keys,'fix')}='';
        trial_trig_idx_dec{ismember(table_keys,'audio_onset')}=[trial_trig_keys_dec(1),trial_trig_keys_dec(2)];
        table_words{ismember(table_keys,'audio_onset')}=transcript;
        trial_trig_idx_dec{ismember(table_keys,'pause_onset')}=[trial_trig_keys_dec(2),trial_trig_end_dec];
        table_words{ismember(table_keys,'pause_onset')}='';
        trial_trig_idx_dec{ismember(table_keys,'question_onset')}=[trial_trig_keys_dec(2),trial_trig_end_dec];
        table_words{ismember(table_keys,'question_onset')}=question;
        trial_trig_idx_dec{ismember(table_keys,'response_onset')}=[trial_trig_keys_dec(2),trial_trig_end_dec];
        table_words{ismember(table_keys,'response_onset')}='';
    elseif strcmp(trial_type{1},'music')
        lookup_table=Stim_table.music;
        lookup_table=lookup_table(ismember(lookup_table.Var1,phase),:);
        lookup_row=lookup_table(ismember(erase(lookup_table.fname,"'"),sprintf('music_%s.wav',trial_id{1})),:);
        assert(size(lookup_row,1)==1);
        transcript =lookup_row.originalPiece{1};
        question = '';
        % add timings 
        % {'audio_start','audio_end','question_start','question_end'};
        trial_trig_idx_dec{ismember(table_keys,'fix')}=[trial_trig_start_dec,trial_trig_keys_dec(1)];
        table_words{ismember(table_keys,'fix')}='';
        trial_trig_idx_dec{ismember(table_keys,'audio_onset')}=[trial_trig_keys_dec(1),trial_trig_keys_dec(2)];
        table_words{ismember(table_keys,'audio_onset')}=transcript;
        trial_trig_idx_dec{ismember(table_keys,'pause_onset')}=[trial_trig_keys_dec(2),trial_trig_end_dec];
        table_words{ismember(table_keys,'pause_onset')}='';
        trial_trig_idx_dec{ismember(table_keys,'question_onset')}=[trial_trig_keys_dec(2),trial_trig_end_dec];
        table_words{ismember(table_keys,'question_onset')}=question;
        trial_trig_idx_dec{ismember(table_keys,'response_onset')}=[trial_trig_keys_dec(2),trial_trig_end_dec];
        table_words{ismember(table_keys,'response_onset')}='';
    else
        error('trial_type is unkown!')
    end
    if strcmp(phase,'awake')
        audio_path=awake_audio_path;
    else
        audio_path=under_audio_path;
    end 
    trial_trig_idx_dec{ismember(table_keys,sprintf('trial_start_end'))}=[trial_trig_start_dec,trial_trig_end_dec];
    table_words{ismember(table_keys,'trial_start_end')}=transcript;
    pairs_times=cell2mat(trial_trig_idx_dec');
    durations=pairs_times(:,2)-pairs_times(:,1);
    durations_ms=durations*Q/P;
    
    audio_fname=strcat(audio_path,filesep,all_audio_filenames{k});
    [orig_audio,orig_fs]=audioread(audio_fname);
    %soundsc(orig_audio,orig_fs)
    orig_aud_ts=length(orig_audio)/orig_fs;
    all_audio_length=[all_audio_length;orig_aud_ts];
    trial_aud_ts=diff([all_trial_timing{k,2}{1:2}])./dataout.MetaTags.SamplingFreq;
    fprintf('timing difference %f ms\n',(orig_aud_ts-trial_aud_ts)*1000)
    all_timing_diff=[all_timing_diff;(orig_aud_ts-trial_aud_ts)*1000];

    % construct fix index
    
    % construct the rest of tags
    
    trial_timing_dec{k,1}=table(table_keys',table_words',pairs_times(:,1),pairs_times(:,2),durations,durations_ms,'VariableNames',{'key','string','start','end','duration','duration_ms'});
    trial_timing_dec{k,2}=trial_type{1};
    trial_timing_dec{k,3}=strrep(phase,' ','_');
end

figure;
histogram(all_timing_diff)
xlabel('ms')
assert(max(abs(all_timing_diff))<50)

end