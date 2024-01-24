function [trial_timing_dec] = get_timing_from_json_files_LangLocAudio(dataout,all_trial_timing,events_table,audio_align_path,with_wavelet)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
all_timing_diff=[];
table_keys={'fix','word_1','word_2','word_3','word_4','word_5','word_6',...
    'word_7','word_8','word_9','word_10','word_11','word_12',...
    'preprobe','probe','extra_probe'};
all_condition={};
%
trial_timing_dec={};
if with_wavelet
    audio_path='/Users/eghbalhosseini/MyCodes/MIT_U01_experiments/U01_langloc_audio_vFeb2021/stimuli_w_wavelet';
    delay_val=0.2;
else
    audio_path='/Users/eghbalhosseini/MyCodes/MIT_U01_experiments/U01_langloc_audio_vFeb2021/stimuli_orig';
    delay_val=0;
end 
all_audio_filenames=events_table.final_audio_filename;
for k=1:size(all_trial_timing,1)    
trial_trig_idx_dec=cell(size(table_keys));
trial_words_idx=cell(size(table_keys));

trial_trig_start=all_trial_timing{k,1}(1);
trial_trig_end=all_trial_timing{k,1}(2);
trial_audio_start=all_trial_timing{k,2}{1};
trial_audio_end=all_trial_timing{k,2}{2};
trial_preprobe_start=all_trial_timing{k,2}{3};
trial_probe_start=all_trial_timing{k,2}{4};
trial_post_probe_start=all_trial_timing{k,2}{5};
[P,Q] = rat(dataout.ops.fsDownsample/dataout.ops.sr);
trial_trig_start_dec=floor(trial_trig_start*P/Q);
trial_trig_end_dec=floor(trial_trig_end*P/Q);
trial_auido_start_dec=floor(trial_audio_start*P/Q);
trial_auido_end_dec=floor(trial_audio_end*P/Q);
trial_preprobe_start_dec=floor(trial_preprobe_start*P/Q);
trial_probe_start_dec=floor(trial_probe_start*P/Q);
trial_post_probe_start_dec=floor(trial_post_probe_start*P/Q);

trial_sent_file=all_audio_filenames{k};
trial_sent_file=erase(trial_sent_file,'.wav');
trial_sent_split=strsplit(trial_sent_file,'_');
sent_id=str2num(trial_sent_split{1});
sent_type=(trial_sent_split{2});
    if strcmp(sent_type,'English')
        trial_type='S';
    elseif strcmp(sent_type,'Nonsense')
        trial_type='N';
    else
        error('trial_type is unkown!')
    end
all_condition{k}=trial_type;
trial_audio_align=strrep(all_audio_filenames{k},'.wav','_handfix.json');

fname=strcat(audio_align_path,filesep,trial_audio_align);
audio_align = jsondecode(fileread(fname));
word_align=audio_align.words;
assert(size(word_align,1)==12)
if not(size(word_align,1)==12)
    disp(fname)
end
if ~isa(word_align,'struct')
    field_names={'case','endOffset','startOffset','word','start','end'};
    tempo=cellfun(@(x) cellfun(@(y) x.(y), field_names,'uni',false),word_align,'uni',false);
    word_align_temp=struct;
    for nn=1:size(tempo,1)
        temp=tempo{nn};
        for nnn=1:length(field_names)
            word_align_temp(nn).(field_names{nnn})=temp{nnn};
        end
    end
    word_align=word_align_temp;
end

for kk = 1:length(word_align)
    word_align(kk).start = word_align(kk).start + delay_val;
    word_align(kk).end = word_align(kk).end + delay_val;
end

audio_fname=strcat(audio_path,filesep,trial_sent_file,'.wav');
[orig_audio,orig_fs]=audioread(audio_fname);
%soundsc(orig_audio,orig_fs)
orig_aud_ts=length(orig_audio)/orig_fs;
trial_aud_ts=diff([all_trial_timing{k,2}{1:2}])./dataout.MetaTags.SamplingFreq;
fprintf('timing difference %f ms\n',(orig_aud_ts-trial_aud_ts)*1000)
[P_aud,Q_aud] = rat(dataout.ops.fsDownsample/orig_fs);
all_timing_diff=[all_timing_diff,(orig_aud_ts-trial_aud_ts)*1000];

for p=1:size(word_align,1)
    
    if ~strcmp(word_align(p).case,'success')
        fprintf('no allignment for word %d in %s \n',p,trial_sent_file)
        
    end
    word_st_ts=word_align(p).start;
    if word_st_ts==0
        word_st_ts=eps;
    end
    word_end_ts=word_align(p).end;
    words_audio{p}=orig_audio(ceil(word_st_ts*orig_fs):floor(word_end_ts*orig_fs));
    %soundsc(words_audio{p},orig_fs)
    %pause(length(words_audio{p})/orig_fs+.5);
    %words{1,p}=word_align(p).alignedWord;
    try
        words{1,p}=word_align(p).word;
    catch err
        words{1,p}=word_align(p).alignedWord;
    end
    word_idx_dec=[ceil(word_st_ts*dataout.ops.fsDownsample),floor(word_end_ts*dataout.ops.fsDownsample)];
    trial_trig_idx_dec{find(ismember(table_keys,sprintf('word_%d',p)))}=[trial_auido_start_dec]+word_idx_dec;
    trial_words_idx{find(ismember(table_keys,sprintf('word_%d',p)))}=words{1,p};
end
% construct fix index
trial_trig_idx_dec{find(ismember(table_keys,sprintf('fix',p)))}=[trial_trig_start_dec,trial_auido_start_dec-1];
trial_words_idx{find(ismember(table_keys,sprintf('fix',p)))}='';
% construct the rest of tags 
trial_trig_idx_dec{find(ismember(table_keys,sprintf('preprobe',p)))}=[trial_auido_end_dec,trial_probe_start_dec-1];
trial_trig_idx_dec{find(ismember(table_keys,sprintf('probe',p)))}=[trial_probe_start_dec,trial_post_probe_start_dec-1];
trial_trig_idx_dec{find(ismember(table_keys,sprintf('extra_probe',p)))}=[trial_post_probe_start_dec,trial_trig_end_dec];
% 
trial_words_idx{find(ismember(table_keys,sprintf('preprobe',p)))}='';
trial_words_idx{find(ismember(table_keys,sprintf('probe',p)))}='';
trial_words_idx{find(ismember(table_keys,sprintf('extra_probe',p)))}='';


pairs_times=cell2mat(trial_trig_idx_dec');
trial_timing_dec{k,1}=table(table_keys',trial_words_idx',pairs_times(:,1),pairs_times(:,2),'VariableNames',{'key','string','start','end'});
trial_timing_dec{k,2}=trial_type;
end 

figure;
histogram(all_timing_diff)
xlabel('ms')
assert(max(abs(all_timing_diff))<50)

end