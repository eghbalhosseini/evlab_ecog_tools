function [trial_timing_dec] = get_timing_from_json_files_ANNSet1(dataout,all_trial_timing,events_table,audio_align_path,with_wavelet)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
all_timing_diff=[];
initial_table_keys={'trial_start_end',};
all_condition={};
all_audio_length=[];
%
trial_timing_dec={};
if with_wavelet
    audio_path='/Users/eghbalhosseini/MyCodes/MIT_U01_experiments/U01_Expt6_ANNsentSET1/stimuli_w_wavelet';
    delay_val=0.2;
else
    audio_path='/Users/eghbalhosseini/MyCodes/MIT_U01_experiments/U01_Expt6_ANNsentSET1/stimuli_orig';
    delay_val=0;
end 
all_audio_filenames=events_table.final_audio_filename;
for k=1:size(all_trial_timing,1)
    table_keys=initial_table_keys;
    trial_trig_idx_dec=cell(size(initial_table_keys));
    trial_trig_start=all_trial_timing{k,1}(1);
    trial_trig_end=all_trial_timing{k,1}(2);
    [P,Q] = rat(dataout.ops.fsDownsample/dataout.ops.sr);
    trial_trig_start_dec=floor(trial_trig_start*P/Q);
    trial_trig_end_dec=floor(trial_trig_end*P/Q);
    trial_sent_file=all_audio_filenames{k};
    trial_sent_file=erase(trial_sent_file,'.wav');
    sent_id=extract(trial_sent_file,digitsPattern);
    sent_type=extract(trial_sent_file,lettersPattern);
    trial_audio_align=strrep(all_audio_filenames{k},'.wav','_handfix.json');
    fname=strcat(audio_align_path,filesep,trial_audio_align);
    audio_align = jsondecode(fileread(fname));
    word_align=audio_align.words;
    trinscript=strsplit(audio_align.transcript);
    if strcmp(sent_type{1},'sentence')
        trial_type='S';
    elseif strcmp(sent_type{1},'nonword')
        trial_type='N';
    else
        error('trial_type is unkown!')
    end
    all_condition{k}=trial_type;
    assert(length(word_align)==length(trinscript))
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
    table_words={};
    table_words{1}=audio_align.transcript;
    table_words=horzcat(table_words,arrayfun(@(x) word_align(x).word,1:size(word_align,1),'uni',false));
    table_keys=horzcat(table_keys,arrayfun(@(x) ['word_',num2str(x)],1:size(word_align,1),'uni',false));
    audio_fname=strcat(audio_path,filesep,'norm_endfix_filt_',trial_sent_file,'.wav');
    [orig_audio,orig_fs]=audioread(audio_fname);
    %soundsc(orig_audio,orig_fs)
    orig_aud_ts=length(orig_audio)/orig_fs;
    all_audio_length=[all_audio_length;orig_aud_ts];
    trial_aud_ts=diff([all_trial_timing{k,1}])./dataout.MetaTags.SamplingFreq;
    fprintf('timing difference %f ms\n',(orig_aud_ts-trial_aud_ts)*1000)
    all_timing_diff=[all_timing_diff;(orig_aud_ts-trial_aud_ts)*1000];
    [P_aud,Q_aud] = rat(dataout.ops.fsDownsample/orig_fs);
    for p=1:size(word_align,1)
        if ~strcmp(word_align(p).case,'success')
            fprintf('no allignment for word %d in %s \n',p,trial_sent_file) 
        end
        word_st_ts=word_align(p).start+delay_val;
        if word_st_ts==0
            word_st_ts=eps;
        end
        word_end_ts=word_align(p).end+delay_val;
        words_audio{p}=orig_audio(ceil(word_st_ts*orig_fs):floor(word_end_ts*orig_fs));
        try
            words{1,p}=word_align(p).word;
        catch err
            words{1,p}=word_align(p).alignedWord;
        end
        assert(word_end_ts>word_st_ts);
        word_idx_dec=[ceil(word_st_ts*dataout.ops.fsDownsample),floor(word_end_ts*dataout.ops.fsDownsample)];
        trial_trig_idx_dec{find(ismember(table_keys,sprintf('word_%d',p)))}=[trial_trig_start_dec]+word_idx_dec;
    end
    % construct fix index
    trial_trig_idx_dec{find(ismember(table_keys,sprintf('trial_start_end',p)))}=[trial_trig_start_dec,trial_trig_end_dec];
    % construct the rest of tags
    pairs_times=cell2mat(trial_trig_idx_dec');
    trial_timing_dec{k,1}=table(table_keys',table_words',pairs_times(:,1),pairs_times(:,2),'VariableNames',{'key','string','start','end'});
    trial_timing_dec{k,2}=trial_type;
end

figure;
histogram(all_timing_diff)
xlabel('ms')
assert(max(abs(all_timing_diff))<50)

end