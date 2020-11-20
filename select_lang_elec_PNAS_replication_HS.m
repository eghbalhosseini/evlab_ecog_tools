% selection of language specific electrodes by comparing electrode response
% via the method in PNAS 2016 paper

%for each electrode the script, 
%1. gets the response amplitudes to condition1 (or pretrial)
%2. get the response amplitudes to condition2 (or pretrial)


% - first computed the mean of the signal envelope for each of the eight word positions (time-locked to the onset of each word/ nonword and averaging over the presentation window) in each trial for each condition in each electrode, using data from the odd-numbered runs
% only (see SI Appendix, Part F for sample EOIs without within-position averaging).
% - We then computed the mean across the eight word positions in each
% trial for each condition in each electrode.
% - Finally, we correlated the trial means with a vector of condition labels (sentences = 1, nonword-lists = −1). The resulting Spearman’s ρ provided a benchmark against which to test the significance of any positive correlations.
% - The condition labels vector was randomly reordered (via a permutation test without replacement) and a new Spearman’s ρ was computed, and this process was repeated 1,000 times.
% - The fraction of correlations from randomly assigned labels that produced a
% higher ρ than the benchmark correlation became our P value. Electrodes with P ≤ 0.01 and a positive ρ were included in step 2.
%if significant, adds the electrode to language responsive electrodes

%% clean up workspace
clear all
close all

%% specify parameters
subject_id = 'AMC096'; %must be string
experiment_name = 'MITLangloc'; %must be string
condition1 = [{'Sentences'}]; %must be a list
condition2 = [{'Jabberwocky'}]; %must be a list
conditions = {condition1,condition2};
use_pretrials = [{false},{false}]; %use same indices as conditions, must be logicals
target_word_idxs = [{1:12},{1:12}]; %use same indices as conditions, word(s) to get the average response amplitude to, must be a list

data_to_use = "hilbert_zs"; %must be string
verbose = true;

num_permutations = 1000;
p_threshold = 0.01;

global testing;
testing = false;

% subject_id = 'AMC096';
% experiment_name = 'MITNLengthSentences';
% condition1 = [{'3sents_8words_intact'},{'6sents_4words_intact'}, {'1sent_24words_intact'}];
% condition2 = [{'3sents_8words_intact'},{'6sents_4words_intact'}, {'1sent_24words_intact'}];
% %condition1 = [{'3sents_8words_scrambled'},{'6sents_4words_scrambled'}, {'1sent_24words_scrambled'}];
% %condition2 = [{'3sents_8words_scrambled'},{'6sents_4words_scrambled'}, {'1sent_24words_scrambled'}];
% cond1_use_pretrial = true;
% cond2_use_pretrial = false;
% cond1_target_word_idx = [0];
% cond2_target_word_idx = [8];
% data_to_use = "hilbert_zs"; %options: hilbert_zs;
% verbose = true;



%% specify where the data is
addpath(genpath('~/GitHub/evlab_ecog_tools'));
addpath(genpath('~/GitHub/evlab_ecog_tools/ecog-filters/'));
addpath(genpath('~/GitHub/evlab_ecog_tools/albany_mex_files'));
addpath(genpath('~/GitHub/evlab_matlab_tools/Colormaps'));
code_path='~/GitHub/evlab_ecog_tools';
ecog_path = '~/Desktop/ECOG';

crunched_data_path = [ecog_path filesep 'crunched' filesep experiment_name filesep]; %save it into an experiment specific folder

%path specifically for expt sub_op_info
expt_sub_op_info_path = [crunched_data_path 'sub_op_info_' experiment_name filesep];

%% extract subject's data for the given experiment
d_data= dir(strcat(crunched_data_path,filesep,subject_id,'*_crunched.mat'));
fprintf(' %d .mat files were found \n', length(d_data));
d_data=arrayfun(@(x) strcat(d_data(x).folder,filesep,d_data(x).name),[1:length(d_data)]','uni',false);

check_parameters(conditions, use_pretrials, target_word_idxs, verbose);
if(verbose)fprintf("Parameters checked. \n");end;
if(verbose)fprintf("Using subject %s's %s data from %s for language electrode selection \n", subject_id,data_to_use, experiment_name);end


[mean_amplitudes] = extract_subj_data(d_data, ...
                            data_to_use,... 
                            conditions,use_pretrials,target_word_idxs,...
                            verbose);

%% get language electrodes

[lang_electrodes_list_permutations, lang_electrodes_list_t_test, electrodes] = determine_lang_electrodes(mean_amplitudes{1,4},mean_amplitudes{2,4},p_threshold, num_permutations, verbose);


%%
function [] = check_parameters(conditions,use_pretrials, target_word_idxs,verbose)
if ~(iscell(conditions))
    error('conditions must be a cell array of lists of strings, with each entry contains a list with the conditions to include for each condition (can be just one condition)');
end

if ~iscell(use_pretrials)
    error('use_pretrials must be a logical cell array, specifying whether to use the pretrial for the corresponding condition in conditions');
end

if ~iscell(target_word_idxs)
    error('target_word_idxs must be a cell array, with each entry containing a list of word indices to extract average responses to for the corresponding condition in conditions');
end

if ~islogical(verbose)
    error('verbose must be a logical. true to print out updates, false to suppress outputs.');
end

%check pretrial
for i=1:length(use_pretrials)
    if(use_pretrials{i})
        if(length(target_word_idxs{i})>1)
            error('if you are using the pretrials of a condition, you may only specify one word index (does not matter which one, it will not be used, only the pretrial period before the entire trial will be used). if you specify more than one word index, it will average the pretrials across the trials and cause duplicate values')
        end
    end
end

end
function [list, list_t_test,lang_electrodes] = determine_lang_electrodes(cond1,cond2, p_threshold, num_permutations, verbose)
%each row of cond1 and cond2 contain one electrode's avg response
%amplitudes over all the relevant trials

%t-test comparing the rows of cond1 and cond2
num_electrodes = size(cond1,1);
list_t_test = [];
t_test_results = [1:num_electrodes];
if(verbose)fprintf("Conducting t-tests comparing cond1 and cond2 responses for %d electrodes\n",num_electrodes);end
for i=1:num_electrodes
    cond1_vector = cond1(i,:);
    cond2_vector = cond2(i,:);
    result = ttest2(cond1_vector,cond2_vector, 'Alpha', p_threshold);
    t_test_results(i) = result;
    if (result==1)
        list_t_test = [list_t_test i]; %save the electrode number in the list
    end
end

%lang_electrodes = t_test_results;


%correlated the trial means with a vector of condition labels (sentences = 1, nonword-lists = −1). 
%The resulting Spearman’s ρ provided a benchmark against which to test the significance of any positive correlations. 

trial_amplitudes = cat(2, cond1, cond2);
condition_labels = cat(2, cond1*0+1, cond2*0-1);
trial_amplitudes_transposed = trial_amplitudes';
condition_labels_transposed = condition_labels';

[rho, p_value] = corr(trial_amplitudes_transposed,condition_labels_transposed,'Type','Spearman');
rho_base = rho(:,1);
%The condition labels vector was randomly reordered (via a permutation test without replacement) 
%and a new Spearman’s ρ was computed, and this process was repeated 1,000 times.
if(verbose)fprintf("Performing %d permutations\n", num_permutations);end
rho_permutations = [];
for i=1:num_permutations
    random_index=randperm(size(condition_labels,2));
    random_labels = condition_labels(:,random_index);
    random_labels_transposed = random_labels';
    
    [rho, p_value] = corr(trial_amplitudes_transposed,random_labels_transposed,'Type','Spearman');
    %save the rho for comparison to unscrambled labels
    rho_permutations = cat(2,rho_permutations,rho(:,1));
end

%The fraction of correlations from randomly assigned labels that produced a
%higher ρ than the benchmark correlation became our P value. 
%Electrodes with P ≤ 0.01 and a positive ρ were selected

p_fraction = sum(rho_permutations>repmat(rho_base,[1,num_permutations]),2)./size(rho_permutations,2);
p_significant = p_fraction<p_threshold;
rho_positive = rho_base > 0;

lang_electrodes = [p_significant&rho_positive]';

list = 1:length(lang_electrodes);
list = list(lang_electrodes);

end
function [cond_mean_trial_amplitudes] = extract_subj_data(d_data,data_to_use, conditions, use_pretrials, target_word_idxs,verbose)
global testing;
num_cond = length(conditions);
cond_amplitudes = cell(num_cond,1);
cond_mean_trial_amplitudes = cell(num_cond);

num_sessions = length(d_data);
if(testing)
    fprintf("WARNING: IN TESTING MODE, ONLY USING SESSION 1 DATA\n")
    num_sessions = 1;
end
%compile condition trial data over all of the sessions
for k=1:num_sessions
    if(verbose)fprintf('Session %d \n', k);end
    subj=load(d_data{k});
    for i=1:num_cond
        amplitudes = extract_condition_amplitudes(subj,data_to_use, conditions{i}, use_pretrials{i},target_word_idxs{i}, verbose);
        cond_amplitudes{i} = cat(1,cond_amplitudes{i},amplitudes);
    end
    
end

for i=1:num_cond
    current_cond_amplitudes = cond_amplitudes{i};

    %put data into easy to use format
    cond_word_amplitudes = cell2mat(reshape(current_cond_amplitudes,1,[]));

    num_words = length(target_word_idxs{i});
    num_trials = size(cond_word_amplitudes,2)/num_words;

    %compute the mean across the word positions in each trial for each electrode. 
    if(verbose)fprintf('Computing average across the word positions for each trial in each electrode in condition %d...\n',i);end;
    num_electrodes = size(cond_word_amplitudes,1);
    cond_trial_amplitudes = zeros(num_electrodes,num_trials);
    for k=1:num_electrodes
        for j=1:num_trials
            start_idx=1;
            if(j>1)
                start_idx = start_idx + (j-1)*num_words;
            end
            end_idx = start_idx+num_words-1;
            cond_trial_amplitudes(k,j) = mean(cond_word_amplitudes(k, start_idx:end_idx));
        end
    end
    
    pretrial_message = "not_pretrial";
    if(use_pretrials{i})
        pretrial_message = "pretrial";
    end
    
    cond_mean_trial_amplitudes{i,1} = conditions{i};
    cond_mean_trial_amplitudes{i,2} = pretrial_message;
    cond_mean_trial_amplitudes{i,3} = target_word_idxs{i};
    cond_mean_trial_amplitudes{i,4} = cond_trial_amplitudes;
end

end
function [avg_amplitudes] = extract_condition_amplitudes(subj,data_to_use, conditions,use_pretrial, target_word_idx,verbose)
avg_amplitudes = [];
avg_word_amplitudes = [];
if(verbose)fprintf('extracting average response amplitudes to... \n');end
for i=1:length(conditions)
    current_condition = string(conditions(i));
    for j=1:length(target_word_idx)
        extracted_amplitudes = get_avg_amplitudes(subj,data_to_use, current_condition, use_pretrial, target_word_idx(j), verbose);
        avg_word_amplitudes = cat(1,avg_word_amplitudes, extracted_amplitudes);
    end
    avg_amplitudes = cat(1,avg_amplitudes, avg_word_amplitudes);
end
end
% subj: struct, a subj's crunched materials, containing info and data structs
% data_to_use: string, which type of data from the data struct to use
% cond: string, the condition to extract avg amplitudes from
% use_pretrial: boolean, specifies whether you should extract the amplitude from the pretrials in the relevant trials
% target_word_idx: int, only used if use_pretrial is false, specifies the word to collect avg amplitude from
function [avg_amplitudes] = get_avg_amplitudes(subj, data_to_use, cond, use_pretrial,target_word_idx, verbose)
clear info;
clear data;
subj_sess_id=fieldnames(subj);
subj=subj.(subj_sess_id{1}); %get rid of top layer of struct
info = subj.info;
data = subj.data;

%get strings of the condition names on each trial
trial_type_str = cellfun(@string, info.condition_name,'UniformOutput',false); 
%convert the logical arrays (one per word) into 1 dimension (one per trial)
trial_type_unique = cellfun(@unique,trial_type_str,'UniformOutput',false);
%get the indices of the relevant trials
relevant_trials = cellfun(@(x) x==cond, trial_type_unique);

relevant_trial_data=data(relevant_trials);
if(use_pretrial)
    %get the average pretrial response amplitude for all relevant trials
    data_name = strcat("signal_ave_pre_trial_",data_to_use, "_downsample");
    if(verbose)fprintf('pretrial before %s\n',cond);end
    ave_electrodes=cellfun(@(x) x.signal_ave_pre_trial_hilbert_zs_downsample, relevant_trial_data,'UniformOutput',false);
else
    %get the average response amplitude to the target word for all relevant
    %trials
    data_name = strcat("signal_ave_", data_to_use, "_downsample_parsed");
    if(verbose)fprintf('word %d in %s\n', target_word_idx, cond);end
    ave_electrodes=cellfun(@(x) x.signal_ave_hilbert_zs_downsample_parsed{target_word_idx,1}, relevant_trial_data,'UniformOutput',false);
end

avg_amplitudes = ave_electrodes;
end




      