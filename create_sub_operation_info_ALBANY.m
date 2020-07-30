% function for extraction of subject information for each experiments, 
% Eghbal Hosseini, 3/07/2020
% changelog: 
% HannaH 
% Eghbal hosseini: updates on saving, and sub_op generations. 
function ops_out=create_sub_operation_info_ALBANY(varargin)
p=inputParser();
addParameter(p, 'experiment', 'MITNLengthSentences');
addParameter(p, 'overwrite', true);
addParameter(p, 'save_path', '~/sub_op_info');
parse(p, varargin{:});
ops = p.Results;
ops_out=ops;
% create op info dir if id doesnt exist 
if exist(ops.save_path,'dir')==0,mkdir(ops.save_path);end
% library of experiment names and subjects 
exp_dat=[...
        ['MITNLengthSentences',{arrayfun(@(x) sprintf('AMC%03d',x),[82,83,86,91,92,96,97,99],'uniformoutput',false)}];...
        ['MITLangloc',{arrayfun(@(x) sprintf('AMC%03d',x),[82,86,88,91,92,96,97,99],'uniformoutput',false)}];...
        ['constituent-bounds',{arrayfun(@(x) sprintf('AMC%03d',x),[83,85,86,91,92,96,97,99,100],'uniformoutput',false)}];...
        ]; 

exp_subs=[exp_dat{contains(exp_dat(:,1),ops.experiment),2}]';
% generate all subject info 
all_sub_info=generate_all_sub_info();
% find subject with info file 
exp_sub_with_info=intersect(exp_subs,fields(all_sub_info));
% create and save sub_op_info
for k=1:length(exp_sub_with_info)
    sub_cell=all_sub_info.(exp_sub_with_info{k});
    op_info=cell2struct(sub_cell(:,2),string([sub_cell{:,1}]));
    eval(strcat(op_info.sub_id,'_op_info','=op_info'));
    eval(sprintf('%s_op_info=op_info',op_info.sub_id));
    eval(sprintf('ops_out.%s_op_info=%s_op_info',op_info.sub_id,op_info.sub_id))
    if ops.overwrite
        save(sprintf('%s/%s_op_info.mat',ops.save_path,op_info.sub_id),sprintf('%s_op_info',op_info.sub_id));
    end 
end 
ops_out.exp_subs=exp_subs;
ops_out.exp_subs_wth_info=exp_sub_with_info;
end 

% routine for extracting experiments for a given subject
function exp_info=get_exp_for_sub(sub_id)

end

% creating a dictionarry of all the subjects and their information 
function all_sub_info=generate_all_sub_info()
    all_sub_info=struct;    
    %% AMC082
    sub_dat={{'sub_id'},                ['AMC082'];...
            {'GND'},                    [127];...
            {'REF'},                    [128];...
            {'num_channels'},           [128];...
            {'bad_channels'},           [];...
            {'ecog_channels'},          [1:126];...
            {'skull_eeg_channels'},     [127:128];...
            {'microphone_channels'},    [192,129];...
            {'EMG_channels'},           [129:134];...
            {'visual_trigger'},         ['DigIO_1'];...
            {'button_trigger'},         ['DigIO_2'];...
            {'buzzer_trigger'},         ['DigIO_3'];...
            {'audio_trigger'},          ['DigIO_4'];...
            {'channel_label'},          [arrayfun(@(x) sprintf('LG_%d',x),[1:64]','uniformoutput',false);...
                                        arrayfun(@(x) sprintf('LFM_%d',x),[1:16]','uniformoutput',false);...
                                        arrayfun(@(x) sprintf('LFS_%d',x),[1:16]','uniformoutput',false);...
                                        arrayfun(@(x) sprintf('LFA_%d',x),[1:4]','uniformoutput',false);...
                                        arrayfun(@(x) sprintf('LIP_%d',x),[1:12]','uniformoutput',false);...
                                        arrayfun(@(x) sprintf('LIA_%d',x),[1:10]','uniformoutput',false);...
                                        arrayfun(@(x) sprintf('RIH_%d',x),[1:4]','uniformoutput',false);...
                                        'GND';'REF'];...
            {'visually_inspected'},       [0];...
            {'visually_inspected_by'},  ['N/A']};      
       
     % add subject to all fub definition
     all_sub_info.(sub_dat{1,2})=sub_dat;
     
     %% AMC083
     sub_dat_template={{'sub_id'},      ['AMC083'];...
            {'GND'},                    [];...
            {'REF'},                    [];...
            {'num_channels'},           [];...
            {'bad_channels'},           [];...
            {'ecog_channels'},          [];...
            {'skull_eeg_channels'},     [];...
            {'microphone_channels'},    [];...
            {'EMG_channels'},           [];...
            {'visual_trigger'},         [];...
            {'button_trigger'},         [];...
            {'buzzer_trigger'},         [];...
            {'audio_trigger'},          [];...
            {'channel_label'},          [];...
            {'visually_inspected'},       [0];...
            {'visually_inspected_by'},  ['N/A']};
        % add subject to all fub definition
     all_sub_info.(sub_dat_template{1,2})=sub_dat_template;
    %% AMC086
    sub_dat_template={{'sub_id'},      ['AMC086'];...
            {'GND'},                    [];...
            {'REF'},                    [];...
            {'num_channels'},             [];...
            {'bad_channels'},           [];...
            {'ecog_channels'},          [];...
            {'skull_eeg_channels'},     [];...
            {'microphone_channels'},    [];...
            {'EMG_channels'},           [];...
            {'visual_trigger'},         [];...
            {'button_trigger'},         [];...
            {'buzzer_trigger'},         [];...
            {'audio_trigger'},          [];...
            {'channel_label'},          [];...
            {'visually_inspected'},       [0];...
            {'visually_inspected_by'},  ['N/A']};
        % add subject to all fub definition
     all_sub_info.(sub_dat_template{1,2})=sub_dat_template;
     %% AMC088
    sub_dat_template={{'sub_id'},      ['AMC088'];...
            {'GND'},                    [];...
            {'REF'},                    [];...
            {'num_channels'},             [];...
            {'bad_channels'},           [];...
            {'ecog_channels'},          [];...
            {'skull_eeg_channels'},     [];...
            {'microphone_channels'},    [];...
            {'EMG_channels'},           [];...
            {'visual_trigger'},         [];...
            {'button_trigger'},         [];...
            {'buzzer_trigger'},         [];...
            {'audio_trigger'},          [];...
            {'channel_label'},          [];...
            {'visually_inspected'},       [0];...
            {'visually_inspected_by'},  ['N/A']};
        % add subject to all fub definition
     all_sub_info.(sub_dat_template{1,2})=sub_dat_template;
        %% AMC091
     sub_dat_template={{'sub_id'},      ['AMC091'];...
            {'GND'},                    [];...
            {'REF'},                    [];...
            {'num_channels'},             [];...
            {'bad_channels'},           [];...
            {'ecog_channels'},          [];...
            {'skull_eeg_channels'},     [];...
            {'microphone_channels'},    [];...
            {'EMG_channels'},           [];...
            {'visual_trigger'},         [];...
            {'button_trigger'},         [];...
            {'buzzer_trigger'},         [];...
            {'audio_trigger'},          [];...
            {'channel_label'},          [];...
            {'visually_inspected'},       [0];...
            {'visually_inspected_by'},  ['N/A']};
        % add subject to all fub definition
     all_sub_info.(sub_dat_template{1,2})=sub_dat_template;
     %% AMC092
     sub_dat_template={{'sub_id'},      ['AMC092'];...
            {'GND'},                    [];...
            {'REF'},                    [];...
            {'num_channels'},             [];...
            {'bad_channels'},           [];...
            {'ecog_channels'},          [];...
            {'skull_eeg_channels'},     [];...
            {'microphone_channels'},    [];...
            {'EMG_channels'},           [];...
            {'visual_trigger'},         [];...
            {'button_trigger'},         [];...
            {'buzzer_trigger'},         [];...
            {'audio_trigger'},          [];...
            {'channel_label'},          [];...
            {'visually_inspected'},       [0];...
            {'visually_inspected_by'},  ['N/A']};
        % add subject to all fub definition
     all_sub_info.(sub_dat_template{1,2})=sub_dat_template;
     %% AMC096
     sub_dat_template={{'sub_id'},      ['AMC096'];...
            {'GND'},                    [];...
            {'REF'},                    [];...
            {'num_channels'},             [];...
            {'bad_channels'},           [];...
            {'ecog_channels'},          [];...
            {'skull_eeg_channels'},     [];...
            {'microphone_channels'},    [];...
            {'EMG_channels'},           [];...
            {'visual_trigger'},         [];...
            {'button_trigger'},         [];...
            {'buzzer_trigger'},         [];...
            {'audio_trigger'},          [];...
            {'channel_label'},          [];...
            {'visually_inspected'},       [0];...
            {'visually_inspected_by'},  ['N/A']};
        % add subject to all fub definition
     all_sub_info.(sub_dat_template{1,2})=sub_dat_template;
          %% AMC097
     sub_dat_template={{'sub_id'},      ['AMC097'];...
            {'GND'},                    [];...
            {'REF'},                    [];...
            {'num_channels'},             [];...
            {'bad_channels'},           [];...
            {'ecog_channels'},          [];...
            {'skull_eeg_channels'},     [];...
            {'microphone_channels'},    [];...
            {'EMG_channels'},           [];...
            {'visual_trigger'},         [];...
            {'button_trigger'},         [];...
            {'buzzer_trigger'},         [];...
            {'audio_trigger'},          [];...
            {'channel_label'},          [];...
            {'visually_inspected'},       [0];...
            {'visually_inspected_by'},  ['N/A']};
        % add subject to all fub definition
     all_sub_info.(sub_dat_template{1,2})=sub_dat_template;
     
     %% AMC099
     sub_dat_template={{'sub_id'},      ['AMC099'];...
            {'GND'},                    [];...
            {'REF'},                    [];...
            {'num_channels'},             [];...
            {'bad_channels'},           [];...
            {'ecog_channels'},          [];...
            {'skull_eeg_channels'},     [];...
            {'microphone_channels'},    [];...
            {'EMG_channels'},           [];...
            {'visual_trigger'},         [];...
            {'button_trigger'},         [];...
            {'buzzer_trigger'},         [];...
            {'audio_trigger'},          [];...
            {'channel_label'},          [];...
            {'visually_inspected'},       [0];...
            {'visually_inspected_by'},  ['N/A']};
        % add subject to all fub definition
     all_sub_info.(sub_dat_template{1,2})=sub_dat_template;
     
     

end 




%% TODO - don't save op info if it already exists!!
function trash_can()
%% AMC 26
op_info=struct;
op_info.subject_name='AMC026';

op_info.GND=[47];
op_info.Ref=1;
op_info.num_channels=128;
op_info.bad_channels = [108 121:128];

eval(strcat(op_info.subject_name,'_op','.op_info=op_info'));
save(strcat(save_path,op_info.subject_name,'_operation_info.mat'),strcat(op_info.subject_name,'_op'),'-v7.3');
%% 
%% AMC 29
op_info=struct;
op_info.subject_name='AMC029';

op_info.GND=[7];
op_info.Ref=10;
op_info.num_channels=128;
op_info.bad_channels = [70:75 78:84 87:93 96:102 105:111 114:120 124];

eval(strcat(op_info.subject_name,'_op','.op_info=op_info'));
save(strcat(save_path,op_info.subject_name,'_operation_info.mat'),strcat(op_info.subject_name,'_op'),'-v7.3');

%% 
%% AMC 30
op_info=struct;
op_info.subject_name='AMC030';

op_info.GND=[];
op_info.Ref=[];
op_info.num_channels=[];
op_info.bad_channels = [];

eval(strcat(op_info.subject_name,'_op','.op_info=op_info'));
save(strcat(save_path,op_info.subject_name,'_operation_info.mat'),strcat(op_info.subject_name,'_op'),'-v7.3');
%% 
op_info=struct;
op_info.subject_name='AMC031';

op_info.GND=[78];
op_info.Ref=[97,99];
op_info.num_channels=112;
op_info.bad_channels = [29:31 37:38 61 72:73 77:80 97:112];

eval(strcat(op_info.subject_name,'_op','.op_info=op_info'));
save(strcat(save_path,op_info.subject_name,'_operation_info.mat'),strcat(op_info.subject_name,'_op'),'-v7.3');
%% 
%% 
op_info=struct;
op_info.subject_name='AMC037';

op_info.GND=[41];
op_info.Ref=[42];
op_info.num_channels=134;
op_info.bad_channels = [3 5 7 67 80 99 124 133];

eval(strcat(op_info.subject_name,'_op','.op_info=op_info'));
save(strcat(save_path,op_info.subject_name,'_operation_info.mat'),strcat(op_info.subject_name,'_op'),'-v7.3');
%% 
op_info=struct;
op_info.subject_name='AMC038';

op_info.GND=[16];
op_info.Ref=[1];
op_info.num_channels=98;
op_info.bad_channels = [21 26  53 56 60 84 92];

eval(strcat(op_info.subject_name,'_op','.op_info=op_info'));
save(strcat(save_path,op_info.subject_name,'_operation_info.mat'),strcat(op_info.subject_name,'_op'),'-v7.3');
%% 
op_info=struct;
op_info.subject_name='AMC044';

op_info.GND=[];
op_info.Ref=[];
op_info.num_channels=68;
op_info.bad_channels = [5 7 10 37:68]; 

eval(strcat(op_info.subject_name,'_op','.op_info=op_info'));
save(strcat(save_path,op_info.subject_name,'_operation_info.mat'),strcat(op_info.subject_name,'_op'),'-v7.3');

%%
%% AMC 82
op_info=struct;
op_info.subject_name='AMC082';

op_info.GND=[127];
op_info.Ref=[128];
op_info.num_channels=128;
op_info.ecog_channels=1:126;
op_info.skull_eeg_channels=127:128;
op_info.microphone_channels=[192,129];
op_info.EMG_channels=[129:134];
op_info.visual_trigger='DigIO_1';
op_info.button_trigger='DigIO_2';
op_info.buzzer_trigger='DigIO_3';
op_info.audio_trigger='DigIO_4';
op_info.bad_channels = [];
op_info.task_file_location='/DAY6/MITNaturalisticStoriesTask/ECOG001/';
labels= [arrayfun(@(x) sprintf('LG_%d',x),[1:64]','uniformoutput',false);...
        arrayfun(@(x) sprintf('LFM_%d',x),[1:16]','uniformoutput',false);...
        arrayfun(@(x) sprintf('LFS_%d',x),[1:16]','uniformoutput',false);...
        arrayfun(@(x) sprintf('LFA_%d',x),[1:4]','uniformoutput',false);...
        arrayfun(@(x) sprintf('LIP_%d',x),[1:12]','uniformoutput',false);...
        arrayfun(@(x) sprintf('LIA_%d',x),[1:10]','uniformoutput',false);...
        arrayfun(@(x) sprintf('RIH_%d',x),[1:4]','uniformoutput',false);...
        'GND';'GND'];

op_info.channel_labels=labels;
op_info.visually_inspected = 0;

op_info.seizure_channels=[121,122,105,106,87,95,96];
eval(strcat(op_info.subject_name,'_op','.op_info=op_info'));
save(strcat('~/Desktop/ECOG/subject_op_info_MASTER/',op_info.subject_name,'_operation_info.mat'),strcat(op_info.subject_name,'_op'),'-v7.3');
% op_info=struct;
% op_info.subject_name='AMC082';
% 
% op_info.GND=[];
% op_info.Ref=[];
% op_info.num_channels=128;
% op_info.bad_channels = []; 
% 
% eval(strcat(op_info.subject_name,'_op','.op_info=op_info'));
% save(strcat(save_path,op_info.subject_name,'_operation_info.mat'),strcat(op_info.subject_name,'_op'),'-v7.3');
end
