% function for extraction of subject information for each experiments, 


% Eghbal Hosseini, 24/06/2020
% changelog: 
% HannaH 
% TODO: fix the routine for making tree
function ops_out=create_sub_operation_info_ALBANY(varargin)
p=inputParser();
addParameter(p, 'experiment', 'nlength');
addParameter(p, 'overwrite_all', false);
%addParameter(p, 'n_feat', 28*28);
%addParameter(p, 'beta', 0.01);
%addParameter(p, 'structure', 'partition'); % options : partition , tree
%addParameter(p, 'sigma', 1.5);
%addParameter(p,'norm',true);
%addParameter(p, 'save_path', '~/');
parse(p, varargin{:});
ops = p.Results;
% library of experiment names and subjects 
albany_exp=struct;
exp_dat={...
        {'nlength',arrayfun(@(x) sprintf('AMC%03d',x),[82,83,86,91,92,96,97,99],'uniformoutput',false)};...
        {'langloc',arrayfun(@(x) sprintf('AMC%03d',x),[82,86,88,91,92,96,97,99],'uniformoutput',false)};...
        }; 
end 
% routine for extracting the subject info based on subject id 
function sub_op=extract_sub_info(sub_id)
    % this function does 2 things, 
    % looks at the directory of all subject info load the data, if it was
    % saved and analyzed, then return it, otherwise create it from
    % 
    % generate_all_sub_info output. 
end 
% routine for extracting experiments for a given subject
function exp_info=get_exp_for_sub(sub_id)
    
end 
% creating a dictionarry of all the subjects and their information 
function all_sub_info=genetera_all_sub_info()
    %% AMC026
    all_sub_data=struct;
    sub_dat={{'sub_id'},{'AMC082'};...
            {'GND'},{[47]};...
            {'REF'},{[1]};...
            {'num_chanel'},{[1]};...
            {'bad_channels'},{[1]};...
            {'ecog_channels'},{[1]};...
            {'skull_eeg_channels'},{[1]};...
            {'microphone_channels'},{[1]};...
            {'EMG_channels'},{[1]};...
            {'visual_trigger'},{[1]};...
            {'button_trigger'},{[1]};...
            {'buzzer_trigger'},{[1]};...
            {'audio_trigger'},{[1]};...
            {'audio_trigger'},{[1]};...
            {'channel_label'},[arrayfun(@(x) sprintf('LG_%d',x),[1:64]','uniformoutput',false);...
                              arrayfun(@(x) sprintf('LFM_%d',x),[1:16]','uniformoutput',false);...
                              arrayfun(@(x) sprintf('LFS_%d',x),[1:16]','uniformoutput',false);...
                              arrayfun(@(x) sprintf('LFA_%d',x),[1:4]','uniformoutput',false);...
                              arrayfun(@(x) sprintf('LIP_%d',x),[1:12]','uniformoutput',false);...
                              arrayfun(@(x) sprintf('LIA_%d',x),[1:10]','uniformoutput',false);...
                              arrayfun(@(x) sprintf('RIH_%d',x),[1:4]','uniformoutput',false);...
                              'GND';'REF']
                               };
     all_sub_data.(sub_dat{1,2}{1})=sub_dat;
     % template 
     sub_dat_template={{'sub_id'},{'JaneDoe'};...
            {'GND'},{[]};...
            {'REF'},{[]};...
            {'num_chanel'},{[]};...
            {'bad_channels'},{[]};...
            {'ecog_channels'},{[]};...
            {'skull_eeg_channels'},{[]};...
            {'microphone_channels'},{[]};...
            {'EMG_channels'},{[]};...
            {'visual_trigger'},{[]};...
            {'button_trigger'},{[]};...
            {'buzzer_trigger'},{[]};...
            {'audio_trigger'},{[]};...
            {'audio_trigger'},{[]};...
            {'channel_label'},{[]}
            {'analyzed_by_user'},{[0]}
            {'analyzed_by_user_name'},{['N/A']}};
     all_sub_data.(sub_dat_template{1,2}{1})=sub_dat_template;


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
op_info.analyzed_by_user=0;
op_info.analyzed_by_user_name='Lana';

op_info.seizure_channels=[121,122,105,106,87,95,96];
op_info.exp_name='';
eval(strcat(op_info.subject_name,'_op','.op_info=op_info'));
save(strcat(save_path,op_info.subject_name,'_operation_info.mat'),strcat(op_info.subject_name,'_op'),'-v7.3');
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
