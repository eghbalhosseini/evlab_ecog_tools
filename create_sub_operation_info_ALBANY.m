% create subject info for langloc experiment 
function create_sub_operation_info_ALBANY(save_path)
if ~exist(save_path)
    mkdir(save_path)
end
%% TODO - don't save op info if it already exists!!

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
