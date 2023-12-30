function [dataout]=filter_channels_MGH_NSX(datafile,subject_op_info)
% varargin : 
%1 subject operation info


ecog.param.filter_car  = 1;     % 0 == off % 1 == on
ecog.param.filter_type = 'IIR'; % IIR | FIR                             

% --- gerv's bands ---
ecog.param.f_low_stop  = [ 4,15,30, 64,220];
ecog.param.f_low_pass  = [ 8,18,35, 70,250];
ecog.param.f_high_pass = [12,25,50,170,350];
ecog.param.f_high_stop = [16,30,55,200,380];
ecog.param.topo_size   = [2,3];


% --- highpass filter --- 
ecog.param.highpass.Wp = 0.50; % Hz
ecog.param.highpass.Ws = 0.05; % Hz
ecog.param.highpass.Rp = 3;    % dB
ecog.param.highpass.Rs = 30;   % dB

% --- lowpass filter --- 
ecog.param.lowpass.Wp = 100;     % Hz
ecog.param.lowpass.Ws = 110;     % Hz
ecog.param.lowpass.Rp = 3;     % dB
ecog.param.lowpass.Rs = 30;    % dB

%---- downsampling rate ----
ecog.param.samplingrate = 500; % decimation sampling rate 

audio.param.filter_type = 'IIR'; % IIR | FIR
% --- speech ---
audio.param.f_low_stop  = [  40];
audio.param.f_low_pass  = [  80];
audio.param.f_high_pass = [6000];
audio.param.f_high_stop = [8000];
audio.param.samplingrate = 500;
op_ecog_labels=subject_op_info.op_info.ecog_labels;
op_dig_loc=cell2mat(regexp(op_ecog_labels,'w*_'))-1;
op_ch_tags=arrayfun(@(x) op_ecog_labels{x}(1:op_dig_loc(x)),1:length(op_dig_loc),'uni',false);


unique_op_tags=unique(op_ch_tags,'stable');
ch_names=datafile.ElectrodesInfo.Label;
param=ecog.param;

param.channels_ID=datafile.ElectrodesInfo.Label;
dig_loc=regexp(ch_names,'w*_');
% 
valid_label=find(cell2mat(cellfun(@(x) ~isempty(x),dig_loc,'uni',false)));
ch_tag=arrayfun(@(x) ch_names{x}(1:dig_loc{x}-1),valid_label,'uni',false);
ch_num=arrayfun(@(x) ch_names{x}((dig_loc{x}+1):end),valid_label,'uni',false);
unique_tags=unique(ch_tag,'stable');
ecog_tags=intersect(unique_tags,unique_op_tags,'stable')';
% 
ecog_loc_for_tag=cellfun(@(x) find(...
    cell2mat(...
    cellfun(@(y) ~isempty(regexp(y,[x,'_*\d'])),ch_names,'uni',false))==1)',ecog_tags,'uni',false);
param.ecog_channels=cell2mat(ecog_loc_for_tag);
param.ecog_channels_labels=ch_names(param.ecog_channels)';

% bipolar versions 


if isempty(param.ecog_channels)
    disp('no viable electrode in the dataset, aborting')
    dataout=[]
else 


%% 
parameters.SamplingRate=datafile.MetaTags.SamplingFreq;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE HIGH-PASS FILTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define passband, stopband and attenuation
highpass{1}.Wp = param.highpass.Wp/(parameters.SamplingRate/2); 
highpass{1}.Ws = param.highpass.Ws/(parameters.SamplingRate/2);
highpass{1}.Rp = param.highpass.Rp; 
highpass{1}.Rs = param.highpass.Rs;
% calculate the minimum filter order for butterworth filter
[highpass{1}.n,highpass{1}.Wn] = buttord(highpass{1}.Wp,highpass{1}.Ws,highpass{1}.Rp,highpass{1}.Rs);
highpass{1}.n = highpass{1}.n + rem(highpass{1}.n,2);
% caclulate the filter coefficients in Zero-Pole-Gain design
[highpass{1}.z,highpass{1}.p,highpass{1}.k] = butter(highpass{1}.n,highpass{1}.Wn,'high');
[highpass{1}.sos,highpass{1}.g]=zp2sos(highpass{1}.z,highpass{1}.p,highpass{1}.k);
highpass{1}.h=dfilt.df2sos(highpass{1}.sos,highpass{1}.g);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE BAND-PASS FILTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(param.filter_type,'FIR'), 
    % FIR bandpass definition using Kaiser filter
    for idx=1:length(param.f_low_stop),
        % define passband, stopband and attenuation
        bandpass{idx}.fcuts = [param.f_low_stop(idx) param.f_low_pass(idx) param.f_high_pass(idx) param.f_high_stop(idx)];
        bandpass{idx}.mags  = [0 1 0];
        bandpass{idx}.devs  = [0.01 0.05 0.01];        
        % calculate the minimum filter order        
        [bandpass{idx}.n,bandpass{idx}.Wn,bandpass{idx}.beta,bandpass{idx}.ftype] = kaiserord(bandpass{idx}.fcuts,bandpass{idx}.mags,bandpass{idx}.devs,parameters.SamplingRate);
        bandpass{idx}.n = bandpass{idx}.n + rem(bandpass{idx}.n,2);        
        % caclulate the filter coefficients in a,b design        
        bandpass{idx}.b = fir1(bandpass{idx}.n,bandpass{idx}.Wn,bandpass{idx}.ftype,kaiser(bandpass{idx}.n+1,bandpass{idx}.beta),'noscale');
        [bandpass{idx}.H,bandpass{idx}.f] = freqz(bandpass{idx}.b,1,1024,parameters.SamplingRate);

    end   
elseif strcmp(param.filter_type,'IIR'),
    % IIR bandpass definition using Butterworth filter
    for idx=1:length(param.f_low_stop),
        % define passband, stopband and ripples
        bandpass{idx}.Wp = [param.f_low_pass(idx) param.f_high_pass(idx)]/(parameters.SamplingRate/2); 
        bandpass{idx}.Ws = [param.f_low_stop(idx) param.f_high_stop(idx)]/(parameters.SamplingRate/2);
        bandpass{idx}.Rp = 3; 
        bandpass{idx}.Rs = 30;
        % calculate the minimum filter order
        [bandpass{idx}.n,bandpass{idx}.Wn] = buttord(bandpass{idx}.Wp,bandpass{idx}.Ws,bandpass{idx}.Rp,bandpass{idx}.Rs);
        bandpass{idx}.n = bandpass{idx}.n + rem(bandpass{idx}.n,2);
        % caclulate the filter coefficients in Zero-Pole-Gain design
        [bandpass{idx}.z, bandpass{idx}.p, bandpass{idx}.k] = butter(bandpass{idx}.n,bandpass{idx}.Wn,'bandpass');
        [bandpass{idx}.sos,bandpass{idx}.g]=zp2sos(bandpass{idx}.z,bandpass{idx}.p,bandpass{idx}.k);
        bandpass{idx}.h=dfilt.df2sos(bandpass{idx}.sos,bandpass{idx}.g);
    end     
else 
    error('error: invalid filter type in param.filter_type'); 
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE LOW-PASS FILTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define passband, stopband and ripples 
lowpass{1}.Wp = param.lowpass.Wp/(parameters.SamplingRate/2); 
lowpass{1}.Ws = param.lowpass.Ws/(parameters.SamplingRate/2);
lowpass{1}.Rp = param.lowpass.Rp; 
lowpass{1}.Rs = param.lowpass.Rs;

% calculate the minimum filter order
[lowpass{1}.n,lowpass{1}.Wn] = buttord(lowpass{1}.Wp,lowpass{1}.Ws,lowpass{1}.Rp,lowpass{1}.Rs);
lowpass{1}.n = lowpass{1}.n + rem(lowpass{1}.n,2);

% caclulate the filter coefficients in Zero-Pole-Gain design
[lowpass{1}.z,lowpass{1}.p,lowpass{1}.k] = butter(lowpass{1}.n,lowpass{1}.Wn,'low');
[lowpass{1}.sos,lowpass{1}.g]=zp2sos(lowpass{1}.z,lowpass{1}.p,lowpass{1}.k);
lowpass{1}.h=dfilt.df2sos(lowpass{1}.sos,lowpass{1}.g);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE LINE-NOISE FILTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
peak.fcenter = 60;
peak.bw      = 0.001;

% calculate the IIR-peak filter coefficients in a,b format 
peak.wo = peak.fcenter/(parameters.SamplingRate/2);  
peak.bw = peak.bw;
[peak.b,peak.a] = iirpeak(peak.wo,peak.bw);  

% define the harmonics of line noise frequency
param.filter.notch.fcenter = [60,120,180,240];
param.filter.notch.bw      = ones(1,length(param.filter.notch.fcenter)).*0.001;

% calculate the IIR-peak filter coefficients in a,b format 
for idx = 1:length(param.filter.notch.fcenter),
    notch{idx}.wo = param.filter.notch.fcenter(idx)/(parameters.SamplingRate/2);  
    notch{idx}.bw = param.filter.notch.bw(idx);
    [notch{idx}.b,notch{idx}.a] = iirnotch(notch{idx}.wo,notch{idx}.bw);  
end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD AND CONCATENATE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signal = double([]);
num_of_chan=length(param.ecog_channels);
% go through all data files
%% highpass filter
signal_loop=(datafile.Data)';
fprintf(1,'[');
% this reorders the data into param.ecog_channels format 
for idx=1:num_of_chan
    warning('off', 'signal:filtfilt:ParseSOS');
    idx_channel=param.ecog_channels(idx);
    fprintf('%d --> %d \n',idx_channel,idx)
    assert(idx_channel==idx)
    signal(:,idx) = (filtfilt(highpass{1}.sos,highpass{1}.g,double(signal_loop(:,idx_channel)))); %#ok<PFBNS>
    fprintf(1,'.');
end
fprintf(1,'] done\n');

buffer_before = parameters.SamplingRate * 1; % 1 sec
buffer_after  = parameters.SamplingRate * 2; % 2 sec
clear signal_loop
%% plot the dataset throughout 

% get sections 
t_length=25e4;
kk_total=ceil(size(signal,1)/t_length);
x_cell=mat2cell(signal',ones(1,size(signal,2)));
min_max=mean(cell2mat(cellfun(@(y) [min(y),max(y)],x_cell,'UniformOutput',false)));
x_norm_cell=cellfun(@(x) (x-min_max(1))./(min_max(2)-min_max(1)),x_cell,'UniformOutput',false);
% 
figure(1);
%colors=inferno(size(x_norm_cell,1));
col_inf=inferno(floor(.8*size(x_norm_cell,1)));
col_vir=viridis(floor(.8*size(x_norm_cell,1)));
colors=[col_vir(1:floor(size(x_norm_cell,1)/2),:);col_inf(1:(floor(size(x_norm_cell,1)/2)+1),:)];
for kk=1:kk_total
    clf;
    set(gcf,'color',[.7,.7,.7],'position',[17,1,2454,1337]);
    ax=axes('position',[.02,.02,.95,.95]);
    hold on
    session_end=min([(kk*t_length),size(signal,1)]);
    t_window=((kk-1)*t_length+1):2:(session_end);
    arrayfun(@(x) plot(t_window,x_norm_cell{x}(t_window)+x-.5,'color',colors(x,:)),[1:size(x_norm_cell,1)])
    set(ax,'ytick',[1:size(x_norm_cell,1)]);
    set(ax,'yticklabel','')
    arrayfun(@(x) text(t_window(1),x,num2str(param.ecog_channels(x)),...
        'Color',colors(x,:),'HorizontalAlignment','right','VerticalAlignment','middle','fontsize',6),[1:size(x_norm_cell,1)]);
     arrayfun(@(x) text(t_window(end),x,param.ecog_channels_labels{x},...
        'Color',colors(x,:),'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',6),[1:size(x_norm_cell,1)]);
    
    set(ax,'ylim',[0,size(x_norm_cell,1)+1])
    ax.TickLength=[0.001 0.0015];
    shg;
    title({'pre common source averaging',[num2str(kk),' of ', num2str(kk_total)]})
    pause
end 
close all; 
  
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEASSURE LINE-NOISE POWER BEFORE SIGNAL PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, '> Meassuring 60 Hz noise power before signal processing \n');
fprintf(1,'[');
parfor idx_channel=1:size(signal,2),
    % calculate average root-mean-square of the line-noise
    signal_noise_before(idx_channel) = mean(sqrt(filtfilt(peak.b,peak.a,double(signal(:,idx_channel))).^2),1); %#ok<PFBNS>
    fprintf(1,'.');    
end
fprintf(1,'] done\n');

%% Common averaging 
% chan_group={datafile.ElectrodesInfo(:).ConnectorBank};
% ecog_chan_groups=chan_group(param.ecog_channels);
% ecog_chan_banks=unique(ecog_chan_groups);
% for chan_bk=1:size(ecog_chan_banks,2)
%     group_id=find(cellfun(@(x) strcmp(x,ecog_chan_banks(chan_bk)),ecog_chan_groups));
%     common_signal=signal(:,group_id);
%     common_mean=mean(common_signal,2);
%     signal_demean(:,group_id)=signal(:,group_id)-repmat(common_mean,1,size(group_id,2));
% end 

signal_demean=signal;

%% signal after first common average 
fprintf(1,'[');
parfor idx_channel=1:size(signal_demean,2),
    % calculate average root-mean-square of the line-noise
    signal_noise_after(idx_channel) = mean(sqrt(filtfilt(peak.b,peak.a,double(signal_demean(:,idx_channel))).^2),1); %#ok<PFBNS>
    fprintf(1,'.');    
end
fprintf(1,'] done\n');
param.ecog_channels_noise = param.ecog_channels( signal_noise_after > (mean(signal_noise_after)+1.5*std(signal_noise_after)));
param.channels_selected = setdiff(param.ecog_channels,param.ecog_channels_noise,'stable');

fprintf(1, '> Not using %d channels for common source averaging: ',length(param.ecog_channels_noise));
fprintf(1, '%d ',param.ecog_channels_noise);
fprintf(1, '\n');
%% common source averaging with exclusion of noisy channels 
%% Common averaging 

% ecog_chan_groups=chan_group(param.ecog_channels);
% ecog_chan_banks=unique(ecog_chan_groups);
% for chan_bk=1:size(ecog_chan_banks,2)
%     group_id=find(cellfun(@(x) strcmp(x,ecog_chan_banks(chan_bk)),ecog_chan_groups));
%     no_noise_group_id=setdiff(group_id,param.ecog_channels_noise,'stable');
% 
%     common_signal=signal(:,no_noise_group_id);
%     common_mean=mean(common_signal,2);
%     signal_demean(:,group_id)=signal(:,group_id)-repmat(common_mean,1,size(group_id,2));
% end
%% 
signal=signal_demean;
clear signal_demean
fprintf(1, '> Notch filtering signal \n');
fprintf(1,'[');

% for each channel
parfor idx_channel=1:size(signal,2),
    
    % get the signal for this channel
    signal_preliminary = double(signal(:,idx_channel));
    
    % remove all harmonics of line-noise
    for idx = 1:length(param.filter.notch.fcenter), %#ok<PFBNS>
        signal_preliminary = filtfilt(notch{idx}.b,notch{idx}.a,signal_preliminary); %#ok<PFBNS>
    end 
    
    % return the signal
    signal(:,idx_channel) = (signal_preliminary);
    
    fprintf(1,'.');
end
fprintf(1,'] done\n');
%% 
Cbrewr_map={'Blues','Reds','BuGn','PuBuGn','Oranges','BuPu','GnBu','Greens','Greys',...
    'OrRd','PuBu','PuRd','RdPu','Purples','YlGn',...
    'YlGnBu','YlOrRd','YlOrBr','Blues','Reds','BuGn','PuBuGn','Oranges','BuPu','GnBu','Greens','Greys'};
color_groups=[];
for chan_bk=1:size(ecog_tags,2)
    group_id=find(strcmp(ch_tag,ecog_tags{chan_bk}));
    %colorgroup=flipud(cbrewer('seq',Cbrewr_map{chan_bk},length(group_id)+20));
    colorgroup=flipud(inferno(length(group_id)+20));
    color_groups=[color_groups;colorgroup(1:length(group_id),:)];
end 
%%
%% plot the dataset throughout 
% get sections 
t_length=25e4;
kk_total=ceil(size(signal,1)/t_length);
x_cell=mat2cell(signal',ones(1,size(signal,2)));
valid_channels=param.ecog_channels*nan;
[a,b]=setdiff(param.ecog_channels,param.ecog_channels_noise,'stable');
valid_channels(b)=1;

x_cell=arrayfun(@(x) x_cell{x}*valid_channels(x),1:size(signal,2),'UniformOutput',false);
x_cell=x_cell';
min_max=nanmean(cell2mat(cellfun(@(y) [min(y),max(y)],x_cell,'UniformOutput',false)));
x_norm_cell=cellfun(@(x) (x-min_max(1))./(min_max(2)-min_max(1)),x_cell,'UniformOutput',false);
% 
figure(1);
%colors=inferno(size(x_norm_cell,1));
colors=color_groups;
for kk=1:kk_total
    clf;
    set(gcf,'color',[.7,.7,.7],'position',[17,1,2454,1337])
    ax=axes('position',[.02,.02,.95,.95]);
    hold on
    session_end=min([(kk*t_length),size(signal,1)]);
    t_window=((kk-1)*t_length+1):2:(session_end);
    
    arrayfun(@(x) plot(t_window,x_norm_cell{x}(t_window)+x-.5,'color',colors(x,:)),[1:size(x_norm_cell,1)])
    
    set(ax,'ytick',[1:size(x_norm_cell,1)]);
    set(ax,'yticklabel','')
    arrayfun(@(x) text(t_window(1),x,num2str(param.ecog_channels(x)),...
        'Color',colors(x,:),'HorizontalAlignment','right','VerticalAlignment','middle','fontsize',10),[1:size(x_norm_cell,1)]);
     arrayfun(@(x) text(t_window(end),x,param.ecog_channels_labels{x},...
        'Color',colors(x,:),'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',10),[1:size(x_norm_cell,1)]);
    
    set(ax,'ylim',[0,size(x_norm_cell,1)+1])
    set(ax,'TickLength',[0.001 0.0015])
    shg;
    title({'post common source averaging',[num2str(kk),' of ', num2str(kk_total)]})
    pause
end 
close all; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILTER OUT REMAINING LINE-NOISE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

prompt='additional channels to remove? format :[1,2] ';
x=input(prompt);
if ~isempty(x)
    param.ecog_channels_deselect=sort([param.ecog_channels_noise;(x)']);
    param.ecog_channels_selected = setdiff(param.ecog_channels,param.ecog_channels_deselect,'stable');
    t_length=25e4;
    kk_total=ceil(size(signal,1)/t_length);
    x_cell=mat2cell(signal',ones(1,size(signal,2)));
    valid_channels=param.ecog_channels*nan;
    [a,b]=setdiff(param.ecog_channels,param.ecog_channels_deselect,'stable');
    valid_channels(b)=1;

    x_cell=arrayfun(@(x) x_cell{x}*valid_channels(x),1:size(signal,2),'UniformOutput',false);
    x_cell=x_cell';
    min_max=nanmean(cell2mat(cellfun(@(y) [min(y),max(y)],x_cell,'UniformOutput',false)));
    x_norm_cell=cellfun(@(x) (x-min_max(1))./(.5*(min_max(2)-min_max(1))),x_cell,'UniformOutput',false);
    %
    figure(1);

    for kk=1:kk_total
        clf;
        set(gcf,'color',[.7,.7,.7],'position',[17,1,2454,1337])
        ax=axes('position',[.02,.02,.95,.95]);
        hold on
        session_end=min([(kk*t_length),size(signal,1)]);
        t_window=((kk-1)*t_length+1):2:(session_end);
        arrayfun(@(x) plot(t_window,x_norm_cell{x}(t_window)+x-1,'color',colors(x,:)),find(valid_channels==1))
        
        set(ax,'ytick',[1:size(x_norm_cell,1)]);
        set(ax,'yticklabel','')
        arrayfun(@(x) text(t_window(1),x,num2str(param.ecog_channels(x)),...
        'Color',colors(x,:),'HorizontalAlignment','right','VerticalAlignment','middle','fontsize',10),[1:size(x_norm_cell,1)]);
     arrayfun(@(x) text(t_window(end),x,param.ecog_channels_labels{x},...
        'Color',colors(x,:),'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',10),[1:size(x_norm_cell,1)]);
    
        set(ax,'ylim',[0,size(x_norm_cell,1)+1]);
        set(ax,'TickLength',[0.001 0.0015])
        shg;
        title({'post common source averaging',[num2str(kk),' of ', num2str(kk_total)]})
        pause
    end
end
close all; 
param.ecog_valid_chan=valid_channels;
%% create the bipolar version of the signal before doing hilbert transform 

%% modify the info based on sorted channel label information 
param.ecog_tags=ecog_tags;
param.chan_id_for_ecog_tag_grouped=ecog_loc_for_tag;
param.chan_name_for_ecog_tag_grouped=cellfun(@(x) ch_names(x),ecog_loc_for_tag,'uni',false);
reordered_chan_ids=cell2mat(ecog_loc_for_tag);
assert(issorted(reordered_chan_ids))
param.chan_id_for_ecog_tag=reordered_chan_ids;
param.chan_name_for_ecog_tags=ch_names(reordered_chan_ids);
% filter out non ecog channels 
all_valid_ch=find(param.ecog_valid_chan==1)';
valid_ecog_chan=intersect(all_valid_ch,reordered_chan_ids,'stable');
param.valid_ecog_chan_tag=param.ecog_valid_chan'*nan;
param.valid_ecog_chan_tag(valid_ecog_chan)=1;
valid_chan_id_for_ecog_tag=cellfun(@(x) intersect(valid_ecog_chan,x),ecog_loc_for_tag,'uni',false)';
param.chan_id_for_valid_ecog_tag=cell2mat(valid_chan_id_for_ecog_tag);
param.chan_id_for_valid_ecog_tag_grouped=(valid_chan_id_for_ecog_tag);
param.chan_name_for_valid_ecog_tags=ch_names(cell2mat(valid_chan_id_for_ecog_tag));
param.chan_name_for_valid_ecog_tags_grouped=cellfun(@(x) ch_names(x),valid_chan_id_for_ecog_tag,'uni',false);
valid_bipolar_diffs=cellfun(@(x) [x(find(diff(x)==1)),x(find(diff(x)==1)+1)],valid_chan_id_for_ecog_tag,'uni',false);
param.chan_id_valid_bipolar_diffs_grouped=valid_bipolar_diffs;
reordered_chan_diff_ids=cell2mat(valid_bipolar_diffs);
param.chan_id_valid_bipolar_diffs=reordered_chan_diff_ids;
param.ch_name_for_valid_bipolar_diffs=ch_names(cell2mat(valid_bipolar_diffs));
param.ch_name_for_valid_bipolar_diffs_grouped=cellfun(@(x)ch_names(x), valid_bipolar_diffs,'uni',false);

%% create a biopolar version of the signal 
fprintf(1, '>> creating a bipolar version of the signal:  \n');
    
fprintf(1,'[');
signal_bipolar=[];
for bipolar_id=1:size(reordered_chan_diff_ids,1)
    bipol_ch_1=reordered_chan_diff_ids(bipolar_id,1);
    bipol_ch_2=reordered_chan_diff_ids(bipolar_id,2);
    signal_bipolar(:,bipolar_id)=signal(:,bipol_ch_1)-signal(:,bipol_ch_2);
end 


%% plot the dataset throughout 
% get sections 
t_length=25e4;
kk_total=ceil(size(signal_bipolar,1)/t_length);
x_cell=mat2cell(signal_bipolar',ones(1,size(signal_bipolar,2)));
valid_channels=ones(size(signal_bipolar,2),1);

x_cell=arrayfun(@(x) x_cell{x}*valid_channels(x),1:size(signal_bipolar,2),'UniformOutput',false);
x_cell=x_cell';
min_max=nanmean(cell2mat(cellfun(@(y) [min(y),max(y)],x_cell,'UniformOutput',false)));
x_norm_cell=cellfun(@(x) (x-min_max(1))./(min_max(2)-min_max(1)),x_cell,'UniformOutput',false);
% 
figure(1);
%colors=inferno(size(x_norm_cell,1));
colors=color_groups;
for kk=1:kk_total
    clf;
    set(gcf,'color',[.7,.7,.7],'position',[17,1,2454,1337])
    ax=axes('position',[.02,.02,.95,.95]);
    hold on
    session_end=min([(kk*t_length),size(signal_bipolar,1)]);
    t_window=((kk-1)*t_length+1):2:(session_end);
    
    arrayfun(@(x) plot(t_window,x_norm_cell{x}(t_window)+x-.5,'color',colors(x,:)),[1:size(x_norm_cell,1)])
    
    set(ax,'ytick',[1:size(x_norm_cell,1)]);
    set(ax,'yticklabel','')
    arrayfun(@(x) text(t_window(1),x,num2str(param.ecog_channels(x)),...
        'Color',colors(x,:),'HorizontalAlignment','right','VerticalAlignment','middle','fontsize',10),[1:size(x_norm_cell,1)]);
     arrayfun(@(x) text(t_window(end),x,sprintf('%s-%s',param.ch_name_for_valid_bipolar_diffs{x,1},param.ch_name_for_valid_bipolar_diffs{x,2}),...
        'Color',colors(x,:),'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',10),[1:size(x_norm_cell,1)]);
    
    set(ax,'ylim',[0,size(x_norm_cell,1)+1])
    set(ax,'TickLength',[0.001 0.0015])
    shg;
    title({'post common source averaging',[num2str(kk),' of ', num2str(kk_total)]})
    pause
end 
close all; 
%% high gamma based on chang-lab 
f_gamma_low=70;
f_gamma_high=150;
[cfs,sds]=get_filter_param_chang_lab(f_gamma_low,f_gamma_high);
% add a correction to stay away from 60 hz 
cfs(1)=73.0;
% create a filter bank;


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRACT AND DECIMATE BANDPOWER FOR EACH FREQUENCY BAND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute high gamma based on chang-lab method, see : 
    % Dichter, Benjamin K., Jonathan D. Breshears, Matthew K. Leonard, and Edward F. Chang. 2018.
    % ?The Control of Vocal Pitch in Human Laryngeal Motor Cortex.? Cell 174 (1): 21?31.e9.
fprintf(1, '>> Extracting high gamma envelope based on chang-lab method  \n');
    
fprintf(1,'[');
for idx_channel=1:size(signal,2)
    
    filter_bank={};
    for s=1:length(cfs)
        filter_bank{s}=gaussian_filter(transpose(signal(:,idx_channel)),parameters.SamplingRate,cfs(s),sds(s));
    end
    signal_hilbert_bp = cellfun(@abs,hilbert_transform((transpose(signal(:,idx_channel))),...
        parameters.SamplingRate,...
        filter_bank),'UniformOutput',false);
    % mean across
    signal_hilbert_bp=mean(cell2mat(signal_hilbert_bp),1);
    % zscore the signal with respect to entire experimental block
    % to do
    means = mean(signal_hilbert_bp);
    stds = std(signal_hilbert_bp);
    signal_hilbert_zs = (signal_hilbert_bp - means) / stds;
    
    %
    signal_hilbert_loop(:,idx_channel)=signal_hilbert_bp;
    signal_hilbert_zs_loop(:,idx_channel)=signal_hilbert_zs;
    %
    fprintf(1,'.');
end
fprintf(1,'] done\n');
    
     %% down-sample signal
fprintf(1, '>> Down-sampling signal \n');
[P,Q] = rat(500/parameters.SamplingRate);
     
     % calculate decimation factor for 12 Hz sampling rate
decimation_factor = floor(parameters.SamplingRate / ecog.param.samplingrate);
parameters.decimation_factor=decimation_factor;
     
signal_hilbert_zs_decimated_loop = zeros(size(decimate(double(signal_hilbert_zs_loop(:,1)),decimation_factor),1),size(signal_hilbert_zs_loop,2),'double');
signal_hilbert_decimated_loop = zeros(size(decimate(double(signal_hilbert_loop(:,1)),decimation_factor),1),size(signal_hilbert_loop,2),'double');
     %    %
     %signal_hilbert_pca_decimated_loop = zeros(size(decimate(double(signal_hilbert_zs_loop(:,1)),decimation_factor),1),size(signal_hilbert_zs_loop,2),'single');
     %signal_hilbert_pca_zs_decimated_loop = zeros(size(decimate(double(signal_hilbert_zs_loop(:,1)),decimation_factor),1),size(signal_hilbert_zs_loop,2),'single');
%    
fprintf(1,'[');
for idx_channel=1:size(signal_hilbert_zs_loop,2),
    fprintf('channel %d \n',idx_channel)
    %         % decimate band-power envelope to 12 Hz
    signal_hilbert_zs_decimated_loop(:,idx_channel) = (decimate(double(signal_hilbert_zs_loop(:,idx_channel)),decimation_factor));
    signal_hilbert_decimated_loop(:,idx_channel) = (decimate(double(signal_hilbert_loop(:,idx_channel)),decimation_factor));   
    %signal_hilbert_pca_zs_decimated_loop(:,idx_channel) = (decimate(double(signal_hilbert_pca_zs_loop(:,idx_channel)),decimation_factor));
    %signal_hilbert_pca_decimated_loop(:,idx_channel) = (decimate(double(signal_hilbert_pca_loop(:,idx_channel)),decimation_factor));
    fprintf(1,'.');
end
fprintf(1,'] done\n');

%% high gamma filtering for bipolar version of the signal 

fprintf(1, '>> Extracting high gamma envelope based on chang-lab method  \n');
signal_bipolar_hilbert_loop=[];
signal_bipolar_hilbert_zs_loop=[];

fprintf(1,'[');
for idx_channel=1:size(signal_bipolar,2)
    
    filter_bank={};
    for s=1:length(cfs)
        filter_bank{s}=gaussian_filter(transpose(signal_bipolar(:,idx_channel)),parameters.SamplingRate,cfs(s),sds(s));
    end
    signal_bipolar_hilbert_bp = cellfun(@abs,hilbert_transform((transpose(signal_bipolar(:,idx_channel))),...
        parameters.SamplingRate,...
        filter_bank),'UniformOutput',false);
    % mean across
    signal_bipolar_hilbert_bp=mean(cell2mat(signal_bipolar_hilbert_bp),1);
    % zscore the signal with respect to entire experimental block
    % to do
    means = mean(signal_bipolar_hilbert_bp);
    stds = std(signal_bipolar_hilbert_bp);
    signal_bipolar_hilbert_zs = (signal_bipolar_hilbert_bp - means) / stds;
    
    %
    signal_bipolar_hilbert_loop(:,idx_channel)=signal_bipolar_hilbert_bp;
    signal_bipolar_hilbert_zs_loop(:,idx_channel)=signal_bipolar_hilbert_zs;
    %
    fprintf(1,'.');
end
fprintf(1,'] done\n');


%% downsample the signal bipolar 

fprintf(1, '>> Down-sampling signal \n');

     
     % calculate decimation factor for 12 Hz sampling rate
     
signal_bipolar_hilbert_zs_decimated_loop = zeros(size(decimate(double(signal_bipolar_hilbert_zs_loop(:,1)),decimation_factor),1),size(signal_bipolar_hilbert_zs_loop,2),'double');
signal_bioplar_hilbert_decimated_loop = zeros(size(decimate(double(signal_bipolar_hilbert_loop(:,1)),decimation_factor),1),size(signal_bipolar_hilbert_zs_loop,2),'double');
     %    %
     %signal_hilbert_pca_decimated_loop = zeros(size(decimate(double(signal_hilbert_zs_loop(:,1)),decimation_factor),1),size(signal_hilbert_zs_loop,2),'single');
     %signal_hilbert_pca_zs_decimated_loop = zeros(size(decimate(double(signal_hilbert_zs_loop(:,1)),decimation_factor),1),size(signal_hilbert_zs_loop,2),'single');
%    
fprintf(1,'[');
for idx_channel=1:size(signal_bipolar_hilbert_zs_loop,2),
    fprintf('channel %d \n',idx_channel)
    %         % decimate band-power envelope to 12 Hz
    signal_bipolar_hilbert_zs_decimated_loop(:,idx_channel) = (decimate(double(signal_bipolar_hilbert_zs_loop(:,idx_channel)),decimation_factor));
    signal_bioplar_hilbert_decimated_loop(:,idx_channel) = (decimate(double(signal_bipolar_hilbert_loop(:,idx_channel)),decimation_factor));   
    %signal_hilbert_pca_zs_decimated_loop(:,idx_channel) = (decimate(double(signal_hilbert_pca_zs_loop(:,idx_channel)),decimation_factor));
    %signal_hilbert_pca_decimated_loop(:,idx_channel) = (decimate(double(signal_hilbert_pca_loop(:,idx_channel)),decimation_factor));
    fprintf(1,'.');
end
fprintf(1,'] done\n')

 
%% save data out 
dataout=datafile;
dataout.signal_hilbert_zs_decimated=signal_hilbert_zs_decimated_loop';
dataout.signal_hilbert_decimated=signal_hilbert_decimated_loop';
dataout.signal_bipolar_hilbert_zs_decimated=signal_bipolar_hilbert_zs_decimated_loop';
dataout.signal_bioplar_hilbert_decimated=signal_bioplar_hilbert_decimated_loop';
% 
dataout.signal_hilbert_zs=signal_hilbert_zs_loop';
dataout.signal_hilbert=signal_hilbert_loop';
dataout.signal_bipolar_hilbert_zs=signal_bipolar_hilbert_zs_loop';
dataout.signal_bioplar_hilbert=signal_bipolar_hilbert_loop';


dataout=rmfield(dataout,'Data');

%dataout.filtered_signal=(signal');
dataout.MetaTags.DownSampleFrequency=parameters.SamplingRate/decimation_factor;
param_fields=fields(param);
for k=1:length(param_fields)
    dataout.MetaTags.(param_fields{k})=param.(param_fields{k});
end 
subj_file_name=strcat(dataout.MetaTags.FilePath,filesep,...
    dataout.MetaTags.Filename,'_filtered_v2.mat');
save(subj_file_name,'dataout','-v7.3');

% create a shorter copy 
dataout=rmfield(dataout,'signal_hilbert');
dataout=rmfield(dataout,'signal_hilbert_zs');
dataout=rmfield(dataout,'signal_bipolar_hilbert_zs');
dataout=rmfield(dataout,'signal_bioplar_hilbert');

subj_file_name=strcat(dataout.MetaTags.FilePath,filesep,...
    dataout.MetaTags.Filename,'_filtered_compressed_v2.mat');
save(subj_file_name,'dataout','-v7.3');


end 

%% 
end 