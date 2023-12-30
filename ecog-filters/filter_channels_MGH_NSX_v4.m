function [dataout]=filter_channels_MGH_NSX_v4(datafile,varargin)
% varargin :
%1 subject operation info


p=inputParser();
addParameter(p, 'op_info', struct);
addParameter(p, 'aux_data', []);
addParameter(p, 'elecFromAux', true);
addParameter(p, 'timing', []);
addParameter(p,'do_IED',false)
addParameter(p,'do_globelMean',false)
addParameter(p,'do_shankMean',false)
parse(p, varargin{:});
ops = p.Results;

P = preprocessing_parameters


ecog.param.filter_car  = 1;     % 0 == off % 1 == on
ecog.param.filter_type = 'IIR'; % IIR | FIR

% --- highpass filter ---
ecog.param.highpass.Wp = 0.50; % Hz
ecog.param.highpass.Ws = 0.05; % Hz
ecog.param.highpass.Rp = 3;    % dB
ecog.param.highpass.Rs = 30;   % dB
param=ecog.param;
%---- downsampling rate ----
param.downsamplingrate = 128; % decimation sampling rate
param.sig_sampling_rate=datafile.MetaTags.SamplingFreq;
param.ecog_channels_labels=datafile.ElectrodesInfo.Label';
param.ecog_channels=[1:length(param.ecog_channels_labels)]';
param.preproc_step={};
% bipolar versions
%% do some assertion if the electrodes are selected from aux dataset
search_param=false;
no_aux=false;
if ~isempty(ops.aux_data)
    try
        assert(all(cellfun(@(x,y) strcmp(x,erase(y,'_')),param.ecog_channels_labels ,ops.aux_data.elec_ch_label)))
        for pp=1:size(param.ecog_channels_labels,1)
            assert(strcmp(param.ecog_channels_labels{pp},erase(ops.aux_data.elec_ch_label{pp},'_')))
        end
    catch err
        warning('difference electrode configuration between aux data and datafile; using search to find out electrode correspondance')
        search_param=true;
    end
    assert(param.sig_sampling_rate==ops.aux_data.info.sig_sampling_rate)
else
    no_aux=true;
end
%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE HIGH-PASS FILTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define passband, stopband and attenuation
highpass{1}.Wp = param.highpass.Wp/(param.sig_sampling_rate/2);
highpass{1}.Ws = param.highpass.Ws/(param.sig_sampling_rate/2);
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
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE LINE-NOISE FILTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
peak.fcenter = 60;
peak.bw      = 0.001;

% calculate the IIR-peak filter coefficients in a,b format
peak.wo = peak.fcenter/(param.sig_sampling_rate/2);
peak.bw = peak.bw;
[peak.b,peak.a] = iirpeak(peak.wo,peak.bw);

% define the harmonics of line noise frequency
param.filter.notch.fcenter = [60,120,180,240];
param.filter.notch.bw      = ones(1,length(param.filter.notch.fcenter)).*0.001;

% calculate the IIR-peak filter coefficients in a,b format
for idx = 1:length(param.filter.notch.fcenter),
    notch{idx}.wo = param.filter.notch.fcenter(idx)/(param.sig_sampling_rate/2);
    notch{idx}.bw = param.filter.notch.bw(idx);
    [notch{idx}.b,notch{idx}.a] = iirnotch(notch{idx}.wo,notch{idx}.bw);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD AND CONCATENATE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signal = double([]);
num_of_chan=numel(param.ecog_channels);
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
    signal(:,idx) = (filtfilt(highpass{1}.sos,highpass{1}.g,double(signal_loop(:,idx_channel))));
    fprintf(1,'.');
end
fprintf(1,'] done\n');

P.cutoff=P.hp_filt_cutoff_in_Hz;

signal_1 = hp_filt(signal_loop,P , param.sig_sampling_rate);

clear signal_loop
param.preproc_step=[param.preproc_step;'high_pass'];

%% 

%% test epileptic spike detection
if (~isempty(ops.aux_data) & ops.elecFromAux)
    if ~search_param
        param.ecog_channels_IED_deselected=ops.aux_data.info.ecog_channels_IED_deselected;
        param.IED_results=ops.aux_data.info.IED_results;
    else
        aux_chan_deselect=ops.aux_data.info.ecog_channels_labels(ops.aux_data.info.ecog_channels_IED_deselected)
        where_in_datafile=cell2mat(cellfun(@(x) find(ismember(param.ecog_channels_labels,erase(x,'_'))),aux_chan_deselect,'uni',false));
        param.ecog_channels_IED_deselected=where_in_datafile;
        param.IED_results=ops.aux_data.info.IED_results;
    end
    
else
    detectionIEDs          = []; % output from Janca et al. script - 1 cell array per segment
    detectionIEDs.settings = '-k1 3.65 -h 60 -dec 200 -dt 0.005 -pt 0.12 -ti 1'; % if you change "-dec 200" here, do not forget to change in selectChannels_Using ... below
    detectionIEDs.segments = [];
    
    detectionIEDs=automaticSpikeDetection_UsingJancaMethod(signal, param.sig_sampling_rate, detectionIEDs.settings);
    
    % plot the detections
    % plot IED events
    t_len=100;
    t_length=t_len*200; % 20 sec * 200 downsample rate
    kk_total=ceil(size(detectionIEDs.envelope,1)/t_length);
    x_cell=mat2cell(detectionIEDs.envelope',ones(1,size(detectionIEDs.envelope,2)));
    min_max=mean(cell2mat(cellfun(@(y) prctile(y,[1 99]),x_cell,'UniformOutput',false)));
    x_norm_cell=cellfun(@(x) (x-min_max(1))./(min_max(2)-min_max(1)),x_cell,'UniformOutput',false);
    
    
    col_inf=inferno(floor(.8*size(x_norm_cell,1)));
    col_vir=viridis(floor(.8*size(x_norm_cell,1)));
    colors=[col_vir(1:floor(size(x_norm_cell,1)/2),:);col_inf(1:(floor(size(x_norm_cell,1)/2)+1),:)];
    figure(1);
    clf;
    set(gcf,'position',[31,1,1713,1010]);
    ax=axes('position',[.05,.1,.93,.88]);
    hold on
    time_stamps=[1:size(detectionIEDs.envelope,1)]/200;
    hold on
    H=arrayfun(@(x) plot(time_stamps,x_norm_cell{x}+x,'color',colors(x,:),'tag',sprintf('ch %d, tag %s',x,param.ecog_channels_labels{x})),[1:size(x_norm_cell,1)]);
    %H=arrayfun(@(x) plot(time_stamps,x_norm_cell{x}+x,'color',colors(x,:),'tag',sprintf('ch %d, tag %s',x,param.ecog_channels_labels{x})),[1:3]);
    spike_chan = detectionIEDs.out.chan;
    [spike_chan_sort,sort_idx]=sort(spike_chan);
    spike_times_sort= detectionIEDs.out.pos(sort_idx);
    
    %spike_dur = detectionIEDs.out.dur;
    yval=cell2mat(arrayfun(@(x) x_norm_cell{spike_chan_sort(x)}(floor(spike_times_sort(x)*200))+spike_chan_sort(x),1:length(spike_chan_sort),'uni',false));
    H1=scatter(spike_times_sort,yval,50,'MarkerEdgeColor',[0 0 0],...
        'MarkerFaceColor',[0 .7 .7],...
        'LineWidth',2,'MarkerFaceAlpha',.5);
    
    set(gcf,'doublebuffer','on');
    set(ax,'ytick',[1:size(x_norm_cell,1)]);
    set(ax,'yticklabel','');
    %set(ax,'ylim',[0,10]);
    set(ax,'ylim',[0,size(x_norm_cell,1)+4]);
    ax.XAxis.TickLength=[0.005,0.01];
    ax.YAxis.TickLength=[0.005,0.01];
    set(ax,'xlim',[0 t_len]);
    pos=get(ax,'position');
    Newpos=[pos(1) pos(2)-0.1 pos(3) 0.05];
    xmax=max(time_stamps);
    S=['set(gca,''xlim'',get(gcbo,''value'')+[0 ' num2str(t_len) '])'];
    
    h=uicontrol('style','slider',...
        'units','normalized','position',Newpos,...
        'callback',S,'min',0,'max',xmax-t_len);
    
    datacursormode on
    dcm = datacursormode(gcf);
    set(dcm,'UpdateFcn',@myupdatefcn)
    
    waitfor(findobj('type','figure','number',1));
    
    %
    
    % From this automatic assessment, selection of the final pool of channels
    detectionIEDs.tableChanSelection = [];  % info regarding chan selection
    detectionIEDs.threshold          = 6.5; % channels with IEDs higher than threshold are removed
    
    
    
    currIEDs.fs             = 200; % default downsampling during automatic detection (-dec 200)
    currIEDs.discharges.MV  = [];
    currIEDs.numSamples     = 0;
    
    
    currIEDs.discharges.MV = detectionIEDs.discharges.MV;
    currIEDs.numSamples= currIEDs.numSamples + size(detectionIEDs.d_decim, 1);
    
    
    % Compute number of detected spike per channel: [c x 1] where c channels
    numSpikes     = []; numSpikes     = sum(currIEDs.discharges.MV==1, 1);
    totalDuration = []; totalDuration = (currIEDs.numSamples / currIEDs.fs) / 60; % in minutes
    numSpikes_min = []; numSpikes_min = numSpikes / totalDuration;
    numSpikes     = transpose(numSpikes);
    numSpikes_min = transpose(numSpikes_min);
    
    % Select channels with IEDs / minute below threshold - [c x 1] where c channels
    indChanSelected = [];
    indChanSelected = find(numSpikes_min < detectionIEDs.threshold);
    tableChanSelection.numSpikesAll           = numSpikes_min;
    tableChanSelection.indChansSelected       = indChanSelected;
    tableChanSelection.indChansDeselected     = setdiff(param.ecog_channels,indChanSelected);
    tableChanSelection.nameChansSelected      = transpose(param.ecog_channels_labels(indChanSelected));
    tableChanSelection.numSpikesChansSelected = numSpikes_min(indChanSelected);
    
    
    param.IED_results=tableChanSelection;
    param.IED_results.threshold=detectionIEDs.threshold;
    %
    if ops.do_IED
        param.ecog_channels_IED_deselected=tableChanSelection.indChansDeselected;
        param.preproc_step=[param.preproc_step;'IED_removal'];
    else 
        param.ecog_channels_IED_deselected=[];
    end 

end
%%
[good_channels,bad_channels, rms60Hz] = ...
    channel_selection_from_60Hz_noise(signal, param.sig_sampling_rate)
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


%%
fprintf(1, '> Notch filtering signal \n');
fprintf(1,'[');
% for each channel
parfor idx_channel=1:size(signal,2)
    
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


fprintf(1, '> Meassuring 60 Hz noise power after signal processing \n');
fprintf(1,'[');
signal_noise_after=[];
parfor idx_channel=1:size(signal,2),
    % calculate average root-mean-square of the line-noise
    signal_noise_after(idx_channel) = mean(sqrt(filtfilt(peak.b,peak.a,double(signal(:,idx_channel))).^2),1); %#ok<PFBNS>
    fprintf(1,'.');
end
fprintf(1,'] done\n');
param.preproc_step=[param.preproc_step;'removed_line_noise'];
%%

if (~isempty(ops.aux_data) & ops.elecFromAux)
    if ~search_param
        param.ecog_channels_noise_5std_deselected=ops.aux_data.info.ecog_channels_noise_5std_deselected;
        param.ecog_channels_selected=ops.aux_data.info.ecog_channels_selected;
    else
        aux_chan_deselect=ops.aux_data.info.ecog_channels_labels(ops.aux_data.info.ecog_channels_noise_5std_deselected)
        where_in_datafile=cell2mat(cellfun(@(x) find(ismember(param.ecog_channels_labels,erase(x,'_'))),aux_chan_deselect,'uni',false));
        param.ecog_channels_noise_5std_deselected=where_in_datafile;
        param.ecog_channels_selected = setdiff(param.ecog_channels,union(param.ecog_channels_noise_5std_deselected,param.ecog_channels_IED_deselected),'stable');
        
    end
else
    
    param.ecog_channels_noise_5std_deselected = param.ecog_channels( signal_noise_after > (mean(signal_noise_after)+5*std(signal_noise_after)));
    param.ecog_channels_selected = setdiff(param.ecog_channels,union(param.ecog_channels_noise_5std_deselected,param.ecog_channels_IED_deselected),'stable');
end
%
valid_channels=param.ecog_channels*nan;
valid_channels(param.ecog_channels_selected)=1;

%% 
col_inf=inferno(floor(.8*size(signal,2)));
col_vir=viridis(floor(.8*size(signal,2)));
colors=[col_vir(1:floor(size(signal,2)/2),:);col_inf(1:(floor(size(signal,2)/2)+1),:)];
plot_channels(signal,param.sig_sampling_rate,colors,param.ecog_channels_labels,valid_channels)
waitfor(findobj('type','figure','number',1));

%%
if (~isempty(ops.aux_data) & ops.elecFromAux)
    param.ecog_valid_chan=valid_channels;
else
    prompt='additional channels to remove? format :[1,2] ';
    x=input(prompt);
    if ~isempty(x)
        param.ecog_channels_user_deselect=x;
        param.ecog_channels_deselect=sort(union(union(param.ecog_channels_noise_5std_deselected,param.ecog_channels_IED_deselected,'stable'),x,'stable'));
        param.ecog_channels_selected = setdiff(param.ecog_channels_selected,param.ecog_channels_deselect,'stable');
        
        valid_channels=param.ecog_channels*nan;
        valid_channels(param.ecog_channels_selected)=1;
    else
        param.ecog_channels_user_deselect=[];
    end
    
    param.ecog_valid_chan=valid_channels;
end
param.preproc_step=[param.preproc_step;'ch_removal_by_visual_inspection'];

%% do global mean removal (not sure if we should do this!)
if ops.do_globelMean
    overall_mean=mean(signal(:,param.ecog_channels),2);
    %overall_mean=mean(signal,2);
    %figure;plot(overall_mean);
    signal_glob=signal-repmat(overall_mean,1,size(signal,2));
    signal=signal_glob;
    clear signal_glob;
    param.preproc_step=[param.preproc_step;'global_mean_removal'];
end

%% do common source averaging based on tags (not sure if this should be done in this case)
op_ecog_labels=param.ecog_channels_labels;
ch_num=cellfun(@str2num,extract(op_ecog_labels,digitsPattern));
ch_tags=extract(op_ecog_labels,lettersPattern);
[ecog_tags,~,~]=unique(ch_tags,'stable');
ecog_loc_for_tag=cellfun(@(x) find(contains(op_ecog_labels,x))', ecog_tags,'uni',false);

% 
signal_ca=nan*zeros(size(signal));
for kk=1:size(ecog_loc_for_tag,1)
    same_shank=ecog_loc_for_tag{kk};
    same_shank_clean=intersect(same_shank,param.ecog_channels_selected);
    if ~isempty(same_shank_clean)
        com_ave=mean(signal(:,same_shank_clean),2);
        %com_ave_glob=mean(signal_glob(:,same_shank_clean),2);
    else 
        warning('no clean electrode for common source average')
        com_ave=zeros(size(signal,1),1);
        %com_ave_glob=zeros(size(signal,1),1);
    end 
    signal_ca(:,same_shank)=signal(:,same_shank)-com_ave;
    %signal_ca_glob(:,same_shank)=signal_glob(:,same_shank)-com_ave_glob;
end
if ops.do_shankMean
    signal=signal_ca;
    param.preproc_step=[param.preproc_step;'shank_common_source_removal'];
end 

%% 
color_sequence={{'Blues'   },{'YlOrRd'  },{'BuGn'    },{'Greys'   },{'BuPu'    },{'GnBu'    },{'Oranges' },{'Greens'  },{'OrRd'    },{'Blues'   },{'PuBu'    },{'PuBuGn'  },...
{'PuRd'    },{'Blues'   },{'Purples' },{'Greens'  },{'RdPu'    },{'YlGn'    },{'Reds'    },{'YlGnBu'  },{'YlOrBr'  },{'YlOrRd'  },...
{'BuGn'    },{'BuPu'    },{'GnBu'    },{'Greys'   },{'Oranges' },{'OrRd'    },{'PuBu'    },{'PuBuGn'  },{'PuRd'    },{'Purples' },{'RdPu'    },{'Reds'    },{'YlGn'    },{'YlGnBu'  },...
{'YlOrBr'  }};
%color_sequence=color_sequence(randperm(length(color_sequence)));
colors_map=(arrayfun(@(x) flipud(cbrewer('seq',color_sequence{x}{1},floor(1.5*length(ecog_loc_for_tag{x})),'linear')),[1:size(ecog_loc_for_tag,1)]','uni',false));
colors=cell2mat(cellfun(@(x,y) y(1:size(x,2),:),ecog_loc_for_tag,colors_map,'uni',false));

plot_channels(signal,param.sig_sampling_rate,colors,param.ecog_channels_labels,valid_channels)
waitfor(findobj('type','figure','number',1));

%% modify the info based on sorted channel label information
if (~isempty(ops.aux_data) & ops.elecFromAux)
    if ~search_param
        param.chan_id_valid_bipolar_diffs_grouped=ops.aux_data.info.chan_id_valid_bipolar_diffs_grouped;
        param.chan_id_valid_bipolar_diffs=ops.aux_data.info.chan_id_valid_bipolar_diffs;
        param.ch_name_for_valid_bipolar_diffs=ops.aux_data.info.ch_name_for_valid_bipolar_diffs;
        param.ch_name_for_valid_bipolar_diffs_grouped=ops.aux_data.info.ch_name_for_valid_bipolar_diffs_grouped;
        reordered_chan_diff_ids=param.chan_id_valid_bipolar_diffs;
        op_ecog_labels=param.ecog_channels_labels;
    else
        where_in_datafile=cell2mat(cellfun(@(x) sum(ismember(param.ecog_channels_labels,erase(x,'_'))==1),ops.aux_data.info.ch_name_for_valid_bipolar_diffs,'uni',false));
        % find rows that have empty slot
        pairs_to_maintain=find(sum(where_in_datafile,2)==2);
        param.ch_name_for_valid_bipolar_diffs=ops.aux_data.info.ch_name_for_valid_bipolar_diffs(pairs_to_maintain,:);
        param.chan_id_valid_bipolar_diffs=cell2mat(cellfun(@(x) find(ismember(param.ecog_channels_labels,erase(x,'_'))),param.ch_name_for_valid_bipolar_diffs,'uni',false));
        op_ecog_labels=param.ecog_channels_labels;
        ch_num=cellfun(@str2num,extract(op_ecog_labels,digitsPattern));
        ch_tags=extract(op_ecog_labels,lettersPattern);
        [ecog_tags,ia,ic]=unique(ch_tags,'stable');
        param.ch_name_for_valid_bipolar_diffs_grouped=cellfun(@(x) reshape(param.ch_name_for_valid_bipolar_diffs(contains(param.ch_name_for_valid_bipolar_diffs,x)),[],2), ecog_tags,'uni',false);
        param.chan_id_valid_bipolar_diffs_grouped=cellfun(@(x) reshape(param.chan_id_valid_bipolar_diffs(contains(param.ch_name_for_valid_bipolar_diffs,x)),[],2), ecog_tags,'uni',false);
        % do an assertion
        for tt=1:size(param.ch_name_for_valid_bipolar_diffs_grouped)
            xx=param.ch_name_for_valid_bipolar_diffs_grouped{tt};
            yy=cell2mat(cellfun(@(y) find(strcmp(param.ecog_channels_labels,y)),erase(xx,'_'),'uni',false));
            assert(isequal(yy,param.chan_id_valid_bipolar_diffs_grouped{tt}));
        end
        reordered_chan_diff_ids=param.chan_id_valid_bipolar_diffs;
        op_ecog_labels=param.ecog_channels_labels;
    end
else
    op_ecog_labels=param.ecog_channels_labels;
    ch_num=cellfun(@str2num,extract(op_ecog_labels,digitsPattern));
    ch_tags=extract(op_ecog_labels,lettersPattern);
    [ecog_tags,ia,ic]=unique(ch_tags,'stable');
    
    ecog_loc_for_tag=cellfun(@(x) find(contains(op_ecog_labels,x))', ecog_tags,'uni',false);
    % remove unvalid from ecog loc
    valid_chan_id_for_ecog_tag=cellfun(@(x) intersect(x,param.ecog_channels_selected,'stable'),ecog_loc_for_tag,'uni',false);
    valid_chan_num_for_ecog_tag=cellfun(@(x) ch_num(x),valid_chan_id_for_ecog_tag,'uni',false);
    valid_bipolar_diffs=cellfun(@(x,y) [y(find(diff(x)==1))+1,y(find(diff(x)==1))],valid_chan_num_for_ecog_tag,valid_chan_id_for_ecog_tag,'uni',false);
    param.chan_id_valid_bipolar_diffs_grouped=valid_bipolar_diffs';
    reordered_chan_diff_ids=cell2mat(param.chan_id_valid_bipolar_diffs_grouped');
    param.chan_id_valid_bipolar_diffs=reordered_chan_diff_ids;
    param.ch_name_for_valid_bipolar_diffs=op_ecog_labels(reordered_chan_diff_ids);
    param.ch_name_for_valid_bipolar_diffs_grouped=cellfun(@(x)op_ecog_labels(x), param.chan_id_valid_bipolar_diffs_grouped,'uni',false);
end


%% create a biopolar version of the signal
fprintf(1, '>> creating a bipolar version of the signal:  \n');
signal_bipolar= double([]);
%signal_bipolar_ca= double([]);
for bipolar_id=1:size(reordered_chan_diff_ids,1)
    bipol_ch_1=reordered_chan_diff_ids(bipolar_id,1);
    bipol_ch_2=reordered_chan_diff_ids(bipolar_id,2);
    bip_ch_name2=op_ecog_labels{bipol_ch_2};
    bip_ch_name1=op_ecog_labels{bipol_ch_1};
    B = cell2mat(extract(bip_ch_name2,lettersPattern));
    A =cell2mat(extract(bip_ch_name1,lettersPattern));
    assert(all(A==B));
    B =str2num(cell2mat(extract(bip_ch_name2,digitsPattern)));
    A =str2num(cell2mat(extract(bip_ch_name1,digitsPattern)));
    assert(A-B==1)
    signal_bipolar(:,bipolar_id)=signal(:,bipol_ch_1)-signal(:,bipol_ch_2);
    %signal_bipolar_ca(:,bipolar_id)=signal_ca(:,bipol_ch_1)-signal_ca(:,bipol_ch_2);
end


% if max(reshape(abs(signal_bipolar_ca-signal_bipolar),1,[]))>1e-5
%     warning('difference between bipolar and bipolar with common referencing');
% else
%     signal=signal_ca;
% end 
param.preproc_step=[param.preproc_step;'created_bipolar_signal'];
%% creat color map for bipolar 
param.chan_id_valid_bipolar_diffs;
colors_map_bip=(arrayfun(@(x) flipud(cbrewer('seq',color_sequence{x}{1},floor(1.5*size(param.chan_id_valid_bipolar_diffs_grouped{x},1)),'linear')),[1:max(size(param.chan_id_valid_bipolar_diffs_grouped))]','uni',false));
if size(param.chan_id_valid_bipolar_diffs_grouped,2)>1
    param.chan_id_valid_bipolar_diffs_grouped=param.chan_id_valid_bipolar_diffs_grouped';
end 
colors_bip=cell2mat(cellfun(@(x,y) y(1:size(x,1),:),param.chan_id_valid_bipolar_diffs_grouped,colors_map_bip,'uni',false));

bip_chan_labels=arrayfun(@(x) sprintf('%s-%s',param.ch_name_for_valid_bipolar_diffs{x,1},param.ch_name_for_valid_bipolar_diffs{x,2}),1:size(signal_bipolar,2),'uni',false);
bip_chan_labels=erase(bip_chan_labels,'_');


plot_channels(signal_bipolar,param.sig_sampling_rate,colors,bip_chan_labels,ones(size(bip_chan_labels,2)),'timing',ops.timing)
waitfor(findobj('type','figure','number',1));


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
% https://github.com/bendichter/process_ecog/blob/master/ecog/signal_processing/hilbert_transform.py
fprintf(1, '>> Extracting high gamma envelope based on chang-lab method  \n');

fprintf(1,'[');
for idx_channel=1:size(signal,2)
    
    filter_bank={};
    for s=1:length(cfs)
        filter_bank{s}=gaussian_filter(transpose(signal(:,idx_channel)),param.sig_sampling_rate,cfs(s),sds(s));
    end
    signal_hilbert_bp_all = cellfun(@abs,hilbert_transform((transpose(signal(:,idx_channel))),...
        param.sig_sampling_rate,...
        filter_bank),'UniformOutput',false);
    % mean across
    signal_hilbert_bp=mean(cell2mat(signal_hilbert_bp_all),1);
    % zscore the signal with respect to entire experimental block
    % to do
    means = mean(signal_hilbert_bp);
    stds = std(signal_hilbert_bp);
    zscore(signal_hilbert_bp);
    signal_hilbert_zs = (signal_hilbert_bp - means) / stds;
    assert(all(signal_hilbert_zs==zscore(signal_hilbert_bp)))
    
    
    %
    signal_hilbert_loop(:,idx_channel)=signal_hilbert_bp;
    signal_hilbert_zs_loop(:,idx_channel)=signal_hilbert_zs;
    %
    fprintf(1,'.\n');
end
fprintf(1,'] done\n');

param.preproc_step=[param.preproc_step;'extracting_unipolar_gamma_envelope'];
%% do some checks 
plot_channels(cell2mat(signal_hilbert_bp_all)',param.sig_sampling_rate,colors,param.ecog_channels_labels,1:8)

%% 

plot_channels(signal_hilbert_loop,param.sig_sampling_rate,colors,param.ecog_channels_labels,valid_channels)
waitfor(findobj('type','figure','number',1));
%%
down_sample_rate=8;
t_len=60;
t_length=t_len*param.sig_sampling_rate/down_sample_rate; % 20 sec * 200 downsample rate
x_cell=mat2cell(signal_hilbert_loop',ones(1,size(signal_hilbert_loop,2)));
x_cell=cellfun(@(x) downsample(x,down_sample_rate),x_cell,'uni',false);
min_max=mean(cell2mat(cellfun(@(y) prctile(y,[5 95]),x_cell,'UniformOutput',false)));
x_norm_cell=cellfun(@(x) (x-min_max(1))./(min_max(2)-min_max(1)),x_cell,'UniformOutput',false);

%col_inf=inferno(floor(.8*size(x_norm_cell,1)));
%col_vir=viridis(floor(.8*size(x_norm_cell,1)));
%colors=[col_vir(1:floor(size(x_norm_cell,1)/2),:);col_inf(1:(floor(size(x_norm_cell,1)/2)+1),:)];
figure(1);
clf;
set(gcf,'position',[31,1,1713,1010]);
ax=axes('position',[.05,.1,.93,.88]);
time_stamps=([1:size(x_cell{1},2)]/param.sig_sampling_rate)*down_sample_rate;
hold on
H=arrayfun(@(x) plot(time_stamps,x_norm_cell{x}+x,'color',colors(x,:),'linewidth',1.5,'tag',sprintf('ch %d, tag %s',x,param.ecog_channels_labels{x})),[1:size(x_norm_cell,1)]);

if ~isempty(ops.timing)
    timings_ds=cell2mat(cellfun(@(x) floor(x./param.sig_sampling_rate),ops.timing,'uni',false));
     for kk=1:length(timings_ds)
        event=timings_ds(kk,:);
        corners_x=[min(event),max(event),max(event),min(event)];
        corners_y=[ax.YLim(1),ax.YLim(1),ax.YLim(2),ax.YLim(2)];
        pat=patch(corners_x,corners_y,[1,.4,.4]);
        pat.FaceAlpha=.1;
    end 
end 
set(gcf,'doublebuffer','on');
set(ax,'ytick',[1:size(x_norm_cell,1)]);
set(ax,'yticklabel','');
set(ax,'children',flipud(get(gca,'children')))

%set(ax,'ylim',[0,10]);
set(ax,'ylim',[0,size(x_norm_cell,1)+4]);
ax.XAxis.TickLength=[0.005,0.01];
ax.YAxis.TickLength=[0.005,0.01];
set(ax,'xlim',[0 t_len]);
pos=get(ax,'position');
Newpos=[pos(1) pos(2)-0.1 pos(3) 0.05];
xmax=max(time_stamps);
S=['set(gca,''xlim'',get(gcbo,''value'')+[0 ' num2str(t_len) '])'];

h=uicontrol('style','slider',...
    'units','normalized','position',Newpos,...
    'callback',S,'min',0,'max',xmax-t_len);

datacursormode on
dcm = datacursormode(gcf);
set(dcm,'UpdateFcn',@myupdatefcn)

waitfor(findobj('type','figure','number',1));

%% down-sample signal
fprintf(1, '>> downsampling signal \n');


decimation_factor = floor(param.sig_sampling_rate / param.downsamplingrate);
param.decimation_factor=decimation_factor;

signal_hilbert_zs_decimated_loop = [];
signal_hilbert_decimated_loop = [];
%    %
fprintf(1,'[');
for idx_channel=1:size(signal_hilbert_zs_loop,2),
    fprintf('channel %d \n',idx_channel)
    signal_hilbert_zs_decimated_loop(:,idx_channel) = (decimate((signal_hilbert_zs_loop(:,idx_channel)),decimation_factor));
    signal_hilbert_decimated_loop(:,idx_channel) = (decimate((signal_hilbert_loop(:,idx_channel)),decimation_factor));
    fprintf(1,'.');
end
fprintf(1,'] done\n');
param.preproc_step=[param.preproc_step;'downsampling_unipolar_gamma_envelope'];
%% high gamma filtering for bipolar version of the signal

fprintf(1, '>> Extracting high gamma envelope based on chang-lab method  \n');
signal_bipolar_hilbert_loop=[];
signal_bipolar_hilbert_zs_loop=[];

fprintf(1,'[');
for idx_channel=1:size(signal_bipolar,2)
    
    filter_bank={};
    for s=1:length(cfs)
        filter_bank{s}=gaussian_filter(transpose(signal_bipolar(:,idx_channel)),param.sig_sampling_rate,cfs(s),sds(s));
    end
    signal_bipolar_hilbert_bp = cellfun(@abs,hilbert_transform((transpose(signal_bipolar(:,idx_channel))),...
        param.sig_sampling_rate,...
        filter_bank),'UniformOutput',false);
    % mean across
    signal_bipolar_hilbert_bp=mean(cell2mat(signal_bipolar_hilbert_bp),1);
    % zscore the signal with respect to entire experimental block
    % to do
    means = mean(signal_bipolar_hilbert_bp);
    stds = std(signal_bipolar_hilbert_bp);
    signal_bipolar_hilbert_zs = (signal_bipolar_hilbert_bp - means) / stds;
    assert(all(signal_bipolar_hilbert_zs==zscore(signal_bipolar_hilbert_bp)))
    %
    signal_bipolar_hilbert_loop(:,idx_channel)=signal_bipolar_hilbert_bp;
    signal_bipolar_hilbert_zs_loop(:,idx_channel)=signal_bipolar_hilbert_zs;
    %
    fprintf(1,'.\n');
end
fprintf(1,'] done\n');
param.preproc_step=[param.preproc_step;'extracting_bipolar_gamma_envelope'];
%% plot throughout 
down_sample_rate=8;
t_len=60;
t_length=t_len*param.sig_sampling_rate/down_sample_rate; % 20 sec * 200 downsample rate
x_cell=mat2cell(signal_bipolar_hilbert_loop',ones(1,size(signal_bipolar_hilbert_loop,2)));
x_cell=cellfun(@(x) downsample(x,down_sample_rate),x_cell,'uni',false);
min_max=mean(cell2mat(cellfun(@(y) prctile(y,[5 95]),x_cell,'UniformOutput',false)));
x_norm_cell=cellfun(@(x) (x-min_max(1))./(min_max(2)-min_max(1)),x_cell,'UniformOutput',false);

%col_inf=inferno(floor(.8*size(x_norm_cell,1)));
%col_vir=viridis(floor(.8*size(x_norm_cell,1)));
%colors=[col_vir(1:floor(size(x_norm_cell,1)/2),:);col_inf(1:(floor(size(x_norm_cell,1)/2)+1),:)];
figure(1);
clf;
set(gcf,'position',[31,1,1713,1010]);
ax=axes('position',[.05,.1,.93,.88]);
time_stamps=([1:size(x_cell{1},2)]/param.sig_sampling_rate)*down_sample_rate;
hold on
H=arrayfun(@(x) plot(time_stamps,x_norm_cell{x}+x,'color',colors_bip(x,:),'linewidth',1.5,'tag',...
    sprintf('ch %d, tag %s-%s',x,param.ch_name_for_valid_bipolar_diffs{x,1},param.ch_name_for_valid_bipolar_diffs{x,2})),[1:size(x_norm_cell,1)]);

if ~isempty(ops.timing)
    timings_ds=cell2mat(cellfun(@(x) floor(x./param.sig_sampling_rate),ops.timing,'uni',false));
     for kk=1:length(timings_ds)
        event=timings_ds(kk,:);
        corners_x=[min(event),max(event),max(event),min(event)];
        corners_y=[ax.YLim(1),ax.YLim(1),ax.YLim(2),ax.YLim(2)];
        pat=patch(corners_x,corners_y,[1,.4,.4]);
        pat.FaceAlpha=.1;
    end 
end 
set(gcf,'doublebuffer','on');
set(ax,'ytick',[1:size(x_norm_cell,1)]);
set(ax,'yticklabel','');

set(ax,'children',flipud(get(gca,'children')))
%set(ax,'ylim',[0,10]);
set(ax,'ylim',[0,size(x_norm_cell,1)+4]);
ax.XAxis.TickLength=[0.005,0.01];
ax.YAxis.TickLength=[0.005,0.01];
set(ax,'xlim',[0 t_len]);
pos=get(ax,'position');
Newpos=[pos(1) pos(2)-0.1 pos(3) 0.05];
xmax=max(time_stamps);
S=['set(gca,''xlim'',get(gcbo,''value'')+[0 ' num2str(t_len) '])'];

h=uicontrol('style','slider',...
    'units','normalized','position',Newpos,...
    'callback',S,'min',0,'max',xmax-t_len);

datacursormode on
dcm = datacursormode(gcf);
set(dcm,'UpdateFcn',@myupdatefcn)

waitfor(findobj('type','figure','number',1));
%% downsample the signal bipolar
fprintf(1, '>> Down-sampling signal \n');
signal_bipolar_hilbert_zs_decimated_loop = [];
signal_bioplar_hilbert_decimated_loop = [];
fprintf(1,'[');
for idx_channel=1:size(signal_bipolar_hilbert_zs_loop,2),
    fprintf('channel %d \n',idx_channel)
    %         % decimate band-power envelope to 12 Hz
    signal_bipolar_hilbert_zs_decimated_loop(:,idx_channel) = (decimate((signal_bipolar_hilbert_zs_loop(:,idx_channel)),decimation_factor));
    signal_bioplar_hilbert_decimated_loop(:,idx_channel) = (decimate((signal_bipolar_hilbert_loop(:,idx_channel)),decimation_factor));
       fprintf(1,'.');
end
fprintf(1,'] done\n')
param.preproc_step=[param.preproc_step;'downsampling_bipolar_gamma_envelope'];
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
dataout.ElectrodesInfo.IEDoutput= param.IED_results;

dataout=rmfield(dataout,'Data');

%dataout.filtered_signal=(signal');
dataout.MetaTags.DownSampleFrequency=param.downsamplingrate;
param_fields=fields(param);

for k=1:length(param_fields)
    dataout.MetaTags.(param_fields{k})=param.(param_fields{k});
end
%save(strrep(dataout.MetaTags.fullPath,'v3','v4'),'dataout','-v7.3');

% create a shorter copy
dataout=rmfield(dataout,'signal_hilbert');
dataout=rmfield(dataout,'signal_hilbert_zs');
dataout=rmfield(dataout,'signal_bipolar_hilbert_zs');
dataout=rmfield(dataout,'signal_bioplar_hilbert');

save(strrep(dataout.MetaTags.fullPath,'v3.mat','v4_compressed.mat'),'dataout','-v7.3');
end
