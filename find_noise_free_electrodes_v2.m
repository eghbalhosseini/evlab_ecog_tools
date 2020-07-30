function ops_out=find_noise_free_electrodes_v2(varargin)
p=inputParser();
addParameter(p, 'datafile', 'test');
addParameter(p, 'op_info', struct);
addParameter(p, 'exp_name', 'test');
addParameter(p, 'save_plots', false);
addParameter(p, 'plot_save_path', '~/ECOG/crunched/plots/');
parse(p, varargin{:});
ops = p.Results; 
save('ops_here','ops')

if ~isempty(ops.op_info)
    f=fields(ops.op_info)';
    for i =f ,eval(sprintf('%s=ops.op_info.%s;',i{1},i{1})); end 
else
    error('error: this function requires subj_op_info to run!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  CHECK DATA FILE DEPENDENCIES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters={};
fprintf(1, '> Checking data file dependencies \n');
for idx=1:length(ops.datafile)
    % load just the header of this data file
    [ ~, ~, parameter ] = load_bcidat(ops.datafile{idx},[0 0]); %#ok<*SAGROW>.
    parameters{idx}=parameter;
end
% compare sampling rates
s_rates=cellfun(@(x) x.SamplingRate.NumericValue,parameters);
assert(all(s_rates==s_rates(1)),'error: Data files dont have the same SamplingRate!');
sampling_rate=unique(s_rates);
% compare source channels
src_ch=cellfun(@(x) x.SourceCh.NumericValue,parameters);
assert(all(src_ch==src_ch(1)),'error: Data files dont have the same SourceNumber!');
if ~isempty(ecog_channels)
    param.channels=ecog_channels;
else 
    param.channels=1:unique(src_ch);
end 

%%
param.filter_car  = 1;     % 0 == off % 1 == on
param.filter_type = 'IIR'; % IIR | FIR                             
% --- highpass filter --- 
param.highpass.Wp = 0.50; % Hz
param.highpass.Ws = 0.05; % Hz
param.highpass.Rp = 3;    % dB
param.highpass.Rs = 30;   % dB
% --- lowpass filter --- 
param.lowpass.Wp = 100;     % Hz
param.lowpass.Ws = 110;     % Hz
param.lowpass.Rp = 3;     % dB
param.lowpass.Rs = 30;    % dB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE HIGH-PASS FILTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define passband, stopband and attenuation
highpass{1}.Wp = param.highpass.Wp/(sampling_rate/2); 
highpass{1}.Ws = param.highpass.Ws/(sampling_rate/2);
highpass{1}.Rp = param.highpass.Rp; 
highpass{1}.Rs = param.highpass.Rs;
% calculate the minimum filter order for butterworth filter
[highpass{1}.n,highpass{1}.Wn] = buttord(highpass{1}.Wp,highpass{1}.Ws,highpass{1}.Rp,highpass{1}.Rs);
highpass{1}.n = highpass{1}.n + rem(highpass{1}.n,2);
% caclulate the filter coefficients in Zero-Pole-Gain design
[highpass{1}.z,highpass{1}.p,highpass{1}.k] = butter(highpass{1}.n,highpass{1}.Wn,'high');
[highpass{1}.sos,highpass{1}.g]=zp2sos(highpass{1}.z,highpass{1}.p,highpass{1}.k);
highpass{1}.h=dfilt.df2sos(highpass{1}.sos,highpass{1}.g);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE LINE-NOISE FILTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
peak.fcenter = 60;
peak.bw      = 0.001;

% calculate the IIR-peak filter coefficients in a,b format 
peak.wo = peak.fcenter/(sampling_rate/2);  
peak.bw = peak.bw;
[peak.b,peak.a] = iirpeak(peak.wo,peak.bw);  

% define the harmonics of line noise frequency
param.filter.notch.fcenter = [60,120,180,240];
param.filter.notch.bw      = ones(1,length(param.filter.notch.fcenter)).*0.001;

% calculate the IIR-peak filter coefficients in a,b format 
for idx = 1:length(param.filter.notch.fcenter),
    notch{idx}.wo = param.filter.notch.fcenter(idx)/(sampling_rate/2);  
    notch{idx}.bw = param.filter.notch.bw(idx);
    [notch{idx}.b,notch{idx}.a] = iirnotch(notch{idx}.wo,notch{idx}.bw);  
end
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD AND CONCATENATE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signal = [];
session_length=[];
buffer_before = sampling_rate * 1; % 1 sec 
buffer_after  = sampling_rate * 2; % 2 sec
    
% go through all data files
for idx_file = 1:length(ops.datafile)
    datafile=ops.datafile{idx_file};
    fprintf(1, '> Processing data file %d/%d \n', idx_file,length(ops.datafile));
    % load data file
    [ signal_loop, ~, ~ ] = load_bcidat(datafile);
    signal_loop=double(signal_loop(:,param.channels));
    session_length=[session_length;size(signal_loop,1)];
    % highpass filter
    fprintf(1, '>> Highpass filtering signal \n');
    fprintf(1,'[');
    parfor idx_channel=1:size(signal_loop,2) 
        warning('off', 'signal:filtfilt:ParseSOS');
        signal_loop(:,idx_channel) = filtfilt(highpass{1}.sos,highpass{1}.g,(signal_loop(:,idx_channel))); %#ok<PFBNS>
        fprintf(1,'.');    
    end
    fprintf(1,'] done\n');
    %
    % concatenate states and signals from data files
    signal              = [signal;signal_loop]; %#ok<AGROW>
    clear signal_loop
end

%% plot the dataset throughout 
session_length=[1;session_length];
session_start=cumsum(session_length);
param.channels_deselect = []; %no channels deselected yet!
fprintf('Visually inspect the pre- line noise removal and common source averaging signal \n')
print_figure = false;
figure_label = "";
%plot_signal_over_time(signal,param,session_start,'pre-common source averaging', ops.save_plots,print_figure, ops.plot_save_path, figure_label,ops.op_info.sub_id)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEASURE LINE-NOISE POWER BEFORE SIGNAL PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, '> Measuring 60 Hz noise power before signal processing \n');
fprintf(1,'[');
parfor idx_channel=1:size(signal,2),
    % calculate average root-mean-square of the line-noise
    signal_noise_before(idx_channel) = mean(sqrt(filtfilt(peak.b,peak.a,double(signal(:,idx_channel))).^2),1); %#ok<PFBNS>
    fprintf(1,'.');    
end
fprintf(1,'] done\n');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND DATA CHANNELS WITH SIGNIFICANT LINE-NOISE POWER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1, '> Searching for channels with significant 60 Hz noise \n');

fprintf(1,'[');
if isfield(parameters{1},'DeviceIDMaster')
    signal_noise=zeros(size(signal,2),1);
    % determine the number of g.USBamps where each amp has 16 channels
    num_amps = round(size(signal,2) / 16);

    % for each of these amps
    for idx_amp = 1:num_amps,
        % the the set of channels
        idx_low  = (idx_amp-1)*16+1;
        idx_high = (idx_amp-0)*16+0;

        % get the channels that actually had signals 
        list_channels = intersect(idx_low:idx_high,param.channels);
        % check if any channels are left and 
        if ~isempty(list_channels),
            % calculate the common average reference signal 
            signal_mean = mean(signal(:,list_channels),2);
        else
            % if no channels is left then the common average reference signal is zero
            signal_mean = zeros(size(signal,1),1);
        end
            
        % for each channel on this amp
        for idx_ch=1:length(list_channels),
            % subtract the common average signal from each channel of this amp          
            signal_preliminary = double(signal(:,list_channels(idx_ch)));% - double(signal_mean);
            
            % calculate the residual line-noise root-mean-square power
            v=mean(sqrt(filtfilt(peak.b,peak.a,signal_preliminary).^2),1);
            signal_noise(list_channels(idx_ch)) = v;
            fprintf(1,'.');
            
        end
    end  
else
    warning('error: common average reference filtering is only supported for g.USBamp data');  %#ok<WNTAG>
end

clear signal_mean
fprintf(1,'] done\n');
%%
% find those channels for which the line-noise power is 1.5 standard deviations higher than the average 

param.channels_noise = param.channels( signal_noise > (mean(signal_noise)+1.5*std(signal_noise)));
%param.channels_deselect = sort(unique([param.channels_noise param.channel_ground param.channel_ref param.channel_seizure]));

param.channels_deselect = sort(unique([reshape(param.channels_noise,1,[]),unique([REF,GND,bad_channels])]));
param.channels_selected = setdiff(param.channels,param.channels_deselect);

fprintf(1, '> Not using %d channels for common source averaging: ',length(param.channels_deselect));
fprintf(1, '%d ',param.channels_deselect);
fprintf(1, '\n');

%% 
fprintf(1, '> Common average filtering signal \n');

% now calculate the common-average reference signal again, without those channels that have significant line-noise
fprintf(1,'[');
if isfield(parameters{1},'DeviceIDMaster') && param.filter_car > 0,
    % determine the number of g.USBamps where each amp has 16 channels
    num_amps = round(size(signal,2) / 16);
    % for each of these amps
    for idx_amp = 1:num_amps,
        idx_low  = (idx_amp-1)*16+1;
        idx_high = (idx_amp-0)*16+0;

        % exclude the channels that had signifiant line-noise
        list_channels = intersect(idx_low:idx_high,param.channels_selected);
        % check if any channels are left and 
        if ~isempty(list_channels) && length(list_channels)>1 
            % calculate the common average reference signal 
            signal_mean = mean(signal(:,list_channels),2);
            % subtract the common average signal from each channel of this amp
            for idx_ch=list_channels
                signal(:,idx_ch) = signal(:,idx_ch) - signal_mean;
                fprintf(1,'.');
            end
        elseif ~isempty(list_channels) && length(list_channels)==1 ,
            param.channels_selected(param.channels_selected==list_channels)=[];
            param.channels_deselect=[param.channels_deselect,list_channels];
            fprintf('remove channel with no additional pair for common reference averaging')
        else
            % if no channel is left then the signal is not filtered
            for idx_ch=list_channels
                signal(:,idx_ch) = signal(:,idx_ch);
                fprintf(1,'.');
            end
        end
    end  
end 
fprintf(1,'] done\n');
%% 
fprintf(1, '> Measuring 60 Hz noise power after signal processing \n');

fprintf(1,'[');
% for each channels calculate the root-mean-square line-noise power
parfor idx_channel=1:size(signal,2)
    signal_noise_after(idx_channel) = mean(sqrt(filtfilt(peak.b,peak.a,(signal(:,idx_channel))).^2),1); %#ok<PFBNS>
    fprintf(1,'.');    
end
fprintf(1,'] done\n');
%% 
fprintf(1, '> Reduced 60 Hz noise from %.2f to %.2f uV',mean(signal_noise_before(param.channels_selected)),mean(signal_noise_after(param.channels_selected)));
fprintf(1, '\n');


%% 
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
    signal(:,idx_channel) = single(signal_preliminary);
    
    fprintf(1,'.');
end
fprintf(1,'] done\n');
%% plot post filtering 
fprintf('Visually inspect the post filtering signal for manual noisy channel removal \n');
print_figure = true;
figure_label = "PRE-REMOVAL_";
plot_signal_over_time(signal,param,session_start,'post 60 Hz noise removal, common source averaging, and notch filtering',ops.save_plots,print_figure, ops.plot_save_path,figure_label,ops.op_info.sub_id)
print_figure = false;
figure_label = "";
%% 
x = 1;
channels_removed_automatically = param.channels_deselect;
channels_removed_manually = [];
num_removals = 0;
while ~isempty(x)
    prompt='Additional channels to remove? Format: [1,2] or [1:10, 20:25, 31] (if you are done removing channels, press Enter)  \n';
    x=input(prompt);
    if ~isempty(x)
         %save the channels that were automatically removed due to noise
        fprintf(['Removing channels... \n\n']);
        num_removals = num_removals + 1;
        channels_removed_manually = [channels_removed_manually x];
        param.channels_deselect=sort([param.channels_deselect,x]); %all channels that are being removed, including the user specified channels and channels removed due to noise
        param.channels_selected = setdiff(param.channels,param.channels_deselect);
        fprintf('visually inspect signal for manual channel removal \n');
        print_figure = true;
        figure_label = strcat("POST-REMOVAL", num2str(num_removals), "_");
        plot_signal_over_time(signal,param,session_start,'post 60 Hz noise removal, common source averaging, and notch filtering and manual removal of channels',ops.save_plots,print_figure,ops.plot_save_path,figure_label,ops.op_info.sub_id);
        print_figure = false;
        figure_label = "";
    end
end

prompt='enter your full name please :) -- ';
user_name = input(prompt, 's');
param.valid_channels=zeros(size(signal,2),1);
param.valid_channels(param.channels_selected)=1;
%% add the results to subject_op_info, experiment specific!!
ops.op_info.visually_inspected_by = user_name;
ops.op_info.visually_inspected = 1;
ops.op_info.visually_inspected_date = datestr(now, 'yyyy/mm/dd-HH:MM');
ops.op_info.experiment_name = ops.exp_name;
ops.op_info.channels_removed_sig_line_noise = param.channels_noise;
ops.op_info.channels_removed_automatically_all = channels_removed_automatically;
ops.op_info.channels_removed_manually= unique(channels_removed_manually);
ops.op_info.transmit_chan=param.channels;
ops.op_info.unselected_channels=param.channels_deselect;
ops.op_info.clean_channels=param.channels_selected;
ops.op_info.valid_channels=param.channels_selected;
ops.op_info.sampling_rate=sampling_rate;
ops.op_info.channel_noise_across_all_sess=signal_noise_before;
ops.op_info.channel_denoise_across_all_sess=signal_noise_after;
ops.op_info.filter_param=param;
ops_out=ops;
fprintf('done \n');
end 

function []= plot_signal_over_time(signal,param,session_start,plot_title,save_plots,print_figure,save_path, figure_label, sub_id)
    filename = "";
    fprintf('Figure is loading... \n');
    fprintf('Figure will advance through time when you press any key.\n')
    t_length=25e4;
    kk_total=ceil(size(signal,1)/t_length);
    x_cell=mat2cell(signal',ones(1,size(signal,2)));
    valid_channels=ones(1,size(signal,2));
    valid_channels(param.channels_deselect)=nan;
    x_cell=arrayfun(@(x) x_cell{x}*valid_channels(x),1:size(signal,2),'UniformOutput',false);
    x_cell=x_cell';
    min_max=nanmean(cell2mat(cellfun(@(y) [min(y),max(y)],x_cell,'UniformOutput',false)));
    x_norm_cell=cellfun(@(x) (x-min_max(1))./(.5*(min_max(2)-min_max(1))),x_cell,'UniformOutput',false);
    %
    figure(1);
   
    
    col_inf=inferno(floor(.8*size(x_norm_cell,1)));
    col_vir=viridis(floor(.8*size(x_norm_cell,1)));
    colors=[col_vir(1:floor(size(x_norm_cell,1)/2),:);col_inf(1:(floor(size(x_norm_cell,1)/2)+1),:)];

    %all_filenames = {};
    for kk=1:kk_total
        clf;
        set(0,'units','pixels');
        screen = get(0,'ScreenSize');
        set(gcf,'color',[.7,.7,.7],'position', screen)
        ax=axes('position',[.02,.02,.95,.95]);
        hold on
        session_end=min([(kk*t_length),size(signal,1)]);
        t_window=((kk-1)*t_length+1):2:(session_end);
        arrayfun(@(x) plot(t_window,x_norm_cell{x}(t_window)+x-1,'color',colors(x,:)),param.channels)
        
        set(ax,'ytick',[1:size(x_norm_cell,1)]);
        set(ax,'yticklabel','')
        arrayfun(@(x) text(t_window(1),x,num2str(x),'Color',colors(x,:),'HorizontalAlignment','right','VerticalAlignment','middle'),param.channels);
        set(ax,'ylim',[0,size(x_norm_cell,1)+1]);
        if any(ismember(t_window,session_start))
        arrayfun(@(x) plot([t_window(x),t_window(x)],[0,size(x_norm_cell,1)+1],'k-','LineWidth',2),find(ismember(t_window,session_start)))
        session_name=num2str(find(ismember(session_start,t_window)));
        arrayfun(@(x) text(t_window(x),size(x_norm_cell,1)+1,['sess: ',session_name],'VerticalAlignment','bottom','fontsize',20),find(ismember(t_window,session_start)))
        end 
        shg;
        title({plot_title,[num2str(kk),' of ', num2str(kk_total)]})
        
        
       
        
        if(print_figure && save_plots)
            if ~exist(strcat(save_path,'channel_plots'))
                mkdir(strcat(save_path,'channel_plots'))
            end
            %filename = strcat(save_path, 'temp', num2str(kk), '.pdf');
            filename = strcat(save_path, 'channel_plots', filesep, sub_id, '_', figure_label, num2str(kk),'.pdf')
            %all_filenames{end+1} = strcat(save_path, filename);
            
            %f=gcf; ?
            %f.Units='Inches';?
            %f.PaperOrientation='portrait';?
            %pos = get(gcf,'Position');?
%             %letter_paper_ratio=11/8.5;?
%             set(gcf,'PaperPositionMode','auto');
%             set(gcf,'PaperUnits','inches');
%             set(gcf,'PaperSize',[8.5,11])?;
%             pos=gcf.Position;?
%             set(gcf,'position',[ pos(1),pos(2),8.5,11])
            set(gcf,'PaperOrientation', 'landscape');
            print(gcf,'-painters','-fillpage', '-dpdf', filename);
            %export_fig(filename, '-append');
          
           % print(gcf,'-painters', '-depsc', strcat(analysis_path,info.subject,'/',info.subject,'_','fig_for_U01','.pdf'));

        end      
        pause
    end

    % set(gcf,'PaperPositionMode','auto'); %used to modify how we save the figure
    %if(print_figure && save_plots)
        %concatenate all pdfs into single file
    %    output_filename = strcat(save_path, sub_id, '_', figure_label,'.pdf')
    %    output_filename
    %    all_filenames
    %    append_pdfs(output_filename, all_filenames{:});
        %print(output_filename)
    %end
    

    
    close all;
end

