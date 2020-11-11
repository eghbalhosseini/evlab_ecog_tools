function [signal_broadband_out,...
    signal_bandpass_out,...
    signal_envelope_out,...
    signal_envelope_downsample_out,...
    signal_hilbert_decimated_loop_out,...
    signal_hilbert_zs_decimated_loop_out,...
    state_out,parameter_out,...
    ecog_param_out,...
    signal_hilbert_pca_decimated_loop_out,...
    signal_hilbert_pca_zs_decimated_loop_out,...
    subject_op_info]=filter_channels_using_schalk_v2(datafile,subject_op_info)
% varargin : 
%1 subject operation info

    if ~isempty(subject_op_info)
    %ecog_ch=subject_op_info.op_info.ecog_ch;
    ref_ch=subject_op_info.op_info.Ref;
    gnd_ch=subject_op_info.op_info.GND;
    bad_ch=subject_op_info.op_info.bad_channels;
    unselected_channels=subject_op_info.op_info.unselected_channels;
    transmit_chan=subject_op_info.op_info.transmit_chan;
    clean_channels=subject_op_info.op_info.clean_channels;
    sampling_rate=subject_op_info.op_info.sampling_rate;
    channel_noise_across_all_sess=subject_op_info.op_info.channel_noise_across_all_sess;
    channel_denoise_across_all_sess=subject_op_info.op_info.channel_denoise_across_all_sess;
    else 
    error('error: this function requires subj_op_info to run!');
    end 

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
ecog.param.samplingrate = 300; % decimation sampling rate 

audio.param.filter_type = 'IIR'; % IIR | FIR
% --- speech ---
audio.param.f_low_stop  = [  40];
audio.param.f_low_pass  = [  80];
audio.param.f_high_pass = [6000];
audio.param.f_high_stop = [8000];
audio.param.samplingrate = 300;

param=ecog.param;
param.channels=transmit_chan;
%% 

[ ~, ~, parameters ] = load_bcidat(datafile,[0 0]);

%%
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
        [bandpass{idx}.n,bandpass{idx}.Wn,bandpass{idx}.beta,bandpass{idx}.ftype] = kaiserord(bandpass{idx}.fcuts,bandpass{idx}.mags,bandpass{idx}.devs,parameters.SamplingRate.NumericValue);
        bandpass{idx}.n = bandpass{idx}.n + rem(bandpass{idx}.n,2);        
        % caclulate the filter coefficients in a,b design        
        bandpass{idx}.b = fir1(bandpass{idx}.n,bandpass{idx}.Wn,bandpass{idx}.ftype,kaiser(bandpass{idx}.n+1,bandpass{idx}.beta),'noscale');
        [bandpass{idx}.H,bandpass{idx}.f] = freqz(bandpass{idx}.b,1,1024,parameters.SamplingRate.NumericValue);

    end   
elseif strcmp(param.filter_type,'IIR'),
    % IIR bandpass definition using Butterworth filter
    for idx=1:length(param.f_low_stop),
        % define passband, stopband and ripples
        bandpass{idx}.Wp = [param.f_low_pass(idx) param.f_high_pass(idx)]/(parameters.SamplingRate.NumericValue/2); 
        bandpass{idx}.Ws = [param.f_low_stop(idx) param.f_high_stop(idx)]/(parameters.SamplingRate.NumericValue/2);
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
lowpass{1}.Wp = param.lowpass.Wp/(parameters.SamplingRate.NumericValue/2); 
lowpass{1}.Ws = param.lowpass.Ws/(parameters.SamplingRate.NumericValue/2);
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
signal = single([]);
states.StimType = [];
states.WordType = [];
states.IsRight  = [];
states.Response = [];
%
% go through all data files


    
    % load data file
    fprintf(1, '>> Loading data file %s \n',datafile);

    [ signal_loop, states_loop, parameters_loop ] = load_bcidat(datafile);
    signal_loop=signal_loop(:,param.channels);
    signal_broadband=signal_loop;
    
    % highpass filter
    fprintf(1, '>> Highpass filtering signal \n');

    fprintf(1,'[');
    parfor idx_channel=1:size(signal_loop,2) 
        warning('off', 'signal:filtfilt:ParseSOS');
        signal_loop(:,idx_channel) = single(filtfilt(highpass{1}.sos,highpass{1}.g,double(signal_loop(:,idx_channel)))); %#ok<PFBNS>
        fprintf(1,'.');    
    end
    fprintf(1,'] done\n');

    index_StimType = find(strcmp(parameters.Stimuli.RowLabels,'StimType'));
    
    index_WordType = find(strcmp(parameters.Stimuli.RowLabels,'Condition'));
    
    if isempty(index_WordType) 
        index_WordType = 7;
    end
    
    
    WordTypeList           = parameters_loop.Stimuli.NumericValue(index_WordType,:);   
    WordTypeList(isnan(WordTypeList)) = 0;
    StimTypeListText       = parameters_loop.Stimuli.Value(index_StimType,:); 
    StimTypeList           = zeros(size(WordTypeList));
    IsRightListText        = parameters_loop.Stimuli.Value(10,:); 
    IsRightList            = zeros(size(WordTypeList));
       
    
    for idx=1:length(StimTypeListText),
        if strcmp(StimTypeListText{idx},'word'),
            StimTypeList(idx) = 1;
        elseif strcmp(StimTypeListText{idx},'posttrial'),
            StimTypeList(idx) = 2;
        elseif strcmp(StimTypeListText{idx},'probe'),
            StimTypeList(idx) = 3;
        elseif strcmp(StimTypeListText{idx},'fixation'),
            StimTypeList(idx) = 4;
        else
            StimTypeList(idx) = -1; 
        end
    end
    
    IsRightLabels = unique(IsRightListText); 
    
    % new data format
    if strcmp(IsRightLabels{2},'0') && strcmp(IsRightLabels{3},'1')
        
        for idx=1:length(IsRightListText),
            if strcmp(IsRightListText{idx},'1'),
                IsRightList(idx) = 1;
            elseif strcmp(IsRightListText{idx},'0'),
                IsRightList(idx) = 2;
            else
                IsRightList(idx) = -1; 
            end
        end            
        
    else
    
        for idx=1:length(IsRightListText),
            if strcmp(IsRightListText{idx},'RIGHT'),
                IsRightList(idx) = 1;
            elseif strcmp(IsRightListText{idx},'WRONG'),
                IsRightList(idx) = 2;
            else
                IsRightList(idx) = -1; 
            end
        end    

    end
    
    index_valid            = find(states_loop.StimulusCode > 0);     
    index_invalid          = find(states_loop.StimulusCode == 0); 
    states_loop.WordType   = ones(size(states_loop.StimulusCode)) .* 0; 
    states_loop.StimType   = ones(size(states_loop.StimulusCode)) .* 0;    
    states_loop.IsRight    = ones(size(states_loop.StimulusCode)) .* 0;    
    states_loop.Response   = ones(size(states_loop.StimulusCode)) .* -1;  
    
    states_loop.WordType(index_valid) = WordTypeList(states_loop.StimulusCode(index_valid));    
    states_loop.StimType(index_valid) = StimTypeList(states_loop.StimulusCode(index_valid));    
    states_loop.IsRight(index_valid)  = IsRightList(states_loop.StimulusCode(index_valid));    
    
    index_isright_start = find(diff(double(states_loop.IsRight > 0)) == 1)+1;
    index_isright_stop  = find(diff(double(states_loop.IsRight > 0)) == -1);
    
    % YES = C = 67 | 99 
    %  NO = M = 77 | 109 
    
    buffer_before = parameters.SamplingRate.NumericValue * 1; % 1 sec 
    buffer_after  = parameters.SamplingRate.NumericValue * 2; % 2 sec
    
    states_loop.Response(index_invalid) = 0;
    
    for idx_trial = 1:length(index_isright_start),
        KeyDown = unique(states_loop.KeyDown((index_isright_start(idx_trial)-buffer_before):(index_isright_stop(idx_trial)+buffer_after)));
        KeyDown = intersect(KeyDown,[67,77,99,109]);
        
        if length(KeyDown) ~= 1                    % too many key's pressed or incorrect response
           TrialResponse = 0;
        elseif KeyDown == 67 || KeyDown == 99      % response is yes (1)
           TrialResponse = 1; 
        elseif KeyDown == 77 || KeyDown == 109     % response is no  (2) 
           TrialResponse = 2; 
        else                                       % incorrect response 
           TrialResponse = 0;
        end
        
        states_loop.Response((index_isright_start(idx_trial)):index_isright_stop(idx_trial)) = TrialResponse; 
        
    end
    
    
    % concatenate states and signals from data files
    signal              = [signal;signal_loop]; %#ok<AGROW>
    states=states_loop;
    states.StimType     = [states.StimType;states_loop.StimType];
    states.WordType     = [states.WordType;states_loop.WordType];    
    states.IsRight      = [states.IsRight;states_loop.IsRight];       
    states.Response     = [states.Response;states_loop.Response];
    clear signal_loop
    clear states_loop
    

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND DATA CHANNELS WITH SIGNIFICANT LINE-NOISE POWER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1, '> Searching for channels with significant 60 Hz noise \n');

fprintf(1,'[');
if isfield(parameters,'DeviceIDMaster')
    signal_noise=zeros(size(signal,2),1);
    % determine the number of g.USBamps where each amp has 16 channels
    num_amps = round(size(signal,2) / 16);

    % for each of these amps
    for idx_amp = 1:num_amps,
        
        % the the set of channels
        idx_low  = (idx_amp-1)*16+1;
        idx_high = (idx_amp-0)*16+0;

        % get the channels that actually had signals 
        list_channels = intersect(idx_low:idx_high,transmit_chan);
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
param.channels_noise = find(signal_noise > (mean(signal_noise(transmit_chan))+1.5*std(signal_noise(transmit_chan))));

%param.channels_noise = intersect(param.channels_noise,param.channels);
%param.channels_deselect = sort(unique([param.channels_noise param.channel_ground param.channel_ref param.channel_seizure]));
param.channels_deselect = sort(unique([reshape(param.channels_noise,1,[]),unique([ref_ch,gnd_ch]) ]));
if isequal(param.channels_deselect,unselected_channels)
    fprintf('the noise is the same for session and across sessions')
else
    fprintf('the noise is NOT the same for session and across sessions, selecting channels from spec files')
    param.channels_deselect=unselected_channels;
end 
param.channels_selected = setdiff(param.channels,param.channels_deselect);

fprintf(1, '> Not using %d channels: ',length(param.channels_deselect));
fprintf(1, '%d ',param.channels_deselect);
fprintf(1, '\n');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REMOVE COMMON NOISE USING A COMMON-AVERAGE REFERENCE FILTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
fprintf(1, '> Common average filtering signal \n');

% now calculate the common-average reference signal again, without those channels that have significant line-noise
fprintf(1,'[');
if isfield(parameters,'DeviceIDMaster') && param.filter_car > 0,
    
    % determine the number of g.USBamps where each amp has 16 channels
    num_amps = round(size(signal,2) / 16);

    % for each of these amps
    for idx_amp = 1:num_amps,
        idx_low  = (idx_amp-1)*16+1;
        idx_high = (idx_amp-0)*16+0;

        % exclude the channels that had signifiant line-noise
        list_channels = intersect(idx_low:idx_high,param.channels_selected);
        % check if any channels are left and 
        if ~isempty(list_channels) && length(list_channels)>1 ,   
            
            % calculate the common average reference signal 
            signal_mean = mean(signal(:,intersect(idx_low:idx_high,param.channels_selected)),2);
        
            % subtract the common average signal from each channel of this amp
            for idx_ch=idx_low:idx_high,
                signal(:,idx_ch) = signal(:,idx_ch) - signal_mean;
                fprintf(1,'.');
            end
        else
            % if no channel is left then the signal is not filtered
            for idx_ch=idx_low:idx_high,
                signal(:,idx_ch) = signal(:,idx_ch);
                fprintf(1,'.');
            end
        end
    end  
end 
fprintf(1,'] done\n');



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEASSURE LINE-NOISE POWER AFTER SIGNAL PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1, '> Meassuring 60 Hz noise power after signal processing \n');

fprintf(1,'[');
% for each channels calculate the root-mean-square line-noise power
parfor idx_channel=1:size(signal,2),
    signal_noise_after(idx_channel) = mean(sqrt(filtfilt(peak.b,peak.a,double(signal(:,idx_channel))).^2),1); %#ok<PFBNS>
    fprintf(1,'.');    
end
fprintf(1,'] done\n');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REPORT LINE-NOISE POWER AFTER SIGNAL PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1, '> Reduced 60 Hz noise from %.2f to %.2f uV',mean(signal_noise_before(param.channels_selected)),mean(signal_noise_after(param.channels_selected)));
fprintf(1, '\n');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILTER OUT REMAINING LINE-NOISE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%% plot results 
%% plot post filtering 
t_length=10e4;
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
for kk=1:kk_total
    clf;
    set(gcf,'color',[.7,.7,.7],'position',[17,1,2454,1337])
    ax=axes('position',[.02,.02,.95,.95]);
    hold on
    session_end=min([(kk*t_length),size(signal,1)]);
    t_window=((kk-1)*t_length+1):2:(session_end);
    
    arrayfun(@(x) plot(t_window,x_norm_cell{x}(t_window)+x-1,'color',colors(x,:)),param.channels)
    
    set(ax,'ytick',[1:size(x_norm_cell,1)]);
    set(ax,'yticklabel','')
    arrayfun(@(x) text(t_window(1),x,num2str(x),'Color',colors(x,:),'HorizontalAlignment','right','VerticalAlignment','middle'),param.channels);
    set(ax,'ylim',[0,size(x_norm_cell,1)+1])
    shg;
    title({'Session : post common source averaging',[num2str(kk),' of ', num2str(kk_total)]})
    pause
end 
close all; 

%%
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
for idx_bandpass = 1:length(bandpass)
    
fprintf(1, '> Extracting Frequency band %d/%d (%d-%d Hz)\n',idx_bandpass,length(bandpass),param.f_low_pass(idx_bandpass),param.f_high_pass(idx_bandpass));       
    %% band-pass filtering signal
    fprintf(1, '>> Band-pass filtering signal %d-%d Hz\n',param.f_low_pass(idx_bandpass),param.f_high_pass(idx_bandpass));

    warning('off', 'signal:filtfilt:ParseSOS');

    fprintf(1,'[');
    parfor idx_channel=1:size(signal,2),

        % check if a,b or zero-pole-gain filter design
        if isfield(bandpass{idx_bandpass},'b'),  %#ok<PFBNS>
            % filter signal using a,b coeffiienct 
            signal_loop(:,idx_channel) = single(filtfilt(bandpass{idx_bandpass}.b,1,double(signal(:,idx_channel))));
        else
            % filter signal using zero-pole-gain coefficients
            signal_loop(:,idx_channel) = single(filtfilt(bandpass{idx_bandpass}.sos,bandpass{idx_bandpass}.g,double(signal(:,idx_channel))));
        end
        fprintf(1,'.');
    end

    %% remove IIR artifact
    if ~isfield(bandpass{idx_bandpass},'b')
        % set first and last few samples of the band-pass filtered signal to zero
        signal_loop(1:bandpass{idx_bandpass}.n*30      ,:) = 0;
        signal_loop(end-bandpass{idx_bandpass}.n*30:end,:) = 0;
    end

    fprintf(1,'] done\n');
    signal_bandpass{idx_bandpass,1}=signal_loop;
    %% root-mean-square filter signal
%     fprintf(1, '>> Root-mean-square filtering signal \n');
% 
%     fprintf(1,'[');
%     parfor idx_channel=1:size(signal_loop,2),
%         % calculate band-power using bay applying the root-mean-square to the band-pass filtered signal
%         signal_loop(:,idx_channel) = sqrt(signal_loop(:,idx_channel).^2);
%         fprintf(1,'.');
%     end
%     fprintf(1,'] done\n');

    %% hilbert filter signal
    fprintf(1, '>> Extracting signal envelope \n');
    
    fprintf(1,'['); 
    parfor idx_channel=1:size(signal_loop,2),
        signal_loop(:,idx_channel) = abs(hilbert(double(signal_loop(:,idx_channel))));
        fprintf(1,'.');        
    end
    fprintf(1,'] done\n');
    
    
    %% low-pass filter signal 
    fprintf(1, '>> Low-pass filtering signal \n');

    warning('off', 'signal:filtfilt:ParseSOS');

    fprintf(1,'[');
    parfor idx_channel=1:size(signal_loop,2),
        % extract the envelope of the band-power by low-pass filtering the root-mean square signal 
        signal_loop(:,idx_channel) = single(filtfilt(lowpass{1}.sos,lowpass{1}.g,double(signal_loop(:,idx_channel)))); %#ok<PFBNS>
        fprintf(1,'.');
    end
    fprintf(1,'] done\n');
    signal_envelope{idx_bandpass,1} = signal_loop;
     %% down-sample signal
     fprintf(1, '>> Down-sampling signal \n');
%        
     % calculate decimation factor for 12 Hz sampling rate
     decimation_factor = parameters.SamplingRate.NumericValue / param.samplingrate;
     parameters.decimation_factor=decimation_factor;
%     % decimate StimulusCode to 12 Hz
     StimType = downsample(double(states.StimType),decimation_factor);
     WordType = downsample(double(states.WordType),decimation_factor);    
     IsRight  = downsample(double(states.IsRight),decimation_factor);    
     Response = downsample(double(states.Response),decimation_factor); 
     states.StimulusCodeDownsample=downsample(double(states.StimulusCode),decimation_factor);
%     % initialize variable for decimating band-power envelope to 12 Hz
     signal_decimated_loop = zeros(size(decimate(double(signal_loop(:,1)),decimation_factor),1),size(signal_loop,2),'single');
%       
     fprintf(1,'[');
     parfor idx_channel=1:size(signal,2),
% 
%         % decimate band-power envelope to 12 Hz
         signal_decimated_loop(:,idx_channel) = single(decimate(double(signal_loop(:,idx_channel)),decimation_factor));
         fprintf(1,'.');
     end
%     
%     
     % return variable for each band-pass
     signal_envelope_downsample{idx_bandpass,1} = signal_decimated_loop;
     
     clear signal_decimated_loop;
     clear signal_loop;
     fprintf(1,'] done\n');

end 
%%    
    %% compute high gamma based on chang-lab method, see : 
    % Dichter, Benjamin K., Jonathan D. Breshears, Matthew K. Leonard, and Edward F. Chang. 2018.
    % ?The Control of Vocal Pitch in Human Laryngeal Motor Cortex.? Cell 174 (1): 21?31.e9.
    fprintf(1, '>> Extracting high gamma envelope based on chang-lab method  \n');
    
    fprintf(1,'['); 
    parfor idx_channel=1:size(signal,2)
        filter_bank={};
        for s=1:length(cfs)
        filter_bank{s}=gaussian_filter(transpose(signal(:,idx_channel)),parameters.SamplingRate.NumericValue,cfs(s),sds(s));
        end 
        signal_hilbert_bp = cellfun(@abs,hilbert_transform(double(transpose(signal(:,idx_channel))),...
                                                  parameters.SamplingRate.NumericValue,...
                                                  filter_bank),'UniformOutput',false);
        % mean across 
        signal_hilbert_bp_cell=signal_hilbert_bp;
        signal_hilbert_bp=mean(cell2mat(signal_hilbert_bp),1);
        % to do : do pca instead 
        % De-mean
        signal_hilbert_bp_mat=transpose(cell2mat(signal_hilbert_bp_cell));
        X = bsxfun(@minus,signal_hilbert_bp_mat,mean(signal_hilbert_bp_mat));
        [coeff,score,latent,tsquared,explained,mu] = pca(X);
        signal_hilbert_pca=transpose(score*coeff(:,1));

        % zscore the signal with respect to entire experimental block 
        % to do 
        means = mean(signal_hilbert_bp);
        stds = std(signal_hilbert_bp);
        signal_hilbert_zs = (signal_hilbert_bp - means) / stds;
        
        means = mean(signal_hilbert_pca);
        stds = std(signal_hilbert_pca);
        signal_hilbert_pca_zs = (signal_hilbert_pca - means) / stds;
        %
        signal_hilbert_loop(:,idx_channel)=signal_hilbert_bp;
        signal_hilbert_zs_loop(:,idx_channel)=signal_hilbert_zs;
        % 
        signal_hilbert_pca_loop(:,idx_channel)=signal_hilbert_pca;
        signal_hilbert_pca_zs_loop(:,idx_channel)=signal_hilbert_pca_zs;
        fprintf(1,'.');        
    end
    fprintf(1,'] done\n');
    
     %% down-sample signal
     fprintf(1, '>> Down-sampling signal \n');
%        
     % calculate decimation factor for 12 Hz sampling rate
     decimation_factor = parameters.SamplingRate.NumericValue / param.samplingrate;
     parameters.decimation_factor=decimation_factor;

     signal_hilbert_zs_decimated_loop = zeros(size(decimate(double(signal_hilbert_zs_loop(:,1)),decimation_factor),1),size(signal_hilbert_zs_loop,2),'single');
     signal_hilbert_decimated_loop = zeros(size(decimate(double(signal_hilbert_zs_loop(:,1)),decimation_factor),1),size(signal_hilbert_zs_loop,2),'single');
%    % 
     signal_hilbert_pca_decimated_loop = zeros(size(decimate(double(signal_hilbert_zs_loop(:,1)),decimation_factor),1),size(signal_hilbert_zs_loop,2),'single');
     signal_hilbert_pca_zs_decimated_loop = zeros(size(decimate(double(signal_hilbert_zs_loop(:,1)),decimation_factor),1),size(signal_hilbert_zs_loop,2),'single');
%    
     fprintf(1,'[');
     for idx_channel=1:size(signal_hilbert_zs_loop,2),
% 
%         % decimate band-power envelope to 12 Hz
         signal_hilbert_zs_decimated_loop(:,idx_channel) = single(decimate(double(signal_hilbert_zs_loop(:,idx_channel)),decimation_factor));
         signal_hilbert_decimated_loop(:,idx_channel) = single(decimate(double(signal_hilbert_loop(:,idx_channel)),decimation_factor));
         
         %signal_hilbert_pca_zs_decimated_loop(:,idx_channel) = single(decimate(double(signal_hilbert_pca_zs_loop(:,idx_channel)),decimation_factor));
         %signal_hilbert_pca_decimated_loop(:,idx_channel) = single(decimate(double(signal_hilbert_pca_loop(:,idx_channel)),decimation_factor));
      
         
         fprintf(1,'.');
     end
     fprintf(1,'] done\n');


%% 

parameter_out=parameters; 
signal_broadband_out=signal_broadband;
state_out=states;
signal_bandpass_out=signal_bandpass;
signal_envelope_out=signal_envelope;
signal_envelope_downsample_out=signal_envelope_downsample;
signal_hilbert_zs_decimated_loop_out=signal_hilbert_zs_decimated_loop;
signal_hilbert_decimated_loop_out=signal_hilbert_decimated_loop;
ecog_param_out=param;
signal_hilbert_pca_zs_decimated_loop_out=signal_hilbert_pca_zs_decimated_loop;
signal_hilbert_pca_decimated_loop_out=signal_hilbert_pca_decimated_loop;
% 
%subject_op_info.op_info.clean_channels=param.channels_selected;
subject_op_info.op_info.(strcat('channel_noise_sess_',datafile(end-6:end-4)))=signal_noise_before;
subject_op_info.op_info.(strcat('channel_denoise_sess_',datafile(end-6:end-4)))=signal_noise_after;

end 