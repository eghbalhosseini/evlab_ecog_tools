function ops_out=filter_channels_using_gaussian(varargin)
% varargin : 
p=inputParser();
addParameter(p, 'datafile', 'test');
addParameter(p, 'op_info', struct);
parse(p, varargin{:});
ops = p.Results;
if ~isempty(ops.op_info)
    f=fields(ops.op_info)';
    for i =f ,eval(sprintf('%s=ops.op_info.%s;',i{1},i{1})); end 
else
    error('error: this function requires subj_op_info to run!');
end
%% 
try param=ops.op_info.filter_param;
catch err
    ecog.param.filter_car  = 1;     % 0 == off % 1 == on
    ecog.param.filter_type = 'IIR'; % IIR | FIR
    
    
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
    
    param=ecog.param;
    peak.fcenter = 60;
    peak.bw      = 0.001;
    
    % calculate the IIR-peak filter coefficients in a,b format
    peak.wo = peak.fcenter/(sampling_rate/2);
    peak.bw = peak.bw;
    [peak.b,peak.a] = iirpeak(peak.wo,peak.bw);
    
    % define the harmonics of line noise frequency
    param.filter.notch.fcenter = [60,120,180,240];
    param.filter.notch.bw      = ones(1,length(param.filter.notch.fcenter)).*0.001;
end
ops.downsamplingrate = 300; % decimation sampling rate 
% 
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
% calculate the IIR-peak filter coefficients in a,b format 
for idx = 1:length(param.filter.notch.fcenter),
    notch{idx}.wo = param.filter.notch.fcenter(idx)/(sampling_rate/2);  
    notch{idx}.bw = param.filter.notch.bw(idx);
    [notch{idx}.b,notch{idx}.a] = iirnotch(notch{idx}.wo,notch{idx}.bw);  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD AND CONCATENATE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signal = double([]);
states.StimType = [];
states.WordType = [];
states.IsRight  = [];
states.Response = [];
%
% go through all data files
% load data file
fprintf(1, '>> Loading data file %s \n',ops.datafile);
[ signal_loop, states, parameters ] = load_bcidat(ops.datafile);
%signal=signal_loop(:,clean_channels);
signal=signal_loop;
% highpass filter
fprintf(1, '>> Highpass filtering signal \n');
fprintf(1,'[');
parfor idx_channel=1:size(signal,2)
    warning('off', 'signal:filtfilt:ParseSOS');
    signal(:,idx_channel) = filtfilt(highpass{1}.sos,highpass{1}.g,double(signal(:,idx_channel)));
    fprintf(1,'.');
end
fprintf(1,'] done\n');

index_StimType = find(strcmp(parameters.Stimuli.RowLabels,'StimType'));
index_WordType = find(strcmp(parameters.Stimuli.RowLabels,'Condition'));
index_IsRight = find(strcmp(parameters.Stimuli.RowLabels,'IsRight') |...
strcmp(parameters.Stimuli.RowLabels,'IsProbeCorrect')) ;
if isempty(index_WordType)
    index_WordType = 7;
end


WordTypeList           = parameters.Stimuli.NumericValue(index_WordType,:);
WordTypeList(isnan(WordTypeList)) = 0;
StimTypeListText       = parameters.Stimuli.Value(index_StimType,:);
StimTypeList           = zeros(size(WordTypeList));
IsRightListText        = parameters.Stimuli.Value(index_IsRight,:);
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

if ~isempty(IsRightLabels) %only if these exist -HS
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
end
    
index_valid            = find(states.StimulusCode > 0);
index_invalid          = find(states.StimulusCode == 0);
states.WordType   = ones(size(states.StimulusCode)) .* 0;
states.StimType   = ones(size(states.StimulusCode)) .* 0;
states.IsRight    = ones(size(states.StimulusCode)) .* 0;
states.Response   = ones(size(states.StimulusCode)) .* -1;
    
states.WordType(index_valid) = WordTypeList(states.StimulusCode(index_valid));
states.StimType(index_valid) = StimTypeList(states.StimulusCode(index_valid));
states.IsRight(index_valid)  = IsRightList(states.StimulusCode(index_valid));

index_isright_start = find(diff(double(states.IsRight > 0)) == 1)+1;
index_isright_stop  = find(diff(double(states.IsRight > 0)) == -1);

% YES = C = 67 | 99
%  NO = M = 77 | 109

buffer_before = parameters.SamplingRate.NumericValue * 1; % 1 sec
buffer_after  = parameters.SamplingRate.NumericValue * 2; % 2 sec
states.Response(index_invalid) = 0;
for idx_trial = 1:length(index_isright_start),
    KeyDown = unique(states.KeyDown((index_isright_start(idx_trial)-buffer_before):(index_isright_stop(idx_trial)+buffer_after)));
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
    
    states.Response((index_isright_start(idx_trial)):index_isright_stop(idx_trial)) = TrialResponse;
    
end

clear signal_loop
clear states_loop

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REMOVE COMMON NOISE USING A COMMON-AVERAGE REFERENCE FILTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
fprintf(1, '> Common average filtering signal \n');
% now calculate the common-average reference signal again, without those channels that have significant line-noise
fprintf(1,'[');
if isfield(parameters,'DeviceIDMaster')
    % determine the number of g.USBamps where each amp has 16 channels
    num_amps = round(size(signal,2) / 16);
    % for each of these amps
    for idx_amp = 1:num_amps,
        idx_low  = (idx_amp-1)*16+1;
        idx_high = (idx_amp-0)*16+0;

        % exclude the channels that had signifiant line-noise
        list_channels = intersect(clean_channels,idx_low:idx_high);
        % check if any channels are left and 
        if ~isempty(list_channels) && length(list_channels)>1
            % calculate the common average reference signal
            signal_mean = mean(signal(:,list_channels),2);
            % subtract the common average signal from each channel of this amp
            for idx_ch=list_channels
                signal(:,idx_ch) = signal(:,idx_ch) - signal_mean;
                fprintf(1,'.');
            end
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILTER OUT REMAINING LINE-NOISE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    signal(:,idx_channel) = signal_preliminary;
    fprintf(1,'.');
end
fprintf(1,'] done\n');

%% gaussian band based on chang-lab 
bands = {'theta', 'alpha', 'beta', 'gamma', 'high_gamma','broad_band'};
min_freqs = [4., 8., 15., 30., 70.,4.];
max_freqs = [7., 14., 29., 58., 150.,150];
gaus_filt_defs=struct;
for k=1:length(min_freqs)
    [cfs,sds]=get_filter_param_chang_lab(min_freqs(k),max_freqs(k));
    gaus_filt_defs(k).band=bands{k};
    gaus_filt_defs(k).cfs=cfs;
    gaus_filt_defs(k).sds=sds;
end 
ops.gaus_filt_defs=gaus_filt_defs;    
%% compute bandpass signal in different bands based on chang-lab method, see :
% Dichter, Benjamin K., Jonathan D. Breshears, Matthew K. Leonard, and Edward F. Chang. 2018.
% ?The Control of Vocal Pitch in Human Laryngeal Motor Cortex.? Cell 174 (1): 21?31.e9.
fprintf(1, '>> Extracting and decimating envelopes based on chang-lab method  \n');
band_idx=find(~strcmp(bands,'broad_band'));
decimation_factor = parameters.SamplingRate.NumericValue / ops.downsamplingrate;
ops.decimation_factor=decimation_factor;

signal_gaus_bands=struct;
for k=band_idx
    cfs=gaus_filt_defs(k).cfs;
    sds=gaus_filt_defs(k).sds;
    filter_bank={};
    for s=1:length(cfs)
        filter_bank{s}=gaussian_filter(transpose(signal(:,1)),parameters.SamplingRate.NumericValue,cfs(s),sds(s));
    end
      %
    signal_hilbert_band=nan*double(signal);
    signal_hilbert_band_zs=nan*double(signal);
    signal_hilbert_band_dec = nan*zeros(size(decimate(double(signal(:,1)),decimation_factor),1),size(signal,2));
    signal_hilbert_band_dec_zs = nan*zeros(size(decimate(double(signal(:,1)),decimation_factor),1),size(signal,2));
    %
    fprintf('>> Extracting %s envelope  \n',gaus_filt_defs(k).band);
    parfor idx_channel=1:size(signal,2)
        signal_hilbert_bp_cell = cellfun(@abs,hilbert_transform(double(transpose(signal(:,idx_channel))),parameters.SamplingRate.NumericValue,filter_bank),'UniformOutput',false);
        % mean across
        signal_hilbert_bp=mean(cell2mat(signal_hilbert_bp_cell),1);
        % zscore the signal with respect to entire experimental block
        signal_hilbert_bp_zs=zscore(signal_hilbert_bp);
        % save 
        signal_hilbert_band(:,idx_channel)=signal_hilbert_bp;
        signal_hilbert_band_zs(:,idx_channel)=signal_hilbert_bp_zs;
        signal_hilbert_band_dec(:,idx_channel) = decimate(signal_hilbert_bp,decimation_factor);
        signal_hilbert_band_dec_zs(:,idx_channel) = decimate(signal_hilbert_bp_zs,decimation_factor);
         % calculate decimation factor for 12 Hz sampling rate    
    end
    signal_gaus_bands(k).band=gaus_filt_defs(k).band;
    signal_gaus_bands(k).hilbert=signal_hilbert_band;
    signal_gaus_bands(k).hilbert_zs=signal_hilbert_band_zs;
    signal_gaus_bands(k).hilbert_dec=signal_hilbert_band_dec;
    signal_gaus_bands(k).hilbert_dec_zs=signal_hilbert_band_dec_zs;
end 
%% calculating broad_band spectrum
fprintf(1, '>> Extracting and decimating broadband envelope based on chang-lab method  \n');
band_idx=find(strcmp(bands,'broad_band'));
signal_broad_bands=struct;
idx=0;
for k=band_idx
    idx=idx+1;
    cfs=gaus_filt_defs(k).cfs;
    sds=gaus_filt_defs(k).sds;
    filter_bank={};
    for s=1:length(cfs)
        filter_bank{s}=gaussian_filter(transpose(signal(:,1)),parameters.SamplingRate.NumericValue,cfs(s),sds(s));
    end
    %
    signal_hiblert_broad_band=cell(size(filter_bank));
    signal_hiblert_broad_band_dec=cell(size(filter_bank));
    %
    fprintf('>> Extracting %s envelope  \n',gaus_filt_defs(k).band);
    parfor idx_channel=1:size(signal,2)
        signal_hilbert_bp_cell = cellfun(@abs,hilbert_transform(double(transpose(signal(:,idx_channel))),parameters.SamplingRate.NumericValue,filter_bank),'UniformOutput',false);
       
        signal_hilbert_broadband=cell2mat(signal_hilbert_bp_cell);
        
        signal_hilbert_bp_dec_cell=cellfun(@(x) decimate(x,decimation_factor),signal_hilbert_bp_cell,'uni',false);
        signal_hilbert_broadband_dec=cell2mat(signal_hilbert_bp_dec_cell);
        % save 
        signal_hiblert_broad_band{:,idx_channel}=signal_hilbert_broadband;
        signal_hiblert_broad_band_dec{:,idx_channel}=signal_hilbert_broadband_dec;
        fprintf(1,'%d ',idx_channel);
    end
    signal_broad_bands(idx).band=gaus_filt_defs(k).band;
    signal_broad_bands(idx).hilbert=signal_hiblert_broad_band;
    signal_broad_bands(idx).hilbert_dec=signal_hiblert_broad_band_dec;
end
%% downsample states 
%% down-sample signal
fprintf(1, '>> Down-sampling states \n');

%     % decimate StimulusCode to 12 Hz
states.StimType_ds = downsample(double(states.StimType),ops.decimation_factor);
states.WordType_ds = downsample(double(states.WordType),ops.decimation_factor);
states.IsRight_ds  = downsample(double(states.IsRight),ops.decimation_factor);
states.Response_ds = downsample(double(states.Response),ops.decimation_factor);
states.StimulusCodeDownsample=downsample(double(states.StimulusCode),ops.decimation_factor);

%% 
ops_out=ops;
ops_out.signal_broad_bands=signal_broad_bands;
ops_out.signal_gaus_bands=signal_gaus_bands;
ops_out.states=states;
ops_out.parameters=parameters;
fprintf('done \n');
end 