function [dataout]=filter_channels_MGH_NSX_v5(datafile,varargin)
p=inputParser();
addParameter(p, 'op_info', struct);
addParameter(p, 'aux_data', []);
addParameter(p, 'timing', []);
addParameter(p,'bw60Hz',0.6);
addParameter(p,'frac',0.5);% see channel_selection_from_60Hz_noise.m
addParameter(p,'min_nchannels',10);
addParameter(p,'thresh60Hz',5);
addParameter(p,'hpcutoff',0.5); % hz --> Rad/(sr/2)
addParameter(p,'hporder',4);
addParameter(p,'notchbw',1);
addParameter(p,'notchfreqs',[60, 120, 180, 240]);
addParameter(p,'fsDownsample',128);
addParameter(p,'steps',{'highpass', 'global_car', 'notch','IEDremovalFromAux','High60HzFromAux','VisualInspectFromAux','removeOutliers'})
parse(p, varargin{:});
ops = p.Results;

% spec for analysis 
num_of_chan=size(datafile.Data,1);
sr=datafile.MetaTags.SamplingFreq;
ops.sr=sr;
ops.subject_id=datafile.MetaTags.subject_id;
ops.session_name=datafile.MetaTags.session_name;
ops.file_name=datafile.MetaTags.Filename;
ops.file_path=datafile.MetaTags.FilePath;
ops.full_path=datafile.MetaTags.fullPath;
ops.electrode_numbers = 1:num_of_chan;
ops.elecnames = datafile.ElectrodesInfo.Label';
ops.elec_clinic_info=datafile.ElectrodesInfo.clinical_info_table;
ops.elecids = [1:length(ops.elecnames)]';
ops.ecog_channels_labels=ops.elecnames;
ops.perfromed_steps={};
ops.step_ops={};
%% 
func_col=@(x) reshape(x,[],1);
%% 
signal = double((datafile.Data)');
% set some color schemas: 
color_sequence={{'Blues'   },{'YlOrRd'  },{'BuGn'    },{'Greys'   },{'BuPu'    },{'GnBu'    },{'Oranges' },{'Greens'  },{'OrRd'    },{'Blues'   },{'PuBu'    },{'PuBuGn'  },...
                {'PuRd'    },{'Blues'   },{'Purples' },{'Greens'  },{'RdPu'    },{'YlGn'    },{'Reds'    },{'YlGnBu'  },{'YlOrBr'  },{'YlOrRd'  },{'BuGn'    },{'BuPu'    },...
                {'GnBu'    },{'Greys'   },{'Oranges' },{'OrRd'    },{'PuBu'    },{'PuBuGn'  },{'PuRd'    },{'Purples' },{'RdPu'    },{'Reds'    },{'YlGn'    },{'YlGnBu'  },{'YlOrBr'  }};
% undifferentiated 
col_inf=inferno(floor(.8*size(signal,2)));
col_vir=viridis(floor(.8*size(signal,2)));
colors_base=[col_vir(1:floor(size(signal,2)/2),:);col_inf(1:(floor(size(signal,2)/2)+1),:)];
% uni polar
ch_tags=extract(ops.elecnames,lettersPattern);
[ecog_tags,~,~]=unique(ch_tags,'stable');
ecog_loc_for_tag=cellfun(@(x) find(contains(ops.elecnames,x))', ecog_tags,'uni',false);
colors_map=(arrayfun(@(x) flipud(cbrewer('seq',color_sequence{x}{1},floor(1.5*length(ecog_loc_for_tag{x})),'linear')),[1:size(ecog_loc_for_tag,1)]','uni',false));
colors_unipolar=cell2mat(cellfun(@(x,y) y(1:size(x,2),:),ecog_loc_for_tag,colors_map,'uni',false));
elec_color_info=arrayfun(@(x) [ops.elecnames(x),{colors_unipolar(x,:)}],1:numel(ops.elecnames),'uni',false)';
elec_color_info=vertcat(elec_color_info{:});
ops.elec_color_info=elec_color_info;
% bipolar 
ch_num=cellfun(@str2num,extract(ops.elecnames,digitsPattern));
ch_pairs_id=[find(diff(ch_num)==1)+1,find(diff(ch_num)==1)];
ch_pairs_tags=ch_tags(ch_pairs_id);
ch_pair_names=ops.elecnames(ch_pairs_id);
assert(all(strcmp(ch_pairs_tags(:,1),ch_pairs_tags(:,2))));
bip_loc_for_tag=cellfun(@(x) find(contains(ch_pairs_tags(:,1),x))', ecog_tags,'uni',false);
colors_map=(arrayfun(@(x) flipud(cbrewer('seq',color_sequence{x}{1},floor(1.5*length(bip_loc_for_tag{x})),'linear')),[1:size(bip_loc_for_tag,1)]','uni',false));
colors_bipolar=cell2mat(cellfun(@(x,y) y(1:size(x,2),:),bip_loc_for_tag,colors_map,'uni',false));
ops.bip_elecnames=arrayfun(@(x) [ch_pair_names{x,1},'-',ch_pair_names{x,2}],1:size(ch_pair_names,1),'uni',false)';
bip_elec_color_info=arrayfun(@(x) [ops.bip_elecnames(x),{colors_bipolar(x,:)}],1:numel(ops.bip_elecnames),'uni',false)';
bip_elec_color_info=vertcat(bip_elec_color_info{:});
ops.bip_elec_color_info=bip_elec_color_info;

%% check whether to do IED chan removal 
if any(strcmp('IEDremovalFromAux',ops.steps))
    fprintf('\n reading IED events and channels from auxillary data...\n');
    aux_chan_deselect=erase(ops.aux_data.filt_ops.ecog_channels_labels(ops.aux_data.filt_ops.ecog_channels_IED_deselected),'_');
    [C,ia,ib]=intersect(ops.elecnames,aux_chan_deselect,'stable');
    assert(all(ia==cell2mat(cellfun(@(x) find(ismember(ops.elecnames,x)),aux_chan_deselect,'uni',false))));
    ops.ecog_channels_IED_deselected=ia;
    try
    ops.IED_results=ops.aux_data.filt_ops.IED_results;
    end 
    ops.perfromed_steps=[ops.perfromed_steps;'IEDremovalFromAux'];
    ops.steps = setdiff(ops.steps, 'IEDremovalFromAux');
elseif any(strcmp('IEDremoval',ops.steps))
    fprintf('\n extracting IED events ...\n');
    ops=find_epileptic_elec(signal,ops);
    ops.perfromed_steps=[ops.perfromed_steps;'IEDremoval'];
    ops.steps = setdiff(ops.steps, 'IEDremoval');
elseif any(strcmp('IEDremovalClinical',ops.steps))
    fprintf('\n reading IED electrodes from clinical data ...\n');
    ops.ecog_channels_IED_deselected=find([0*ops.elec_clinic_info.ictal_activity + ops.elec_clinic_info.interictal_activity]>0);
    ops.perfromed_steps=[ops.perfromed_steps;'IEDremovalClinical'];
    ops.steps = setdiff(ops.steps, 'IEDremovalClinical');
else
    ops.ecog_channels_IED_deselected=[];
    ops.IED_results=[];
end 
%% find bad channels based on 60Hz noise 
if any(strcmp('High60Hz',ops.steps))
    fprintf('\n Detecting good/bad channels with 60 Hz noise...\n');
    [good_channels,bad_channels, High60Hz_ops] = channel_selection_from_60Hz_noise(signal, ops.sr);
    ops.ecog_channels_noise_5std_deselected=bad_channels;
    ops.ecog_channels_noise_5std_selected=good_channels;
    ops.perfromed_steps=[ops.perfromed_steps;'High60Hz'];
    ops.step_ops.('High60Hz')=High60Hz_ops;
    ops.steps = setdiff(ops.steps, 'High60Hz');    
elseif any(strcmp('HighShank60Hz', ops.steps))
    fprintf('\n Detecting good/bad channels with 60 Hz noise on each shank...\n');
elseif any(strcmp('High60HzFromAux',ops.steps))
    fprintf('\n reading good/bad channels with 60 Hz noise from auxiliary data...\n');
    aux_chan_deselect=erase(ops.aux_data.filt_ops.ecog_channels_labels(ops.aux_data.filt_ops.ecog_channels_noise_5std_deselected),'_');
    [C,ia,ib]=intersect(ops.elecnames,aux_chan_deselect,'stable');
    assert(all(ia==cell2mat(cellfun(@(x) find(ismember(ops.elecnames,x)),aux_chan_deselect,'uni',false))));
    ops.ecog_channels_noise_5std_deselected=ia;
    ops.steps = setdiff(ops.steps, 'High60HzFromAux');
    ops.perfromed_steps=[ops.perfromed_steps;'High60HzFromAux'];
end
%% visual inpection 
if any(strcmp('VisualInspectFromAux',ops.steps))
    fprintf('\n reading visually inspected channels from auxilary data...\n');
    aux_chan_user_deselect=erase(ops.aux_data.filt_ops.ecog_channels_labels(ops.aux_data.filt_ops.ecog_channels_user_deselect),'_');
    [C,ia,ib]=intersect(ops.elecnames,aux_chan_user_deselect,'stable');
    assert(all(ia==cell2mat(cellfun(@(x) find(ismember(ops.elecnames,x)),aux_chan_user_deselect,'uni',false))));
    ops.ecog_channels_user_deselect=ia;
    ops.steps = setdiff(ops.steps, 'VisualInspectFromAux');    
    ops.perfromed_steps=[ops.perfromed_steps;'VisualInspectFromAux'];
elseif any(strcmp('VisualInspect',ops.steps))
    fprintf('\n perfoming visual inspection...\n');
    ecog_channels_selected = setdiff(ops.elecids,union(func_col(ops.ecog_channels_noise_5std_deselected),func_col(ops.ecog_channels_IED_deselected)),'stable');
    valid_channels=ops.elecids*nan;
    valid_channels(ecog_channels_selected)=1;
    plot_channels(signal,ops.sr,ops.elec_color_info(:,2),ops.elec_color_info(:,1),valid_channels)
    waitfor(findobj('type','figure','number',1));
    prompt='\n additional channels to remove? format :[1,2] ';
    ops.ecog_channels_user_deselect=input(prompt);
    ops.steps = setdiff(ops.steps, 'VisualInspect');    
else 
    ops.ecog_channels_user_deselect=[];
end 
%% get a sumamry 
ops.ecog_channels_deselected=unique(vertcat(func_col(ops.ecog_channels_noise_5std_deselected),func_col(ops.ecog_channels_IED_deselected),func_col(ops.ecog_channels_user_deselect)));
ops.ecog_channels_selected = setdiff(ops.elecids,ops.ecog_channels_deselected,'stable');
valid_channels=ops.elecids*nan;
valid_channels(ops.ecog_channels_selected)=1;
ops.ecog_valid_chan_ids=valid_channels;
%% go through the rest of the steps 
for s_id=1:numel(ops.steps)
    switch ops.steps{s_id}
        case 'highpass'
            fprintf('\n High-pass filtering the signal...\n');
            [signal,hp_ops] = hp_filt(signal, ops.sr, 'order', ops.hporder, 'cutoff', ops.hpcutoff);
            ops.step_ops.('highpass')=hp_ops;
            ops.perfromed_steps=[ops.perfromed_steps;'hp_filter'];
        case 'global_car'
            fprintf('\n global common source averaging based on good electrodes...\n');
            overall_mean=mean(signal(:,ops.ecog_channels_selected),2);
            signal=signal-repmat(overall_mean,1,size(signal,2));
            ops.perfromed_steps=[ops.perfromed_steps;'global_car'];
        case 'shank_car' % this one needs additional workd 
            fprintf('\n common source averaging based on good electrodes in each shank...\n');
            [signal,ops]=comm_src_ave(signal,ops);
            ops.perfromed_steps=[ops.perfromed_steps;'shank_common_source_removal'];
        case 'notch'
            fprintf('\n Notch filtering line noise...\n');
            [signal,notch_ops] = notch_filt(signal, sr);
            ops.step_ops.('notch')=notch_ops;
            ops.perfromed_steps=[ops.perfromed_steps;'notch'];
        case 'removeOutliers'
            fprintf('');
        otherwise
            error('\n step %s is not recognized \n', ops.steps{s_id})
    end 
end 

plot_channels(signal,ops.sr,ops.elec_color_info(:,2),ops.elec_color_info(:,1),ops.ecog_valid_chan_ids)
waitfor(findobj('type','figure','number',1));

%% create bipolar names
fprintf('\n creating names for bipolar pairs...\n');
ops.valid_elecnames=ops.elecnames(ops.ecog_channels_selected);
valid_ch_num=cellfun(@str2num,extract(ops.valid_elecnames,digitsPattern));
valid_ch_pairs_id=[find(diff(valid_ch_num)==1)+1,find(diff(valid_ch_num)==1)];
ch_pairs_tags=extract(ops.valid_elecnames,lettersPattern);
valid_ch_pairs_tags=ch_pairs_tags(valid_ch_pairs_id);
valid_bip_loc_for_tag=cellfun(@(x) find(contains(valid_ch_pairs_tags(:,1),x))', ecog_tags,'uni',false);
assert(all(strcmp(valid_ch_pairs_tags(:,1),valid_ch_pairs_tags(:,2))));
% fill ops 
ops.bip_ch_label_valid=ops.valid_elecnames(valid_ch_pairs_id);
ops.bip_ch_id_valid=cellfun(@(x) find(ismember(ops.elecnames,x)),ops.bip_ch_label_valid);
ops.bip_ch_label_valid_grp=cellfun(@(x) ops.bip_ch_label_valid(x,:),valid_bip_loc_for_tag,'uni',false);
ops.bip_ch_id_valid_grp=cellfun(@(x) ops.bip_ch_id_valid(x,:),valid_bip_loc_for_tag,'uni',false);
ops.bip_elecnames_valid=arrayfun(@(x) [ops.bip_ch_label_valid{x,1},'-',ops.bip_ch_label_valid{x,2}],1:size(ops.bip_ch_label_valid,1),'uni',false)';

bip_colors_valid=cellfun(@(x) ops.bip_elec_color_info{ismember(ops.bip_elec_color_info(:,1),x),2}, ops.bip_elecnames_valid,'uni',false);
ops.bip_elec_color_info_valid=horzcat(ops.bip_elecnames_valid,bip_colors_valid);
%% create a biopolar version of the signal
fprintf(1, '\n creating a bipolar version of the signal... \n');
signal_bipolar=make_bipolar_signal(signal,ops);
ops.preproc_step=[ops.perfromed_steps;'created_bipolar_signal'];
plot_channels(signal_bipolar,ops.sr,ops.bip_elec_color_info_valid(:,2),ops.bip_elec_color_info_valid(:,1),ones(size(signal_bipolar,2),1),'timing',ops.timing)
waitfor(findobj('type','figure','number',1));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRACT AND DECIMATE BANDPOWER FOR EACH FREQUENCY BAND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute high gamma based on chang-lab method, see :
% Dichter, Benjamin K., Jonathan D. Breshears, Matthew K. Leonard, and Edward F. Chang. 2018.
% ?The Control of Vocal Pitch in Human Laryngeal Motor Cortex.? Cell 174 (1): 21?31.e9.
% https://github.com/bendichter/process_ecog/blob/master/ecog/signal_processing/hilbert_transform.py
fprintf(1, '\n Extracting unipolar high gamma envelope based on gaussian filtering...  \n');
fprintf(1, '\n defining gaussian  filter bank... \n');
f_gamma_low=70;
f_gamma_high=150;
[cfs,sds]=get_filter_param_chang_lab(f_gamma_low,f_gamma_high);
% add a correction to stay away from 60 hz
cfs(1)=73.0;

gaussian_ops.cfs=cfs;
gaussian_ops.sds=sds;
%filter_bank=cell(numel(cfs),1);
filter_bank={};
for s=1:length(cfs)
    filter_bank{s}=gaussian_filter(transpose(signal),ops.sr,cfs(s),sds(s));
end
% unipolar gaussian filtering 
signal_hiblert_bp=nan*signal;
signal_hiblert_zs=nan*signal_hiblert_bp;
fprintf(1, '\n computing hilbert tranform  \n');
pbar=ProgressBar(size(signal,2));
for kk=1:size(signal,2)
    signal_hilbert_bp_all = cellfun(@abs,hilbert_transform((transpose(signal(:,kk))),ops.sr,filter_bank),'UniformOutput',false);
    signal_hilbert_bp_all=cell2mat(signal_hilbert_bp_all);
    signal_hiblert_bp(:,kk)=transpose(mean(signal_hilbert_bp_all,1));
    signal_hilbert_zs(:,kk) = zscore(transpose(mean(signal_hilbert_bp_all,1)));
    pbar.step([],[],[]);
end 

% decimation
fprintf(1, '\n downsampling \n');
signal_hiblert_bp=resample(signal_hiblert_bp,ops.fsDownsample,ops.sr);
signal_hilbert_zs=resample(signal_hilbert_zs,ops.fsDownsample,ops.sr);
plot_channels(signal_hiblert_bp,ops.sr,ops.elec_color_info(:,2),ops.elec_color_info(:,1),valid_channels)
waitfor(findobj('type','figure','number',1));

ops.step_ops.('gaussian_hilbert')=gaussian_ops;
ops.perfromed_steps=[ops.perfromed_steps;'extracting_unipolar_gamma_envelope_gaussain'];
%% SNH method using a bandpass butterworth electrode 
fprintf(1, '\n Extracting unipolar high gamma envelope based on band-pass filtering  \n');
fprintf(1, '\n perfoming bandpass filtering and reampling  \n');
cutoffs = [70, 140];
order = 6;

bandpass_ops.cutoffs=cutoffs;
bandpass_ops.order=order;

envelopes = bandpass_envelopes(signal, ops.sr, cutoffs, order);
fprintf(1, '\n downsampling... \n');
envelopes=resample(envelopes,ops.fsDownsample,ops.sr);
% do some checks 
plot_channels(envelopes,ops.sr,ops.elec_color_info(:,2),ops.elec_color_info(:,1),valid_channels)
waitfor(findobj('type','figure','number',1));

ops.step_ops.('bandpass_hilbert')=bandpass_ops;
ops.perfromed_steps=[ops.perfromed_steps;'extracting_unipolar_gamma_envelope_bandpass'];
%% bipolar gaussian filtering 
fprintf(1, '\n Extracting bipolar high gamma envelope based on gaussian filtering  \n');
fprintf(1, '\n computing hilbert tranform  \n');
signal_bipolar_hiblert_bp=nan*signal_bipolar;
signal_bipolar_hiblert_zs=nan*signal_bipolar_hiblert_bp;

pbar=ProgressBar(size(signal_bipolar,2));
for kk=1:size(signal_bipolar,2)
    signal_bipolar_hilbert_bp_all = cellfun(@abs,hilbert_transform((transpose(signal_bipolar(:,kk))),ops.sr,filter_bank),'UniformOutput',false);
    signal_bipolar_hilbert_bp_all=cell2mat(signal_bipolar_hilbert_bp_all);
    signal_bipolar_hiblert_bp(:,kk)=transpose(mean(signal_bipolar_hilbert_bp_all,1));
    signal_bipolar_hiblert_zs(:,kk) = zscore(transpose(mean(signal_bipolar_hilbert_bp_all,1)));
    pbar.step([],[],[]);
end 

signal_bipolar_hiblert_bp=resample(signal_bipolar_hiblert_bp,ops.fsDownsample,ops.sr);
signal_bipolar_hiblert_zs=resample(signal_bipolar_hiblert_zs,ops.fsDownsample,ops.sr);

% plot throughout 
plot_channels(signal_bipolar_hiblert_zs,ops.sr,ops.bip_elec_color_info_valid(:,2),ops.bip_elec_color_info_valid(:,1),ones(size(signal_bipolar_hiblert_zs,2),1),'timing',ops.timing)
waitfor(findobj('type','figure','number',1));
ops.perfromed_steps=[ops.perfromed_steps;'extracting_bipolar_gamma_envelope_gaussain'];
%% 
fprintf(1, '\n Extracting bipolar high gamma envelope based on band-pass filtering  \n');
fprintf(1, '\n perfoming bandpass filtering and reampling  \n');
envelopes_bipolar = bandpass_envelopes(signal_bipolar, ops.sr, cutoffs, order);
envelopes_bipolar=resample(envelopes_bipolar,ops.fsDownsample,ops.sr);

plot_channels(envelopes_bipolar,ops.sr,ops.bip_elec_color_info_valid(:,2),ops.bip_elec_color_info_valid(:,1),ones(size(signal_bipolar_hiblert_zs,2),1),'timing',ops.timing)
waitfor(findobj('type','figure','number',1));
ops.perfromed_steps=[ops.perfromed_steps;'extracting_unipolar_gamma_envelope_bandpass'];
%% do oulier rejection based on 90%: 

if any(strcmp('removeOutliers',ops.steps))
    fprintf('\n removing outliers in the envelopes and HG hilbert...\n');
    fprintf('\n bandpass envelopes ...\n');
    [envelopes, outlier_ops] = envelope_outliers(envelopes, ops.sr);
    ops.step_ops.('removeOutliers_envelope')=outlier_ops;
    
    [envelopes_bipolar, outlier_ops] = envelope_outliers(envelopes_bipolar, ops.sr);
    ops.step_ops.('removeOutliers_envelopes_bipolar')=outlier_ops;
    
    fprintf('\n gaussain filtered signal ...\n');
    [signal_bipolar_hiblert_bp, outlier_ops] = envelope_outliers(signal_bipolar_hiblert_bp, ops.sr);
    ops.step_ops.('removeOutliers_signal_bipolar_hilbert_bp')=outlier_ops;
    
    [signal_bipolar_hiblert_zs, outlier_ops] = envelope_outliers(signal_bipolar_hiblert_zs, ops.sr);
    ops.step_ops.('removeOutliers_signal_bipolar_hilbert_zs')=outlier_ops;
    
    
    [signal_hiblert_bp, outlier_ops] = envelope_outliers(signal_hiblert_bp, ops.sr);
    ops.step_ops.('removeOutliers_signal_hiblert_bp')=outlier_ops;
    
    [signal_hilbert_zs, outlier_ops] = envelope_outliers(signal_hilbert_zs, ops.sr);
    ops.step_ops.('removeOutliers_signal_bipolar_hilbert_zs')=outlier_ops;
    
    
    ops.perfromed_steps=[ops.perfromed_steps;'removeOutliers'];
    ops.steps = setdiff(ops.steps, 'removeOutliers');

end 

%% save data out
dataout=datafile;
dataout.signal_hilbert_zs_decimated=signal_hilbert_zs';
dataout.signal_hilbert_decimated=signal_hiblert_bp';
dataout.evelopes=envelopes';
dataout.signal_bipolar_hilbert_zs_decimated=signal_bipolar_hiblert_zs';
dataout.signal_bioplar_hilbert_decimated=signal_bipolar_hiblert_bp';
dataout.evelopes_bipolar=envelopes_bipolar';
dataout=rmfield(dataout,'Data');
% remove data from aux
if ~isempty(ops.aux_data)
    aux_fields=fieldnames(ops.aux_data);
    remove_fields=aux_fields(contains(aux_fields,'_dec'));
    for x=remove_fields.'
        ops.aux_data.(x{1})=[];
    end 
end 
dataout.ops= ops;
fprintf('\n saving results ...\n');
save(dataout.MetaTags.fullPath,'dataout','-v7.3');
fprintf('\n done! file is saved at %s \n',dataout.MetaTags.fullPath);
end
