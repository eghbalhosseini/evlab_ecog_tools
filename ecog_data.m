
classdef ecog_data < dynamicprops
    %ECOG_DATA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %% data
        elec_data_dec
        elec_data_zs_dec
        bip_elec_data_dec
        bip_elec_data_zs_dec
        envelopes_dec
        bip_envelopes_dec
        %% timing 
        trial_timing_dec
        trial_timing
        %% info 
        info
        events_table
        subject_id
        file_name
        file_path
        experiment
        % 
        save_name
        save_path
        % general trial info 
        trial_type
        session_name
        modality_type
        % filtering info 
        filt_ops
        % frequency info 
        sample_freq
        decimation_freq
        modality
        modality_list
        %% channel labels 
        elec_ch_label
        elec_ch
        elec_ch_with_IED
        elec_ch_with_noise
        elec_ch_user_deselect
        elec_ch_clean
        elec_ch_valid
        % bipolar channel labels 
        bip_ch_valid_grp
        bip_ch_label_valid_grp
        bip_ch_valid
        bip_ch_label_valid
        %% anatomy
        lh_pial
        rh_pial
        elec_ch_pos_anat
        bip_ch_pos_anat
        % fsaverage 
        elec_ch_pos_fs_ave
        bip_ch_pos_fs_ave
        % mni 
        elec_ch_pos_mni
        bip_ch_pos_mni
        % 
        elec_ch_pos_mni_cat12
        bip_ch_pos_mni_cat12
        % hcp 
        elec_ch_HCP_label
        elec_ch_HCP_weight
        bip_ch_HCP_label
        bip_ch_HCP_weight
        %% processing
        % analysis of lang responsive 
        s_vs_n_sig
        s_vs_n_p_ratio
        s_vs_n_ops
        sn_data_ops
        %% trial data 
        trial_data
    end
    
    methods
%% method for constructing the object 
        function obj = ecog_data(prep_data,...
        trial_timing_dec,events_table) % timing
            %% data 
            obj.elec_data_dec = prep_data.signal_hilbert_decimated;
            obj.elec_data_zs_dec = prep_data.signal_hilbert_zs_decimated;
            obj.bip_elec_data_dec = prep_data.signal_bioplar_hilbert_decimated;
            obj.bip_elec_data_zs_dec = prep_data.signal_bipolar_hilbert_zs_decimated;
            obj.envelopes_dec=prep_data.evelopes;
            obj.bip_envelopes_dec=prep_data.evelopes_bipolar;
            %% timing 
            obj.trial_timing_dec=trial_timing_dec;
            %% channel labels 
            obj.elec_ch_label=prep_data.ops.ecog_channels_labels;
            obj.elec_ch=prep_data.ops.elecids;
            obj.elec_ch_with_IED=prep_data.ops.ecog_channels_IED_deselected;
            obj.elec_ch_with_noise=prep_data.ops.ecog_channels_noise_5std_deselected;
            obj.elec_ch_user_deselect=prep_data.ops.ecog_channels_user_deselect;
            obj.elec_ch_clean=prep_data.ops.ecog_channels_selected;
            obj.elec_ch_valid=prep_data.ops.ecog_valid_chan_ids;
            % bipolar channel labels 
            obj.bip_ch_valid_grp=prep_data.ops.bip_ch_id_valid_grp;
            obj.bip_ch_label_valid_grp=prep_data.ops.bip_ch_label_valid_grp;
            obj.bip_ch_valid=prep_data.ops.bip_ch_id_valid;
            obj.bip_ch_label_valid=prep_data.ops.bip_ch_label_valid;
            %% info 
            obj.events_table=events_table;
            obj.subject_id=prep_data.ops.subject_id;
            obj.file_name=prep_data.ops.file_name;
            obj.file_path=prep_data.ops.file_path;
            obj.trial_type=trial_timing_dec(:,2); 
            obj.session_name=prep_data.ops.session_name;
            % frequency info
            obj.sample_freq=prep_data.ops.sr;
            obj.decimation_freq=prep_data.ops.fsDownsample;
            obj.filt_ops=prep_data.ops;
            
        end
        %% methods for denoising the data
        function obj=make_trials(obj)
            %METHOD1 Summary of this method goes here
            trial_keys={'key','string','elec_data_dec','elec_data_zs_dec','bip_elec_data_dec','bip_elec_data_zs_dec','envelope_dec','bip_envelope_dec'};
            %   Detailed explanation goes here
            fprintf(1, '\n splitting the data into trials \n');
            pbar=ProgressBar(size(obj.trial_timing_dec,1));
            for k=1:size(obj.trial_timing_dec,1)
                trial_time_tbl=obj.trial_timing_dec{k,1};
                trial_elec_data_dec=arrayfun(@(x) obj.elec_data_dec(:,floor(trial_time_tbl(x,:).start):ceil(trial_time_tbl(x,:).end)),1:size(trial_time_tbl,1),'uni',false)';
                trial_elec_data_zs_dec=arrayfun(@(x) obj.elec_data_zs_dec(:,floor(trial_time_tbl(x,:).start):ceil(trial_time_tbl(x,:).end)),1:size(trial_time_tbl,1),'uni',false)';
                trial_bip_elec_data_dec=arrayfun(@(x) obj.bip_elec_data_dec(:,floor(trial_time_tbl(x,:).start):ceil(trial_time_tbl(x,:).end)),1:size(trial_time_tbl,1),'uni',false)';
                trial_bip_elec_data_zs_dec=arrayfun(@(x) obj.bip_elec_data_zs_dec(:,floor(trial_time_tbl(x,:).start):ceil(trial_time_tbl(x,:).end)),1:size(trial_time_tbl,1),'uni',false)';
                envelope_dec=arrayfun(@(x) obj.envelopes_dec(:,floor(trial_time_tbl(x,:).start):ceil(trial_time_tbl(x,:).end)),1:size(trial_time_tbl,1),'uni',false)';
                bip_envelope_dec=arrayfun(@(x) obj.bip_envelopes_dec(:,floor(trial_time_tbl(x,:).start):ceil(trial_time_tbl(x,:).end)),1:size(trial_time_tbl,1),'uni',false)';
                obj.trial_data{k,1}=table(trial_time_tbl.key,trial_time_tbl.string,trial_elec_data_dec,trial_elec_data_zs_dec,trial_bip_elec_data_dec,trial_bip_elec_data_zs_dec,envelope_dec,bip_envelope_dec,'VariableNames',trial_keys);
                %obj.trial_data{k,2}=obj.trial_timing_dec{k,2};
                pbar.step([],[],[]);
            end
            pbar.release()
        end
        %% methods for making anatomical labels
        function obj=connect_anatomy(obj,ch_RAS_tbl)
            % first check the labels are correct 
            ch_RAS_tbl.Properties.VariableNames(contains(ch_RAS_tbl.Properties.VariableNames,'Var1'))={'label'};
            % find index of bipolar chan in table
            bip_ch_=cellfun(@(x) erase(x,'_'),obj.biop_ch_label_valid,'uni',false);
            bip_ch_idx=cell2mat(cellfun(@(x) find(ismember(ch_RAS_tbl.label,x)),bip_ch_,'uni',false));
            % find location of bipolar channels 
            func=@(X) arrayfun(@(x) mean([X(bip_ch_idx(x,2)),X(bip_ch_idx(x,1))]),1:size(bip_ch_idx,1),'uni',false)';
            values=ch_RAS_tbl(:,~ismember(ch_RAS_tbl.Properties.VariableNames,'label'));
            values_comb=cellfun(@(X) func(values.(X)),values.Properties.VariableNames,'uni',false);
            % make bipolar names 
            bip_names=arrayfun(@(x) [bip_ch_{x,1},'-',bip_ch_{x,2}],1:size(obj.biop_ch_label_valid,1),'uni',false)';
            bip_ch_table=cell2table(horzcat(bip_names,[values_comb{:}]),'VariableNames',ch_RAS_tbl.Properties.VariableNames);
            obj.bip_ch_pos_anat=bip_ch_table;
            % do the same for uni-polar 
            %assert(length(obj.elec_ch_label)==size(ch_RAS_tbl,1))
            ch_=cellfun(@(x) erase(x,'_'),obj.elec_ch_label,'uni',false);
            ch_idx=cell2mat(cellfun(@(x) find(ismember(ch_RAS_tbl.label,x)),ch_,'uni',false));
            func=@(X) arrayfun(@(x) X(x), ch_idx,'uni',false);
            values_comb=cellfun(@(X) func(values.(X)),values.Properties.VariableNames,'uni',false);
            ch_table=cell2table(horzcat(ch_,[values_comb{:}]),'VariableNames',ch_RAS_tbl.Properties.VariableNames);
            obj.elec_ch_pos_anat=ch_table;
        end 
        function obj=connect_fsaverage(obj,ch_RAS_tbl_uni,ch_RAS_tbl_bip)
            % first check the labels are correct 
            ch_RAS_tbl_uni.Properties.VariableNames(contains(ch_RAS_tbl_uni.Properties.VariableNames,'Var1'))={'label'};
            ch_RAS_tbl_bip.Properties.VariableNames(contains(ch_RAS_tbl_bip.Properties.VariableNames,'name'))={'label'};
            % find index of bipolar chan in table
            bip_ch_=cellfun(@(x) erase(x,'_'),obj.biop_ch_label_valid,'uni',false);
            bip_ch_diff=cellfun(@(x,y) strcat(x,'-',y),bip_ch_(:,1),bip_ch_(:,2),'uni',false)
            bip_ch_idx=cell2mat(cellfun(@(x) find(ismember(ch_RAS_tbl_bip.label,x)),bip_ch_diff,'uni',false));
            assert(size(bip_ch_idx,1)==size(bip_ch_diff,1))
            % find location of bipolar channels 
            func=@(X) arrayfun(@(x) mean([X(bip_ch_idx(x,2)),X(bip_ch_idx(x,1))]),1:size(bip_ch_idx,1),'uni',false)';
            values=ch_RAS_tbl(:,~ismember(ch_RAS_tbl.Properties.VariableNames,'label'));
            values_comb=cellfun(@(X) func(values.(X)),values.Properties.VariableNames,'uni',false);
            % make bipolar names 
            bip_names=arrayfun(@(x) [bip_ch_{x,1},'-',bip_ch_{x,2}],1:size(obj.biop_ch_label_valid,1),'uni',false)';
            bip_ch_table=cell2table(horzcat(bip_names,[values_comb{:}]),'VariableNames',ch_RAS_tbl.Properties.VariableNames);
            obj.bip_ch_pos_anat=bip_ch_table;
            % do the same for uni-polar 
            %assert(length(obj.elec_ch_label)==size(ch_RAS_tbl,1))
            ch_=cellfun(@(x) erase(x,'_'),obj.elec_ch_label,'uni',false);
            ch_idx=cell2mat(cellfun(@(x) find(ismember(ch_RAS_tbl.label,x)),ch_,'uni',false));
            func=@(X) arrayfun(@(x) X(x), ch_idx,'uni',false);
            values_comb=cellfun(@(X) func(values.(X)),values.Properties.VariableNames,'uni',false);
            ch_table=cell2table(horzcat(ch_,[values_comb{:}]),'VariableNames',ch_RAS_tbl.Properties.VariableNames);
            obj.elec_ch_pos_anat=ch_table;
        end
        function obj=connect_mni_brain(obj)

            fprintf('doing electrode mapping onto mni\n');
            clinical_info=obj.filt_ops.elec_clinic_info;
            missing_elec=setdiff(lower(obj.elec_ch_label),lower(clinical_info.label));
            assert(isempty(missing_elec),sprintf('%s, missing: %s',obj.subject_id,strjoin(missing_elec,', ')));
            assert(all(ismember(lower(obj.elec_ch_label),lower(clinical_info.label))));
            chn_pos_in_mni={};
            chn_label_HCP={};
            chn_weigth_HCP={};
            for k_ch=1:length(obj.elec_ch_label)
                ch_id=obj.elec_ch_label{k_ch};
                chan_loc=find(ismember(lower(clinical_info.label),lower(ch_id)));
                assert(length(chan_loc)==1);
                ch_info=clinical_info(chan_loc,:);
                chn_pos_in_mni{k_ch,1}=[ch_info.mni_linear_x,ch_info.mni_linear_y,ch_info.mni_linear_z];
                chn_label_HCP{k_ch,1}=[ch_info.HCPMMP1_label_1];
                chn_weigth_HCP{k_ch,1}=[ch_info.HCPMMP1_weight_1];
            end 
            obj.elec_ch_pos_mni=chn_pos_in_mni;
            obj.elec_ch_HCP_label=chn_label_HCP;
            obj.elec_ch_HCP_weight=chn_weigth_HCP;
            % bipolar - this is estimate 
            bip_chn_pos_in_mni={};
            bip_chn_label_HCP={};
            bip_chn_weigth_HCP={};
            for k_ch=1:length(obj.bip_ch_label_valid)
                bip_1=obj.bip_ch_label_valid{k_ch,1};
                bip_2=obj.bip_ch_label_valid{k_ch,2};
                bip_1_chan_loc=find(ismember(lower(clinical_info.label),lower(bip_1)));
                assert(length(bip_1_chan_loc)==1);
                bip_2_chan_loc=find(ismember(lower(clinical_info.label),lower(bip_2)));
                assert(length(bip_2_chan_loc)==1);
                bip_1_ch_info=clinical_info(bip_1_chan_loc,:);
                bip_2_ch_info=clinical_info(bip_2_chan_loc,:);
                
                bip_chn_pos_in_mni{k_ch,1}=mean([bip_1_ch_info.mni_linear_x,bip_1_ch_info.mni_linear_y,bip_1_ch_info.mni_linear_z;
                    bip_2_ch_info.mni_linear_x,bip_2_ch_info.mni_linear_y,bip_2_ch_info.mni_linear_z],1);
                if strcmp(bip_1_ch_info.HCPMMP1_label_1,bip_2_ch_info.HCPMMP1_label_1)
                    bip_chn_label_HCP{k_ch,1}=[bip_1_ch_info.HCPMMP1_label_1];
                    bip_chn_weigth_HCP{k_ch,1}=[bip_1_ch_info.HCPMMP1_weight_1];
                else
                    bip_chn_label_HCP{k_ch,1}='N/A';
                    bip_chn_weigth_HCP{k_ch,1}=[nan];
                end 
            end 
            obj.bip_ch_pos_mni=bip_chn_pos_in_mni;
            obj.bip_ch_HCP_label=bip_chn_label_HCP;
            obj.bip_ch_HCP_weight=bip_chn_weigth_HCP;
            % add anatomy 
            icbm_file=load('icmb152_pial.mat');
            if not(isprop(obj,'icmb152_pial'))
                P = addprop(obj,'icmb152_pial');
            end
            obj.icmb152_pial=icbm_file.icmb152_pial;
        end 
        function obj=connect_mni_brain_va_cat12(obj)
            fprintf('doing electrode mapping onto mni via cat 12\n');
            clinical_info=obj.filt_ops.elec_clinic_info_MNI_cat12;
            missing_elec=setdiff(lower(obj.elec_ch_label),lower(clinical_info.label));
            assert(isempty(missing_elec),sprintf('%s, missing: %s',obj.subject_id,strjoin(missing_elec,', ')));
            assert(all(ismember(lower(obj.elec_ch_label),lower(clinical_info.label))));
            chn_pos_in_mni_cat12={};
            for k_ch=1:length(obj.elec_ch_label)
                ch_id=obj.elec_ch_label{k_ch};
                chan_loc=find(ismember(lower(clinical_info.label),lower(ch_id)));
                assert(length(chan_loc)==1);
                ch_info=clinical_info(chan_loc,:);
                chn_pos_in_mni_cat12{k_ch,1}=[ch_info.x,ch_info.y,ch_info.z];
            
            end 
            obj.elec_ch_pos_mni_cat12=chn_pos_in_mni_cat12;
            % bipolar - this is estimate 
            bip_chn_pos_in_mni_cat12={};
            for k_ch=1:length(obj.bip_ch_label_valid)
                bip_1=obj.bip_ch_label_valid{k_ch,1};
                bip_2=obj.bip_ch_label_valid{k_ch,2};
                bip_1_chan_loc=find(ismember(lower(clinical_info.label),lower(bip_1)));
                assert(length(bip_1_chan_loc)==1);
                bip_2_chan_loc=find(ismember(lower(clinical_info.label),lower(bip_2)));
                assert(length(bip_2_chan_loc)==1);
                bip_1_ch_info=clinical_info(bip_1_chan_loc,:);
                bip_2_ch_info=clinical_info(bip_2_chan_loc,:);

                bip_chn_pos_in_mni_cat12{k_ch,1}=mean([bip_1_ch_info.x,bip_1_ch_info.y,bip_1_ch_info.z;
                    bip_2_ch_info.x,bip_2_ch_info.y,bip_2_ch_info.z],1);
                
            end 
            obj.bip_ch_pos_mni_cat12=bip_chn_pos_in_mni_cat12;
        end 
        function obj=align_with_lang_atlas(obj,lang_atlas_dir)
            %fprintf('doing something\n')
            info_parc=niftiinfo([lang_atlas_dir,'/ROIS_NOV2020/Func_Lang_LHRH_SN220/allParcels_language.nii']);
            V_parc=niftiread([lang_atlas_dir,'/ROIS_NOV2020/Func_Lang_LHRH_SN220/allParcels_language.nii']);
            info_prb=niftiinfo([lang_atlas_dir,'/LanA/SPM/LanA_n806.nii']);
            V_prb=niftiread([lang_atlas_dir,'/LanA/SPM/LanA_n806.nii']);
            v_mni=niftiread([lang_atlas_dir,'/mni_icbm152_nlin_asym_09b/mni_icbm152_t1_tal_nlin_asym_09b_hires.nii']);
            info_mni=niftiinfo([lang_atlas_dir,'/mni_icbm152_nlin_asym_09b/mni_icbm152_t1_tal_nlin_asym_09b_hires.nii']);
            corner_x=[90:-2:-90];
            corner_y=[-126:2:90];
            corner_z=[-72:2:108];

            x_y_z_prac=[];
            x_y_z_prb=[];
            x_y_z_idx=[];
            x_y_z_lab={};
            x_y_z_cord=[];


            for id_=1:size(obj.elec_ch_pos_mni,1)
                xyz=obj.elec_ch_pos_mni{id_};
                label=obj.elec_ch_label{id_};
                [c,i_x]=min(abs(corner_x-xyz(1)));
                [c,i_y]=min(abs(corner_y-xyz(2)));
                [c,i_z]=min(abs(corner_z-xyz(3)));
                x_y_z_idx(id_,:)=[i_x,i_y,i_z];
                x_y_z_prac(id_,1)=V_parc(i_x,i_y,i_z);
                x_y_z_prb(id_,1)=V_prb(i_x,i_y,i_z);
                x_y_z_lab{id_}=label;
                x_y_z_cord(id_,:)=[corner_x(i_x),corner_y(i_y),corner_z(i_z)];
            end
            if not(isprop(obj,'elec_parc_idx'))
                P = addprop(obj,'elec_parc_idx');
            end
            if not(isprop(obj,'elec_parc_label'))
                P = addprop(obj,'elec_parc_label');
            end
            if not(isprop(obj,'elec_parc_coords'))
                P = addprop(obj,'elec_parc_coords');
            end
            if not(isprop(obj,'elec_parc'))
                P = addprop(obj,'elec_parc');
            end
            if not(isprop(obj,'elec_lana_prb'))
                P = addprop(obj,'elec_lana_prb');
            end
            %
            obj.elec_parc_label=x_y_z_lab;
            obj.elec_parc_idx=x_y_z_idx;
            obj.elec_parc_coords=x_y_z_cord;
            obj.elec_lana_prb=x_y_z_prb;
            obj.elec_parc=x_y_z_prac;
        end
        %% methods for averaging signal
        function cond_data=get_cond_resp(obj,condition)
            cond_id=find(ismember(obj.trial_type,condition));
            cond_data=obj.trial_data(cond_id);

        end
        function input_d_ave=get_average(obj,input_d,varargin)
            % work on a cell of tables
            p=inputParser();
            addParameter(p, 'dim', 2);
            parse(p, varargin{:});
            ops = p.Results;
            input_d_ave=input_d;
            for k=1:size(input_d,1)
                B=input_d{k};
                keys=B(:,ismember(B.Properties.VariableNames,'key'));
                strings=B(:,ismember(B.Properties.VariableNames,'string'));
                values=B(:,~ismember(B.Properties.VariableNames,{'key','string'}));
                values_ave=varfun( @(x) cellfun(@(y) squeeze(nanmean(y,ops.dim)), x,'uni',false), values,'OutputFormat','table');
                values_ave.Properties.VariableNames=values.Properties.VariableNames;
                input_d_ave{k}=[keys,strings,values_ave];
            end
        end
        function output_d=get_value(obj,input_d,varargin)
            % work on tables and cells
            p=inputParser();
            addParameter(p, 'key', 'word');
            addParameter(p, 'type', 'match'); % match or contain
            parse(p, varargin{:});
            ops = p.Results;
            if strcmp(ops.type,'match')
                func=@(x,y) ismember(x,y);
            else
                func=@(x,y) contains(x,y);
            end
            if istable(input_d)
                output_d=input_d(func(input_d.key,ops.key),:);
            elseif iscell(input_d)
                output_d=cellfun(@(x) x(func(x.key,ops.key),:),input_d,'uni',false);
            end
        end


    end
end

