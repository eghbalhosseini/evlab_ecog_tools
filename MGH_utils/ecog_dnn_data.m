classdef ecog_dnn_data < ecog_data
    %ECOG_DNN_DATA Summary of this class goes here
    properties
        % add properties for average signal for trial 
    end
    methods
        function dnn_obj = ecog_dnn_data(prep_data,trial_timing_dec,events_table) 

            dnn_obj@ecog_data(prep_data,...
                              trial_timing_dec,...
                              events_table)
        end
        function obj=add_sn_results(obj,sn_data)
            obj.s_vs_n_ops=sn_data.s_vs_n_ops;
            % get elec alignments:
            % uni
            sn_ch_label=sn_data.elec_ch_label;
            dnn_ch_label=obj.elec_ch_label;
            [~,~,ib_uni]=intersect(dnn_ch_label,sn_ch_label,'stable');
            % bip
            sn_bip_ch_label=cellfun(@(x,y) [x,'-',y],sn_data.bip_ch_label_valid(:,1),sn_data.bip_ch_label_valid(:,2),'uni',false);
            dnn_bip_ch_label=cellfun(@(x,y) [x,'-',y],obj.bip_ch_label_valid(:,1),obj.bip_ch_label_valid(:,2),'uni',false);
            [~,~,ib_bip]=intersect(dnn_bip_ch_label,sn_bip_ch_label,'stable');
            % modify s_vs_n_p_ratio table
            values=sn_data.s_vs_n_p_ratio(:,~contains(sn_data.s_vs_n_p_ratio.Properties.VariableNames,{'key'}));
            uni_values=values(:,~contains(values.Properties.VariableNames,{'bip'}));
            bip_values=values(:,contains(values.Properties.VariableNames,{'bip'}));
            dnn_uni_values=varfun(@(x) x{1}(ib_uni),uni_values,'outputformat','cell');
            dnn_bip_values=varfun(@(x) x{1}(ib_bip),bip_values,'outputformat','cell');
            obj.s_vs_n_p_ratio=cell2table(horzcat('s_vs_n_p_ratio',dnn_uni_values,dnn_bip_values),'VariableNames',horzcat('key',uni_values.Properties.VariableNames,bip_values.Properties.VariableNames));
            % modify s_vs_n_p_ratio table
            values=sn_data.s_vs_n_sig(:,~contains(sn_data.s_vs_n_sig.Properties.VariableNames,{'key'}));
            uni_values=values(:,~contains(values.Properties.VariableNames,{'bip'}));
            bip_values=values(:,contains(values.Properties.VariableNames,{'bip'}));
            dnn_uni_values=varfun(@(x) x{1}(ib_uni),uni_values,'outputformat','cell');
            dnn_bip_values=varfun(@(x) x{1}(ib_bip),bip_values,'outputformat','cell');
            obj.s_vs_n_sig=cell2table(horzcat('s_vs_n_sig',dnn_uni_values,dnn_bip_values),'VariableNames',horzcat('key',uni_values.Properties.VariableNames,bip_values.Properties.VariableNames));

        end 
%% method for creating regression data for neural_nlp
        function obj=make_stim_resp_mat(obj,varargin)
            p=inputParser();
            addParameter(p, 'words', 'all');
            parse(p, varargin{:});
            ops = p.Results;
            %get sentence id from events 
            aud_name=obj.events_table.final_audio_filename;
            aud_name=erase(aud_name,'.wav');
            aud_transcript=obj.events_table.final_audio_transcript;
            sent_aud_name=aud_name(contains(aud_name,'sentence'));
            sent_aud_transc=aud_transcript(contains(aud_name,'sentence'));
            % 
            condtion_flag='S';
            cond_resp=obj.get_cond_resp(condtion_flag);
            % identify and remove sentence repeats 
            % https://www.mathworks.com/matlabcentral/answers/197246-finding-duplicate-strings-in-a-cell-array-and-their-index
            cond_key='trial_start_end';
            trial_strings=cellfun(@(x) x.string{ismember(x.key,cond_key)},cond_resp,'uni',false);
            [D,P,X] = unique(trial_strings,'stable');
            assert(length(D)==200);
            % get only one of the repeats 
            uniq_cond_resp=cond_resp(P);
            uniq_str=cellfun(@(x) x.string{ismember(x.key,cond_key)},uniq_cond_resp,'uni',false);
            % make sure sentences are unique 
            cellfun(@(x,y) assert(strcmp(x,y)),unique(uniq_str,'stable'),uniq_str);
            % check if transcripts from events matches transcript from
            % trials 
            cellfun(@(x,y) assert(strcmp(x,y)),sent_aud_transc(P),uniq_str);
            % add sentence name to table 
            uniq_sent_aud_name=sent_aud_name(P);
            uniq_cond_resp_mod=arrayfun(@(x) [cell2table(repmat(uniq_sent_aud_name(x),size(uniq_cond_resp{x},1),1),"VariableNames",{'sentence_id'}),uniq_cond_resp{x}], 1:length(uniq_cond_resp),'uni',false)';
            % create average responses
            B=uniq_cond_resp;
            B=obj.get_average(B);
            word_data=obj.get_value(B,'key','word','type','contain');
            % 
            word_dat_mod=arrayfun(@(x) [cell2table(repmat(uniq_sent_aud_name(x),size(word_data{x},1),1),"VariableNames",{'sentence_name'}),word_data{x}], 1:length(word_data),'uni',false)';
            sent_dat=vertcat(word_dat_mod{:});
            % add sentence id for later use 
            sent_id=extract(sent_dat.sentence_name,digitsPattern);
            sent_id=num2cell(cellfun(@(x) str2num(x),sent_id));
            % add word id for later use
            word_id=extract(sent_dat.key,digitsPattern);
            word_id=num2cell(cellfun(@(x) str2num(x),word_id));
            % combine everything int stim_reps
            stim_resp_table=[cell2table([sent_id,word_id],"VariableNames",{'sentence_id','word_id'}),sent_dat];
            assert(size(stim_resp_table,1)==2113)
            % add stim_resp_table to the object 
            if not(isprop(obj,'stim_resp_table'))
                P = addprop(obj,'stim_resp_table');
            end 
            obj.stim_resp_table=stim_resp_table;
            % save a structure version to use in python 
            save_struct=struct;
            save_struct.s_vs_n_struct=table2struct(obj.s_vs_n_sig);
            % add electrode names 
            save_struct.s_vs_n_struct.elec_ch_label=obj.elec_ch_label;
            save_struct.s_vs_n_struct.bip_elec_ch_label=...
                arrayfun(@(X) sprintf('%s-%s',obj.bip_ch_label_valid{X,1},obj.bip_ch_label_valid{X,2}),1:size(obj.bip_ch_label_valid,1),'uni',false);
            save_struct.stim_resp_struct=table2struct(obj.stim_resp_table);
            save_file=sprintf('%s/%s_%s_%s_stim_resp_struct',obj.save_path,obj.subject_id,obj.experiment,obj.modality)
            save(save_file,'save_struct');
            % run python 
            %python_path='/Users/eghbalhosseini/anaconda3/envs/neural_nlp_1/bin/python';
            python_path='/Users/eghbalhosseini/anaconda3/envs/fmri_DNN/bin/python';
            commandStr = sprintf('%s /Users/eghbalhosseini/MyCodes/ecog_DNN/ecog_data_analysis/construct_stim_resp_data.py %s %s_%s',python_path,obj.subject_id,obj.experiment,obj.modality);
            [status, commandOut] = system(commandStr);
            if status==0
                fprintf('python conversion was successful\n');
                fprintf(commandOut)
            else 
                fprintf(commandOut)
            end
            
        end
        function obj=make_extended_stim_resp_mat(obj,varargin)
            p=inputParser();
            addParameter(p, 'words', 'all');
            parse(p, varargin{:});
            ops = p.Results;
            %get sentence id from events 
            aud_name=obj.events_table.final_audio_filename;
            list_id=obj.events_table.final_list;
            trial_id=obj.events_table.trial;
            trial_onset=obj.events_table.trial_onset;
            
            aud_name=erase(aud_name,'.wav');
            aud_transcript=obj.events_table.final_audio_transcript;
            % 
            cond_key='trial_start_end';
            trial_strings=cellfun(@(x) x.string{ismember(x.key,cond_key)},obj.trial_data,'uni',false);
            
            cond_key='word_1';
            trial_abs_onset=[];
            for temp=1:length(obj.trial_timing_dec)
                x=obj.trial_timing_dec{temp};
                x.start(ismember(x.key,cond_key));
                trial_abs_onset(temp,1)=x.start(ismember(x.key,cond_key));
            end 
            assert(all(trial_abs_onset==sort(trial_abs_onset)))
            
            B=obj.get_average(obj.trial_data);
            word_data=obj.get_value(B,'key','word','type','contain');
            word_dat_mod=arrayfun(@(x) [cell2table(repmat(aud_name(x),size(word_data{x},1),1),"VariableNames",{'stim_name'}),word_data{x}], 1:length(word_data),'uni',false)';
            word_dat_mod=arrayfun(@(x) [cell2table(repmat({x},size(word_dat_mod{x},1),1),"VariableNames",{'Trial_abs_id'}),word_dat_mod{x}], 1:length(word_dat_mod),'uni',false)';
            word_dat_mod=arrayfun(@(x) [cell2table(repmat({trial_abs_onset(x)},size(word_dat_mod{x},1),1),"VariableNames",{'Trial_abs_onset'}),word_dat_mod{x}], 1:length(word_dat_mod),'uni',false)';
            word_dat_mod=arrayfun(@(x) [cell2table(repmat({trial_id(x)},size(word_dat_mod{x},1),1),"VariableNames",{'Trial_id'}),word_dat_mod{x}], 1:length(word_dat_mod),'uni',false)';
            word_dat_mod=arrayfun(@(x) [cell2table(repmat({trial_onset(x)},size(word_dat_mod{x},1),1),"VariableNames",{'Trial_onset'}),word_dat_mod{x}], 1:length(word_dat_mod),'uni',false)';
            word_dat_mod=arrayfun(@(x) [cell2table(repmat({list_id(x)},size(word_dat_mod{x},1),1),"VariableNames",{'list_id'}),word_dat_mod{x}], 1:length(word_dat_mod),'uni',false)';
            all_dat=vertcat(word_dat_mod{:});
            stim_id=extract(all_dat.stim_name,digitsPattern);
            stim_id=num2cell(cellfun(@(x) str2num(x),stim_id));
            word_id=extract(all_dat.key,digitsPattern);
            word_id=num2cell(cellfun(@(x) str2num(x),word_id));
            stim_types={'S','N'};
            stim_type=cellfun(@(x) stim_types{isempty(regexp(x,'sentence'))+1},all_dat.stim_name,'uni',false);
            stim_value=cellfun(@(x) aud_transcript{find(ismember(aud_name,x),1,'first')},all_dat.stim_name,'uni',false);
            
            stim_resp_table=[cell2table([stim_id,word_id,stim_type,stim_value],"VariableNames",{'stim_id','word_id','stim_type','stim_value'}),all_dat];
            
            
            % add stim_resp_table to the object 
            if not(isprop(obj,'extended_stim_resp_table'))
                P = addprop(obj,'extended_stim_resp_table');
            end 
            obj.extended_stim_resp_table=stim_resp_table;
            % save a structure version to use in python 
            save_struct=struct;
            save_struct.s_vs_n_struct=table2struct(obj.s_vs_n_sig);
            save_struct.s_vs_n_ratio_struct=table2struct(obj.s_vs_n_p_ratio);
            % add electrode names 
            assert(all(find(obj.elec_ch_valid==1)==obj.elec_ch_clean))
            save_struct.s_vs_n_struct.elec_ch_label=obj.elec_ch_label;
            save_struct.s_vs_n_struct.elec_ch_valid=~isnan(obj.elec_ch_valid);
            save_struct.s_vs_n_struct.elec_ch_clean=obj.elec_ch_clean;
            save_struct.s_vs_n_struct.bip_elec_ch_valid=ones(length(obj.bip_ch_label_valid),1);
            save_struct.s_vs_n_struct.elec_pos_mni=obj.elec_ch_pos_mni;
            save_struct.s_vs_n_struct.bip_elec_pos_mni=obj.bip_ch_pos_mni;
            
            save_struct.s_vs_n_struct.elec_HCP_label=obj.elec_ch_HCP_label;
            save_struct.s_vs_n_struct.bip_HCP_label=obj.bip_ch_HCP_label;
            
            
            save_struct.s_vs_n_struct.bip_elec_ch_label=...
                arrayfun(@(X) sprintf('%s-%s',obj.bip_ch_label_valid{X,1},obj.bip_ch_label_valid{X,2}),1:size(obj.bip_ch_label_valid,1),'uni',false);
            save_struct.stim_resp_struct=table2struct(obj.extended_stim_resp_table);
            save_file=sprintf('%s/%s_%s_%s_extended_stim_resp_struct',obj.save_path,obj.subject_id,obj.experiment,obj.modality)
            save(save_file,'save_struct');
            % run python 
            %python_path='/Users/eghbalhosseini/anaconda3/envs/neural_nlp_1/bin/python';
            python_path='/Users/eghbalhosseini/anaconda3/envs/fmri_DNN/bin/python';
            commandStr = sprintf('%s /Users/eghbalhosseini/MyCodes/ecog_DNN/ecog_data_analysis/construct_stim_resp_data.py %s %s_%s',python_path,obj.subject_id,obj.experiment,[obj.modality]);
            [status, commandOut] = system(commandStr);
            if status==0
                fprintf('python conversion was successful\n');
                fprintf(commandOut)
            else 
                fprintf(commandOut)
            end
            
        end
%% method for doing test retest reliability 
        function obj=perform_test_retest_reliability(obj,varargin)
            p=inputParser();
            addParameter(p, 'aux_data', {});
            addParameter(p, 'words', 'all');
            parse(p, varargin{:});
            ops = p.Results;
            cond_key='trial_start_end';
            condtion_flag='S';
            cond_resp=obj.get_cond_resp(condtion_flag);
            % find repetition locations 
            trial_strings=cellfun(@(x) x.string{ismember(x.key,cond_key)},cond_resp,'uni',false);
            [D,~,~] = unique(trial_strings,'stable');
            repeat_loc=cellfun(@(y) find(ismember(trial_strings,y)),D,'uni',false);
            repeat_loc=repeat_loc(cellfun(@(x) length(x)~=1,repeat_loc));
            % check if strings are the same 
            cellfun(@(x) assert(strcmp(trial_strings{x(1)},trial_strings{x(2)})),repeat_loc);
            %repeat_resp=cellfun(@(x) cond_resp(x),repeat_loc,'uni',false);
            repeat_resp=arrayfun(@(x) cond_resp{x},cell2mat(repeat_loc')','uni',false);
            word_repeat=cellfun(@(x) obj.get_value(x,'key','word','type','contain'),repeat_resp,'uni',false);
            % combine responses along their second dimension 
            values=cellfun(@(x) x(:,~ismember(x.Properties.VariableNames,{'key','string'})), word_repeat,'uni',false);
            keys=cellfun(@(x) x(:,ismember(x.Properties.VariableNames,{'key'})), word_repeat,'uni',false);
            strings=cellfun(@(x) x(:,ismember(x.Properties.VariableNames,{'string'})), word_repeat,'uni',false);
            values_stacked=arrayfun(@(z_id) cellfun(@(x,y) cat(3,x,y),values{z_id,1}.Variables,values{z_id,2}.Variables,'uni',false),[1:size(values,1)]','uni',false);
            variable_names=cellflat({'key','string',setdiff(repeat_resp{1}.Properties.VariableNames, {'key','string'},'stable')});
            repeat_table=cellfun(@(x,y,z) cell2table(horzcat(x.key,y.string,z),'VariableNames',variable_names), keys(:,1),strings(:,1),values_stacked,'uni',false);
            if not(isprop(obj,'stim_repeat_table'))
                P = addprop(obj,'stim_repeat_table');
            end 
            obj.stim_repeat_table=repeat_table;
            % create one for full length of the sentence 
            val_sent_stack={};
            for pp=1:numel(values_stacked)
                val=values_stacked{pp};
                % stack along 2nd dimension 
                val_sent=arrayfun(@(col_id) cat(2,val{:,col_id}),1:size(val,2),'uni',false);
                val_sent_stack=[val_sent_stack;val_sent];
            end 
            
            val_sent_all=arrayfun(@(col_id) cat(2,val_sent_stack{:,col_id}),1:size(val_sent_stack,2),'uni',false);
            val_sent_per_elec=cellfun(@(X) mat2cell(X,ones(1,size(X,1))), val_sent_all,'uni',false) ;
            val_sent_corr=cellfun(@(X) cellfun(@(Y) corrcoef(squeeze(Y),'rows','pairwise'),X,'uni',false), val_sent_per_elec, 'uni',false);
            val_sent_corr=cellfun(@(X) cellfun(@(Y) Y(1,2),X,'uni',false), val_sent_corr, 'uni',false);
            % get averages 
            repeat_ave=obj.get_average(obj.stim_repeat_table);
            values=cellfun(@(x) x(:,~ismember(x.Properties.VariableNames,{'key','string'})),repeat_ave,'uni',false);
            % breakdown the cell into individual electrodes for computing
            % the correlation
            func=@(x) cell2mat(reshape(vertcat(x),1,[]));
            values_comb=cellfun(@(Z) cellfun(@(X) cat(3,Z.(X){:}),Z.Properties.VariableNames,'uni',false),values,'uni',false);
            values_comb=cat(1,values_comb{:});
            values_comb=arrayfun(@(col_id) cat(3,values_comb{:,col_id}),1:size(values_comb,2),'uni',false);
            V=cellfun(@(X) mat2cell(X,ones(1,size(X,1))), values_comb,'uni',false) ;
            %V=cat(1,V{:});
            V_corr=cellfun(@(X) cellfun(@(Y) corr(squeeze(Y)'),X,'uni',false), V, 'uni',false);
            V_corr=cellfun(@(X) cellfun(@(Y) Y(1,2),X,'uni',false), V_corr, 'uni',false);
            V_corr=cellfun(@cell2mat,V_corr,'uni',false);
            %V_corr=arrayfun(@(x) horzcat(V_corr{:,x}), 1:size(V_corr,2),'uni',false);
            names_no_string=~ismember(variable_names,{'string'});
            test_retest_corr=cell2table(horzcat('test_retest_reliablity',V_corr),'VariableNames',variable_names(names_no_string));
            if not(isprop(obj,'test_retest_corr'))
                P = addprop(obj,'test_retest_corr');
            end
            obj.test_retest_corr=test_retest_corr;
            % do some plotting for the data 
            close all
            num_rows=5;
            num_columns=2;  
            nbins=50;
            total_plots=num_rows*num_columns;
            pp=0;
            f = figure;
            set(f,'position',[1123          29        1266        1275]);
            value_names=setdiff(test_retest_corr.Properties.VariableNames,{'key'},'stable');
            value_names_pr=strrep(value_names,'_', ' ');
            analysis_path=strcat(obj.save_path,'analysis/test_retest_reliability/');
            if ~exist(strcat(analysis_path,obj.subject_id))
                            mkdir(strcat(analysis_path,obj.subject_id));
            end
            pp=0;
            var_id=4;
            %
            condtion_flag='S';
            [S_ave_tbl]=obj.get_ave_cond_resp('condition',condtion_flag);
            condtion_flag='N';
            [N_ave_tbl]=obj.get_ave_cond_resp('condition',condtion_flag);
            S_dat=S_ave_tbl.(value_names{var_id});
            N_dat=N_ave_tbl.(value_names{var_id});
            
            if not(isempty(ops.aux_data))
                condtion_flag='S';
                [~,S_table_sn]=sn_data.get_ave_cond_trial('condition',condtion_flag);
                condtion_flag='N';
                [~,N_table_sn]=sn_data.get_ave_cond_trial('condition',condtion_flag);
                S_sn_dat=S_table_sn.(value_names{var_id}){1};
                N_sn_dat=N_table_sn.(value_names{var_id}){1};
            end 

            for kk=1:numel(obj.elec_ch_clean)
                elec_id=obj.elec_ch_clean(kk);
                elec_dat=squeeze(V{:,var_id}{elec_id});
                %elec_dat=cellfun(@(x) squeeze(x{elec_id}), V(:,var_id),'uni',false);
                % plot test retest 
                ax=subplot(num_rows,num_columns,num_columns*(kk-num_rows*fix((kk-1)/num_rows)-1)+1);
                hold on
                colormap=inferno(numel(elec_dat));
                scatter(elec_dat(1,:),elec_dat(2,:),15,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0],'LineWidth',1.5)
                limits=[min([ax.XLim,ax.YLim]),max([ax.XLim,ax.YLim])];
                ax.XLim=limits;
                ax.YLim=limits;
                daspect([1,1,1])
                plot(limits,limits,'k--');

                [R,P]=corrcoef(elec_dat');
                ax.Title.String=sprintf('%s ,%s , \n corr=%.3f p_{ratio}=%0.4f',obj.elec_ch_label{elec_id},value_names_pr{var_id},R(1,2),P(1,2));
                ax.Title.FontSize=8;
               
                % plot s n difference 
                ax=subplot(num_rows,num_columns,num_columns*(kk-num_rows*fix((kk-1)/num_rows)-1)+2);
                % sentences
                s_electrode_resp=cell2mat(cellfun(@(x) x(elec_id,:), S_dat,'uni',false));
                word_pos=repmat(1:size(s_electrode_resp,1),size(s_electrode_resp,2),1)';
                y=s_electrode_resp;
                x=word_pos;
                b1=plot(mean(x,2)+.1,nanmean(y,2),'color',[1,.5,.5],'linewidth',1,'marker','o','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,0],'displayname','S');
                b1.MarkerSize=2;
                hold on
                bl=errorbar(mean(x,2)+.1,nanmean(y,2),nanstd(y,[],2)./sqrt(size(y,2)));   

                bl.LineStyle='none';
                bl.Color=[1,.5,.5];
                bl.LineWidth=1;
                bl.CapSize=1;
                hAnnotation = get(bl,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
                
                ax.XTick=1:max(word_pos(:));
                ax.XLim=[0,max(word_pos(:))+1];
                % nonwords
                n_electrode_resp=cell2mat(cellfun(@(x) x(elec_id,:), N_dat,'uni',false));
                word_pos=repmat(1:size(n_electrode_resp,1),size(n_electrode_resp,2),1)';
                y=n_electrode_resp;
                x=word_pos;
                b1=plot(mean(x,2)+.1,nanmean(y,2),'color',[.5,.5,1],'linewidth',1,'marker','o','MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,1],'displayname','N');
                b1.MarkerSize=2;
                hold on
                bl=errorbar(mean(x,2)+.1,nanmean(y,2),nanstd(y,[],2)./sqrt(size(y,2)));
                bl.LineStyle='none';
                bl.Color=[.5,.5,1];
                bl.LineWidth=1;
                bl.CapSize=1;
                hAnnotation = get(bl,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
                ax.XAxis.Visible = 'on';
                h=get(ax,'children');
                ax.FontSize=6;
                set(ax,'ydir', 'normal','box','off');
                ah =get(ax,'children');
                arrayfun(@(x) set(ah(x),'DisplayName',''),[1:2]);
                ah(2).DisplayName='N';
                ah(4).DisplayName='S';
                set(ax,'children',ah);
                legend(ah([2,4]),'Location','northwest','NumColumns',2)
                xlabel('word position');
                ax.YLabel.String='High Gamma (a.u.)';
                ax.XAxis.LineWidth=2;
                ax.YAxis.LineWidth=2;
                ax.Title.String='S vs. N in DNN experiment'
                
                
                % plot s n difference for sn experiment 
%                 ax=subplot(num_rows,num_columns,num_columns*(kk-num_rows*fix((kk-1)/num_rows)-1)+3);
%                 % sentences
%                 s_electrode_resp=squeeze(S_sn_dat(elec_id,:,:))';
%                 n_electrode_resp=squeeze(N_sn_dat(elec_id,:,:))';
%                 word_pos=repmat(1:size(s_electrode_resp,1),size(s_electrode_resp,2),1)';
%                 y=s_electrode_resp;
%                 x=word_pos;
%                 b1=plot(mean(x,2)+.1,nanmean(y,2),'color',[1,.5,.5],'linewidth',1,'marker','o','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,0],'displayname','S');
%                 b1.MarkerSize=2;
%                 hold on
%                 bl=errorbar(mean(x,2)+.1,nanmean(y,2),nanstd(y,[],2)./sqrt(size(y,2)));   
%                 bl.LineStyle='none';
%                 bl.Color=[1,.5,.5];
%                 bl.LineWidth=1;
%                 bl.CapSize=1;
%                 hAnnotation = get(bl,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
%                 
%                 ax.XTick=1:max(word_pos(:));
%                 ax.XLim=[0,max(word_pos(:))+1];
%                 % nonwords
%                 word_pos=repmat(1:size(n_electrode_resp,1),size(n_electrode_resp,2),1)';
%                 y=n_electrode_resp;
%                 x=word_pos;
%                 b1=plot(mean(x,2)+.1,nanmean(y,2),'color',[.5,.5,1],'linewidth',1,'marker','o','MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,1],'displayname','N');
%                 hold on
%                 b1.MarkerSize=2;
%                 bl=errorbar(mean(x,2)+.1,nanmean(y,2),nanstd(y,[],2)./sqrt(size(y,2)));
%                 bl.LineStyle='none';
%                 bl.Color=[.5,.5,1];
%                 bl.LineWidth=1;
%                 bl.CapSize=1;
%                 hAnnotation = get(bl,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
%                 ax.XAxis.Visible = 'on';
%                 h=get(ax,'children');
%                 ax.FontSize=6;
%                 set(ax,'ydir', 'normal','box','off');
%                 ah =get(ax,'children');
%                 arrayfun(@(x) set(ah(x),'DisplayName',''),[1:2]);
%                 ah(2).DisplayName='N';
%                 ah(4).DisplayName='S';
%                 set(ax,'children',ah);
%                 legend(ah([2,4]),'Location','northwest','NumColumns',2)
%                 xlabel('word position');
%                 ax.YLabel.String='High Gamma (a.u.)';
%                 ax.XAxis.LineWidth=2;
%                 ax.YAxis.LineWidth=2;
%                 ax.Title.String='S vs. N in SN experiment';
                
                    if ~mod(kk,num_rows) | kk==numel(value_names)
                        %legend('show','Location','northeastoutside')
                        pp=pp+1;
                        set(gcf,'PaperPosition',[.25 .25 11 8]);
                        set(gcf,'PaperOrientation','portrait');
                        fname=sprintf('%s_test_retests_pg_%d.pdf',obj.subject_id,pp);
                        print(f, '-bestfit','-dpdf','-opengl', strcat(analysis_path,obj.subject_id,'/',fname));
                        %close(f)
                        set(0, 'CurrentFigure', f);
                        clf reset;
                    end

                 
            end
            close(f);
            
        end 
        %% methods for looking at trial level data 
        function [output_tbl]=get_ave_cond_resp(obj,varargin)
            % input definitions 
            p=inputParser();
            addParameter(p, 'condition', 'S');
            parse(p, varargin{:});
            ops = p.Results;
            % output definitions
            output_tbl=table();
            
            cond_data=obj.get_cond_resp(ops.condition);
            cond_data_ave=obj.get_average(cond_data);
            word_data=obj.get_value(cond_data_ave,'key','word','type','contain');
            % up to here it is the same as ecog_sn however the trials have
            % different sizes so the following steps are different 
            func=@(x) cell2mat(reshape(vertcat(x),1,[]));
            all_keys=cellfun(@(x) x.key, word_data,'uni',false);
            [~,longest_string]=max(cellfun(@(x) length(x), all_keys));
            
            for k=1:numel(all_keys{longest_string})
                each_key=all_keys{longest_string}{k};
                temp=(cellfun(@(x) x(ismember(x.key,each_key),:),word_data,'uni',false));
                % fill empty tables with nans 
                nonempty_cells=cellfun(@(x) ~isempty(x), temp); 
                if any(ismember(nonempty_cells,0))
                    template=temp{find(nonempty_cells,1)};
                    values=template(:,~ismember(template.Properties.VariableNames,{'key','string'}));
                    values_comb=cellfun(@(X) cellfun(@(y) y*nan,values.(X),'uni',false),values.Properties.VariableNames,'uni',false);
                    template_table=cell2table(horzcat(each_key,{'N/A'},values_comb),'VariableNames',cond_tbl.Properties.VariableNames);
                    for idx=find(~nonempty_cells)'
                        temp{idx}=template_table;
                    end
                end
                cond_tbl=vertcat(temp{:});
                strings=cond_tbl(:,ismember(cond_tbl.Properties.VariableNames,'string'));
                values=cond_tbl(:,~ismember(cond_tbl.Properties.VariableNames,{'key','string'}));
                
                values_comb=cellfun(@(X) func(values.(X)),values.Properties.VariableNames,'uni',false);
                temp_table=cell2table(horzcat(each_key,{strings.string},values_comb),'VariableNames',cond_tbl.Properties.VariableNames);
                output_tbl=[output_tbl;temp_table];
            end 
            
            
        end 
        function plot_elec_on_mni(obj,varargin)
            defaultHemi = 'lh';
            expectedHemi = {'lh','rh','both'};
            expectedSignal = {'elec_data_dec','elec_data_zs_dec','envelope_dec','bip_elec_data_dec','bip_elec_data_zs_dec','bip_envelope_dec'};
            defaultSignal = 'elec_data_zs_dec';
            p=inputParser();
            addParameter(p, 'mni_path', '/Applications/freesurfer/7.2.0/subjects/cvs_avg35_inMNI152/surf/');
            addParameter(p, 'mode', 'pial');
            addParameter(p,'hemi',defaultHemi,...
                @(x) any(validatestring(x,expectedHemi)));
            addParameter(p, 'tag_lang_elec', true);
            addParameter(p,'sig_type',defaultSignal,...
                @(x) any(validatestring(x,expectedSignal)));
            parse(p, varargin{:});
            ops = p.Results;
            analysis_path=strcat(obj.save_path,'analysis/elec_in_MNI/',obj.subject_id);
            if ~exist(analysis_path);mkdir(analysis_path);end
            %%
            d_anat_lh= fullfile(strcat(ops.mni_path,'/','lh.pial'));
            if isfile(d_anat_lh),
                [lh_pial.vert, lh_pial.face] = freesurfer_read_surf(d_anat_lh);
                nvertices=size(lh_pial.vert,1);
                lh_color_data = repmat([255,229,204]./255,nvertices,1);
            end
            d_anat_rh= fullfile(strcat(ops.mni_path,'/','rh.pial'));
            if isfile(d_anat_rh),
                [rh_pial.vert, rh_pial.face] = freesurfer_read_surf(d_anat_rh);
                nvertices=size(rh_pial.vert,1);
                rh_color_data = repmat([255,229,204]./255,nvertices,1);
            end
            %
            X=obj.filt_ops.elec_clinic_info.mni_linear_x;
            Y=obj.filt_ops.elec_clinic_info.mni_linear_y;
            Z=obj.filt_ops.elec_clinic_info.mni_linear_z;
            lang_elec=find(obj.s_vs_n_sig.(ops.sig_type){1});
            %
            lh_cfg=struct;
            lh_cfg.xlim=[-70 5];lh_cfg.ylim=[-90 90];lh_cfg.zlim=[-70,70];
            lh_cfg.saggital_lateral={[-90 0],'headlight'};
            lh_cfg.saggital_medial={[90 0],'headlight'};
            lh_cfg.coronal_superior={[90 90],'headlight'};
            lh_cfg.coronal_inferior={[-90 -90],'headlight'};
            lh_cfg.horizontal_anterior={[180 0],'left'};
            lh_cfg.horizontal_posterior={[0 0],'headlight'};
            %
            rh_cfg=struct;
            rh_cfg.xlim=[-5 75];rh_cfg.ylim=[-90 90];rh_cfg.zlim=[-70,70];
            rh_cfg.saggital_lateral={[90 0],'headlight'};
            rh_cfg.saggital_medial={[-90 0],'headlight'};
            rh_cfg.coronal_superior={[90 90],'headlight'};
            rh_cfg.coronal_inferior={[-90 -90],'headlight'};
            rh_cfg.horizontal_anterior={[180 0],'right'};
            rh_cfg.horizontal_posterior={[0 0],'headlight'};
            %
            fig_cfg.pap_ratio=8.5/11;
            axis_cfg.saggital_lateral=[.05,.7,.4,.4*fig_cfg.pap_ratio];
            axis_cfg.saggital_medial=[.5,.7,.4,.4*fig_cfg.pap_ratio];
            axis_cfg.coronal_superior=[.05,.35,.4,.4*fig_cfg.pap_ratio];
            axis_cfg.coronal_inferior=[.5,.35,.4,.4*fig_cfg.pap_ratio];
            axis_cfg.horizontal_anterior=[.05,.05,.4,.4*fig_cfg.pap_ratio];
            axis_cfg.horizontal_posterior=[.5,.05,.4,.4*fig_cfg.pap_ratio];
            %
            orients={'saggital_lateral','saggital_medial','coronal_superior','coronal_inferior','horizontal_anterior','horizontal_posterior'};
            if (strcmp(ops.hemi,'lh') | strcmp(ops.hemi,'both'))
                f=figure(1);
                clf(f)
                f.Units='Inches';
                f.PaperOrientation='portrait';
                f.Position=[9. 1.   8.5   11];
                f.Color=[1,1,1];
                lh_elec=find(X<max(lh_pial.vert(:,1)));
                for k=1:numel(orients)
                    orient=orients{k};
                    ax=axes('position',axis_cfg.(orient));
                    patch_handle = patch('vertices', lh_pial.vert, 'Faces', lh_pial.face, 'FaceVertexCData', lh_color_data,'FaceLighting','gouraud','SpecularStrength',0,'DiffuseStrength',0.7);
                    shading interp;view(lh_cfg.(orient){1});camlight(lh_cfg.(orient){2},'infinite')
                    ax.XLim=lh_cfg.xlim;ax.YLim=lh_cfg.ylim;ax.ZLim=lh_cfg.zlim;
                    ax.XAxis.Visible='off';ax.YAxis.Visible='off';ax.ZAxis.Visible='off';
                    daspect([1,1,1]);patch_handle.FaceAlpha=.4;
                    hold on
                    if ~isempty(lh_elec)
                    h=plot3(X(lh_elec),Y(lh_elec),Z(lh_elec));
                    h.Marker='o';h.MarkerFaceColor=[0,0,0];h.MarkerEdgeColor=[0,0,0];h.MarkerSize=3;h.LineStyle='None';
                    if (ops.tag_lang_elec & ~isempty(intersect(lh_elec,lang_elec)))
                        lh_lang=intersect(lh_elec,lang_elec);
                        h=plot3(X(lh_lang),Y(lh_lang),Z(lh_lang));
                        h.Marker='o';h.MarkerFaceColor=[1,0,0];h.MarkerEdgeColor=[1,0,0];h.MarkerSize=5;h.LineStyle='None';
                    end
                    end 
                    ax.Title.String=sprintf('lh - %s', replace(orient,'_',' '));
                end
                print(f, '-bestfit','-dpdf','-opengl', sprintf('%s/%s_%s_lh.pdf',analysis_path,obj.subject_id,ops.sig_type));
                
            end
            if (strcmp(ops.hemi,'rh') | strcmp(ops.hemi,'both'))
                f=figure(1);
                clf(f)
                f.Units='Inches';
                f.PaperOrientation='portrait';
                f.Position=[9. 1.   8.5   11];
                f.Color=[1,1,1];
                rh_elec=find(X>min(rh_pial.vert(:,1)));
                for k=1:numel(orients)
                    orient=orients{k};
                    ax=axes('position',axis_cfg.(orient));
                    patch_handle = patch('vertices', rh_pial.vert, 'Faces', rh_pial.face, 'FaceVertexCData', rh_color_data,'FaceLighting','gouraud','SpecularStrength',0,'DiffuseStrength',0.7);
                    shading interp;view(rh_cfg.(orient){1});camlight(rh_cfg.(orient){2},'infinite')
                    ax.XLim=rh_cfg.xlim;ax.YLim=rh_cfg.ylim;ax.ZLim=rh_cfg.zlim;
                    ax.XAxis.Visible='off';ax.YAxis.Visible='off';ax.ZAxis.Visible='off';
                    daspect([1,1,1]);patch_handle.FaceAlpha=.4;
                    hold on
                    if ~isempty(rh_elec)
                    h=plot3(X(rh_elec),Y(rh_elec),Z(rh_elec));
                    h.Marker='o';h.MarkerFaceColor=[0,0,0];h.MarkerEdgeColor=[0,0,0];h.MarkerSize=3;h.LineStyle='None';
                    if (ops.tag_lang_elec & ~isempty(intersect(rh_elec,lang_elec)))
                        rh_lang=intersect(rh_elec,lang_elec);
                        h=plot3(X(rh_lang),Y(rh_lang),Z(rh_lang));
                        h.Marker='o';h.MarkerFaceColor=[1,0,0];h.MarkerEdgeColor=[1,0,0];h.MarkerSize=5;h.LineStyle='None';
                    end 
                    end 
                    ax.Title.String=sprintf('rh - %s', replace(orient,'_',' '));
                end
                print(f, '-bestfit','-dpdf','-opengl', sprintf('%s/%s_%s_rh.pdf',analysis_path,obj.subject_id,ops.sig_type));
            end
            
        end 
    end
end

