
classdef ecog_sn_data<ecog_data
    %ECOG_DATA Summary of this class goes here
    %   Detailed explanation goes here
    properties
    end
    
    methods
%% method for constructing the object 
function sn_obj = ecog_sn_data(prep_data,...
                               trial_timing_dec,...
                               events_table) 
    sn_obj@ecog_data(prep_data,...
                     trial_timing_dec,...
                     events_table)
end

%% method for combining trials 
        function output_d=combine_trial_cond(obj,input_d)
            all_keys=cellfun(@(x) x.key,input_d,'uni',false);
            % make sure the keys are the same for all trials 
            [X,Y] = ndgrid(1:numel(all_keys));
            Z = tril(true(numel(all_keys)),-1);
            assert(all(arrayfun(@(x,y)isequal(all_keys{x},all_keys{y}),X(Z),Y(Z))));
            func=@(x) cell2mat(reshape(vertcat(x),1,[]));
            output_d=table();
            for k=1:numel(all_keys{1})
                each_key=all_keys{1}{k};
                temp=(cellfun(@(x) x(ismember(x.key,each_key),:),input_d,'uni',false));
                cond_tbl=vertcat(temp{:});
                strings=cond_tbl(:,ismember(cond_tbl.Properties.VariableNames,'string'));
                values=cond_tbl(:,~ismember(cond_tbl.Properties.VariableNames,{'key','string'}));
                
                values_comb=cellfun(@(X) func(values.(X)),values.Properties.VariableNames,'uni',false);
                temp_table=cell2table(horzcat(each_key,{strings.string},values_comb),'VariableNames',cond_tbl.Properties.VariableNames);
                output_d=[output_d;temp_table];
            end 
            
        end
%% methods for calculating average signal
        function [output_tbl,cond_table]=get_ave_cond_trial(obj,varargin)
            p=inputParser();
            addParameter(p, 'words', 1:12);
            addParameter(p, 'condition', 'S');
            parse(p, varargin{:});
            ops = p.Results;
            
            func=@(x) cell2mat(permute(x,[3,2,1])); % format : electrdoe*trial_id*state
            % get trials condition
            condtion_flag=ops.condition;
            cond_data=obj.get_cond_resp(condtion_flag);
            cond_data_ave=obj.get_average(cond_data);
            word_data=obj.get_value(cond_data_ave,'key','word','type','contain');
            B=obj.combine_trial_cond(word_data);
            keys=B(:,ismember(B.Properties.VariableNames,'key'));
            strings=B(:,ismember(B.Properties.VariableNames,'string'));
            values=B(:,~ismember(B.Properties.VariableNames,{'key','string'}));
            values_comb=cellfun(@(X) func(values.(X)),values.Properties.VariableNames,'uni',false);
            cond_table=cell2table(horzcat(condtion_flag,{strings.string},values_comb),'VariableNames',B.Properties.VariableNames);
            % 
            %% select how many  words are selected : default 1:12 and average over words 
            %create the comparision between the two condition
            func_1=@(x) x{1}(:,:,ops.words); % format : electrdoe*trial_id
            func_2=@(x) nanmean(x,3); % format : electrdoe*trial_id
            B=cond_table;
            keys=B(:,ismember(B.Properties.VariableNames,'key'));
            strings=B(:,ismember(B.Properties.VariableNames,'string'));
            values=B(:,~ismember(B.Properties.VariableNames,{'key','string'}));
            condition_ave=cellfun(@(X) func_2(func_1(values.(X))),values.Properties.VariableNames,'uni',false);
            output_tbl=cell2table(horzcat(condtion_flag,strings.string,condition_ave),'VariableNames',B.Properties.VariableNames);
        end 

%% methods for finding langauge electrodes
        function obj=test_s_vs_n(obj,varargin)
            p=inputParser();
            addParameter(p, 'words', 1:12);
            addParameter(p, 'n_rep', 1000);
            addParameter(p,'corr_type','Spearman');
            addParameter(p, 'threshold', 0.01);
            addParameter(p, 'side', 'both');
            addParameter(p, 'do_plot',false);
            parse(p, varargin{:});
            ops = p.Results;
            obj.s_vs_n_ops=ops;
            %%
            % S condition
            condtion_flag='S';
            [S_ave_tbl,S_table]=obj.get_ave_cond_trial('words',ops.words,'condition',condtion_flag);
            B=S_ave_tbl;
            %keys=B(:,ismember(B.Properties.VariableNames,'key'));
            S_ave=table2cell(B(:,~ismember(B.Properties.VariableNames,{'key','string'})));
            % N condition
            condtion_flag='N';
            [N_ave_tbl,N_table]=obj.get_ave_cond_trial('words',ops.words,'condition',condtion_flag);
            B=N_ave_tbl;
            keys=B(:,ismember(B.Properties.VariableNames,'key'));
            N_ave=table2cell(B(:,~ismember(B.Properties.VariableNames,{'key','string'})));
            %% get correlation
            S_N_comb=cellfun(@(x,y) horzcat(x,y), S_ave,N_ave,'uni',false);
            S_N_flag=cellfun(@(x,y) horzcat(x*0+1,y*0-1), S_ave,N_ave,'uni',false);
            S_N_corr=cellfun(@(x,y) diag(corr(x',y','Type',ops.corr_type)),S_N_comb,S_N_flag,'uni',false);
            names_no_string=~ismember(B.Properties.VariableNames,{'string'});
            S_N_corr_table=cell2table(horzcat('s_vs_n_corr',S_N_corr),'VariableNames',B.Properties.VariableNames(names_no_string));
            %
            % do the permuation
            n_rep=ops.n_rep;
            fprintf('permutation :');
            S_N_corr_rnd_all=cell(size(S_N_corr));
            pbar=ProgressBar(n_rep, ...
            'Title', 'Permutations');
            for k=1:n_rep
                
                S_N_flag_rnd_idx=cellfun(@(x) randperm(size(x,2)),S_N_flag,'uni',false);
                S_N_flag_rnd=cellfun(@(x,y) x(:,y),S_N_flag,S_N_flag_rnd_idx,'uni',false);
                S_N_corr_rnd=cellfun(@(x,y) diag(corr(x',y','Type','Spearman')),S_N_comb,S_N_flag_rnd,'uni',false);
                S_N_corr_rnd_all=cellfun(@(x,y) horzcat(x,y),S_N_corr_rnd_all,S_N_corr_rnd,'UniformOutput',false);
                pbar(1,[],[]);
            end
            pbar.release();
            S_N_corr_rnd_table=cell2table(horzcat('s_vs_n_corr_rnd',S_N_corr_rnd_all),'VariableNames',B.Properties.VariableNames(names_no_string));

            %% get the significanc
            switch ops.side
                case 'left'
                    S_N_p_ratio=cellfun(@(x,y) sum(x<y,2)/size(y,2),S_N_corr,S_N_corr_rnd_all,'uni',false);
                    S_N_p_is_sig=cellfun(@(x) (x<(ops.threshold)),S_N_p_ratio,'uni',false);
                case 'right'
                    S_N_p_ratio=cellfun(@(x,y) sum(x>y,2)/size(y,2),S_N_corr,S_N_corr_rnd_all,'uni',false);
                    S_N_p_is_sig=cellfun(@(x) (x>(1-ops.threshold) ),S_N_p_ratio,'uni',false);
                case 'both'
                    S_N_p_ratio=cellfun(@(x,y) sum(x>y,2)/size(y,2),S_N_corr,S_N_corr_rnd_all,'uni',false);
                    S_N_p_is_sig=cellfun(@(x)  ((x<(ops.threshold)) | (x>(1-ops.threshold) )),S_N_p_ratio,'uni',false);
            end 
            S_N_p_ratio_tbl=cell2table(horzcat('s_vs_n_p_ratio',S_N_p_ratio),'VariableNames',B.Properties.VariableNames(names_no_string));
            
            s_vs_n_sig=cell2table(horzcat('s_vs_n_sig',S_N_p_is_sig),'VariableNames',B.Properties.VariableNames(names_no_string));
            assert(all(find(obj.elec_ch_valid==1)==obj.elec_ch_clean))
            s_vs_n_sig.Properties.VariableNames
            obj.s_vs_n_sig=s_vs_n_sig;
            obj.s_vs_n_p_ratio=S_N_p_ratio_tbl;
            %% plotting
            if ops.do_plot
                analysis_path=strcat(obj.save_path,'analysis/s_n_difference/');
                elec_flag='bip_elec_data_dec';
                S_dat=S_table.(elec_flag){1};
                N_dat=N_table.(elec_flag){1};
                SN_corr_dat=S_N_corr_table.(elec_flag){1};
                SN_cor_rnd_dat=S_N_corr_rnd_table.(elec_flag){1};
                SN_sig_dat=s_vs_n_sig.(elec_flag){1};
                SN_p_ratio_dat=S_N_p_ratio_tbl.(elec_flag){1};
                ecog_colors=cbrewer('seq', 'OrRd', 1024,'linear');
                close all
                num_rows=3;
                num_columns=3;
                nbins=50;
                total_plots=num_rows*num_columns;
                pp=0;
                f = figure;
                set(f,'position',[1123          29        1266        1275])
                for i=1:size(S_dat,1)
                    s_electrode_resp=squeeze(S_dat(i,:,:))';
                    n_electrode_resp=squeeze(N_dat(i,:,:))';
                    s_n_rho=SN_corr_dat(i);
                    s_n_rho_rnd=SN_cor_rnd_dat(i,:);
                    is_sig=SN_sig_dat(i);
                    p_ratio=SN_p_ratio_dat(i);
                    word_pos=repmat(1:size(s_electrode_resp,1),size(s_electrode_resp,2),1)';
                    %f=figure(fix((i-1)/total_plots)+1);
                    
                    %set(f,'position',[1123          29        1266        1275])
                    sup_title=(strcat(obj.bip_ch_label_valid{i,1},'-',obj.bip_ch_label_valid{i,2}));
                    ax=subplot(num_rows,num_columns,3*(i-num_rows*fix((i-1)/num_rows)-1)+1);
                    
                    % plot the distbirtubiotn
                    h1=histogram(s_n_rho_rnd,nbins);
                    hold on
                    xline(s_n_rho,'linewidth',3);
                    ax.Box='off'
                    h1.EdgeColor='w';
                    shg
                    ax.XAxis.LineWidth=2;
                    ax.YAxis.LineWidth=2;
                    ax.Title.String=sprintf('corr=%f,\n p_{ratio}=%0.4f sig=%d',s_n_rho,p_ratio,is_sig);
                    %
                    ax=subplot(num_rows,num_columns,3*(i-num_rows*fix((i-1)/num_rows)-1)+2);
                    %x=.3;y=s_electrode_resp(:);
                    b1=plot(mean(word_pos,2)+.1,mean(s_electrode_resp,2),'color',[1,.5,.5],'linewidth',2,'marker','o','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,0],'displayname','S');
                    hold on
                    b1=plot(mean(word_pos,2)-.1,mean(n_electrode_resp,2),'color',[.5,.5,1],'linewidth',2,'marker','o','MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,1],'displayname','N');
                    y=s_electrode_resp;
                    x=word_pos;
                    bl=errorbar(mean(x,2)+.1,mean(y,2),std(y,[],2)./sqrt(size(y,2)));
                    bl.LineStyle='none';
                    bl.Color=[1,.5,.5];
                    bl.LineWidth=2;
                    bl.CapSize=2;
                    hAnnotation = get(bl,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
                    %
                    y=n_electrode_resp;
                    x=word_pos;
                    bl=errorbar(mean(x,2)-.1,mean(y,2),std(y,[],2)./sqrt(size(y,2)));
                    bl.LineStyle='none';
                    bl.Color=[.5,.5,1];
                    bl.LineWidth=2;
                    bl.CapSize=2;
                    hAnnotation = get(bl,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
                    ax.XAxis.Visible = 'on';
                    ax.XTick=1:max(word_pos(:));
                    ax.XLim=[0,max(word_pos(:))+1];
                    %ax.XTick=1:size(electrode_resp,1);
                    all_points=[s_electrode_resp(:);n_electrode_resp(:)];
                    y_quantile=quantile(all_points,10);
                    h=get(ax,'children');
                    ax.FontSize=12;
                    set(ax,'ydir', 'normal','box','off','ylim',[y_quantile(1),y_quantile(end)]);
                    ah =get(ax,'children');
                    arrayfun(@(x) set(ah(x),'DisplayName',''),[1:2]);
                    ah(3).DisplayName='N';
                    ah(4).DisplayName='S';
                    set(ax,'children',ah);
                    legend(ah(3:4),'Location','northwest','NumColumns',2)
                    xlabel('word position');
                    ax.YLabel.String='High Gamma (a.u.)';
                    ax.XAxis.LineWidth=2;
                    ax.YAxis.LineWidth=2;
                    ax.Title.String=erase(sup_title,'_');
                    if ~mod(i,num_rows) | i==size(S_dat,1)
                        %legend('show','Location','northeastoutside')
                        pp=pp+1;
                        if ~exist(strcat(analysis_path,obj.subject_id))
                            mkdir(strcat(analysis_path,obj.subject_id));
                        end
                        set(gcf,'PaperPosition',[.25 .25 8 6])
                        set(gcf,'PaperOrientation','landscape');
                        fname=sprintf('%s_s_v_n_words-%d-%d_p_%0.2f_%s_%s_%s_bipolar_1.pdf',obj.subject_id,...
                        min(ops.words),...
                        max(ops.words),...
                        ops.threshold,...
                        ops.side,...
                        obj.modality,...
                        num2str(num2str(pp)));
                        print(f, '-bestfit','-dpdf','-opengl', strcat(analysis_path,obj.subject_id,'/',fname));
                        %close(f)
                        set(0, 'CurrentFigure', f);
                        clf reset;
                    end
                    
                end
                close(f);
               
            end
                        %% plotting
            if ops.do_plot
                analysis_path=strcat(obj.save_path,'analysis/s_n_difference/');
                elec_flag='envelope_dec';
                S_dat=S_table.(elec_flag){1};
                N_dat=N_table.(elec_flag){1};
                SN_corr_dat=S_N_corr_table.(elec_flag){1};
                SN_cor_rnd_dat=S_N_corr_rnd_table.(elec_flag){1};
                SN_sig_dat=s_vs_n_sig.(elec_flag){1};
                SN_p_ratio_dat=S_N_p_ratio_tbl.(elec_flag){1};
                ecog_colors=cbrewer('seq', 'OrRd', 1024,'linear');
                close all
                num_rows=3;
                num_columns=3;
                nbins=50;
                total_plots=num_rows*num_columns;
                pp=0;
                f = figure;
                set(f,'position',[1123          29        1266        1275])
                for i=1:size(S_dat,1)
                    s_electrode_resp=squeeze(S_dat(i,:,:))';
                    n_electrode_resp=squeeze(N_dat(i,:,:))';
                    s_n_rho=SN_corr_dat(i);
                    s_n_rho_rnd=SN_cor_rnd_dat(i,:);
                    is_sig=SN_sig_dat(i);
                    p_ratio=SN_p_ratio_dat(i);
                    word_pos=repmat(1:size(s_electrode_resp,1),size(s_electrode_resp,2),1)';
                    %f=figure(fix((i-1)/total_plots)+1);
                    
                    %set(f,'position',[1123          29        1266        1275])
                    sup_title=(strcat(obj.elec_ch_label{i,1}));
                    ax=subplot(num_rows,num_columns,3*(i-num_rows*fix((i-1)/num_rows)-1)+1);
                    
                    % plot the distbirtubiotn
                    h1=histogram(s_n_rho_rnd,nbins);
                    hold on
                    xline(s_n_rho,'linewidth',3);
                    ax.Box='off'
                    h1.EdgeColor='w';
                    shg
                    ax.XAxis.LineWidth=2;
                    ax.YAxis.LineWidth=2;
                    ax.Title.String=sprintf('corr=%f,\n p_{ratio}=%0.4f sig=%d',s_n_rho,p_ratio,is_sig);
                    %
                    ax=subplot(num_rows,num_columns,3*(i-num_rows*fix((i-1)/num_rows)-1)+2);
                    %x=.3;y=s_electrode_resp(:);
                    b1=plot(mean(word_pos,2)+.1,mean(s_electrode_resp,2),'color',[1,.5,.5],'linewidth',2,'marker','o','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,0],'displayname','S');
                    hold on
                    b1=plot(mean(word_pos,2)-.1,mean(n_electrode_resp,2),'color',[.5,.5,1],'linewidth',2,'marker','o','MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,1],'displayname','N');
                    y=s_electrode_resp;
                    x=word_pos;
                    bl=errorbar(mean(x,2)+.1,mean(y,2),std(y,[],2)./sqrt(size(y,2)));
                    bl.LineStyle='none';
                    bl.Color=[1,.5,.5];
                    bl.LineWidth=2;
                    bl.CapSize=2;
                    hAnnotation = get(bl,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
                    %
                    y=n_electrode_resp;
                    x=word_pos;
                    bl=errorbar(mean(x,2)-.1,mean(y,2),std(y,[],2)./sqrt(size(y,2)));
                    bl.LineStyle='none';
                    bl.Color=[.5,.5,1];
                    bl.LineWidth=2;
                    bl.CapSize=2;
                    hAnnotation = get(bl,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
                    ax.XAxis.Visible = 'on';
                    ax.XTick=1:max(word_pos(:));
                    ax.XLim=[0,max(word_pos(:))+1];
                    %ax.XTick=1:size(electrode_resp,1);
                    all_points=[s_electrode_resp(:);n_electrode_resp(:)];
                    y_quantile=quantile(all_points,10);
                    h=get(ax,'children');
                    ax.FontSize=12;
                    set(ax,'ydir', 'normal','box','off','ylim',[y_quantile(1),y_quantile(end)]);
                    ah =get(ax,'children');
                    arrayfun(@(x) set(ah(x),'DisplayName',''),[1:2]);
                    ah(3).DisplayName='N';
                    ah(4).DisplayName='S';
                    set(ax,'children',ah);
                    legend(ah(3:4),'Location','northwest','NumColumns',2)
                    xlabel('word position');
                    ax.YLabel.String='High Gamma (a.u.)';
                    ax.XAxis.LineWidth=2;
                    ax.YAxis.LineWidth=2;
                    ax.Title.String=erase(sup_title,'_');
                    if ~mod(i,num_rows) | i==size(S_dat,1)
                        %legend('show','Location','northeastoutside')
                        pp=pp+1;
                        if ~exist(strcat(analysis_path,obj.subject_id))
                            mkdir(strcat(analysis_path,obj.subject_id));
                        end
                        set(gcf,'PaperPosition',[.25 .25 8 6])
                        set(gcf,'PaperOrientation','landscape');
                        fname=sprintf('%s_s_v_n_words-%d-%d_p_%0.2f_%s_%s_%s_unipolar_1.pdf',obj.subject_id,...
                        min(ops.words),...
                        max(ops.words),...
                        ops.threshold,...
                        ops.side,...
                        obj.modality,...
                        num2str(num2str(pp)));
                        print(f, '-bestfit','-dpdf','-opengl', strcat(analysis_path,obj.subject_id,'/',fname));
                        %close(f)
                        set(0, 'CurrentFigure', f);
                        clf reset;
                    end
                    
                end
                close(f);
               
            end
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
            sent_aud_name=aud_name(contains(aud_name,'English'));
            sent_aud_transc=aud_transcript(contains(aud_name,'English'));
            %
            condtion_flag='S';
            cond_resp=obj.get_cond_resp(condtion_flag);
            % identify and remove sentence repeats 
            % https://www.mathworks.com/matlabcentral/answers/197246-finding-duplicate-strings-in-a-cell-array-and-their-index
            cond_key='word_';
            trial_strings=cellfun(@(x) strjoin({x.string{contains(x.key,cond_key)}}),cond_resp,'uni',false);
            [D,P,X] = unique(trial_strings,'stable');
            assert(length(D)==60);
            % create average responses
            % get only one of the repeats 
            uniq_cond_resp=cond_resp(P);
            uniq_str=cellfun(@(x) strjoin({x.string{contains(x.key,cond_key)}}),uniq_cond_resp,'uni',false);
            cellfun(@(x,y) assert(strcmp(x,y)),unique(uniq_str,'stable'),uniq_str);
         
             % add sentence name to table 
            uniq_sent_aud_name=sent_aud_name(P);
            if obj.modality=='auditory'
                uniq_sent_aud_name=erase(uniq_sent_aud_name,'_48000')
            end 
            assert(sum(cellfun(@(x,y) contains(x,y),erase(sent_aud_transc(P),','),uniq_str))>=58);
            %A=[erase(sent_aud_transc(P),','),uniq_str,cellfun(@(x,y) {contains(x,y)},erase(sent_aud_transc(P),','),uniq_str)]
            % use unique str 
            uniq_cond_resp_mod=arrayfun(@(x) [cell2table(repmat(uniq_sent_aud_name(x),size(uniq_cond_resp{x},1),1),"VariableNames",{'sentence_name'}),uniq_cond_resp{x}], 1:length(uniq_cond_resp),'uni',false)';
            
             % create average responses
            B=uniq_cond_resp;
            B=obj.get_average(B);
            B=arrayfun(@(x) [cell2table(repmat(uniq_sent_aud_name(x),size(B{x},1),1),"VariableNames",{'sentence_file'}),B{x}], 1:length(B),'uni',false)';
            word_data=obj.get_value(B,'key','word','type','contain');
            word_dat_mod=arrayfun(@(x) [cell2table(repmat(uniq_sent_aud_name(x),size(word_data{x},1),1),"VariableNames",{'sentence_name'}),word_data{x}], 1:length(word_data),'uni',false)';
            sent_dat=vertcat(word_dat_mod{:});
            sent_id=extract(sent_dat.sentence_name,digitsPattern);
            sent_id=num2cell(cellfun(@(x) str2num(x),sent_id));
            word_id=extract(sent_dat.key,digitsPattern);
            word_id=num2cell(cellfun(@(x) str2num(x),word_id));
            stim_resp_table=[cell2table([sent_id,word_id],"VariableNames",{'sentence_id','word_id'}),sent_dat];
            % make some assertion 
            aa=stim_resp_table(stim_resp_table.sentence_id==13,'string');
            assert(strcmp(aa.string{10},'would'));
            
            aa=stim_resp_table(stim_resp_table.sentence_id==31,'string');
            assert(strcmp(aa.string{9},'driver'));
            
            
            if not(isprop(obj,'stim_resp_table'))
                P = addprop(obj,'stim_resp_table');
            end 
            obj.stim_resp_table=stim_resp_table;
            % save a structure version to use in python 
            save_struct=struct;
            save_struct.s_vs_n_struct=table2struct(obj.s_vs_n_sig);
            save_struct.s_vs_n_struct.elec_ch_label=obj.elec_ch_label;
            save_struct.s_vs_n_struct.bip_elec_ch_label=...
                arrayfun(@(X) sprintf('%s-%s',obj.bip_ch_label_valid{X,1},obj.bip_ch_label_valid{X,2}),1:size(obj.bip_ch_label_valid,1),'uni',false);
            save_struct.stim_resp_struct=table2struct(obj.stim_resp_table);
            save_file=sprintf('%s/%s_%s_%s_stim_resp_struct',obj.save_path,obj.subject_id,'SN',obj.modality)
            save(save_file,'save_struct');
            
            python_path='/Users/eghbalhosseini/anaconda3/envs/neural_nlp_1/bin/python';
             commandStr = sprintf('%s /Users/eghbalhosseini/MyCodes/ecog_SN/construct_stim_resp_data.py %s %s_%s',python_path,obj.subject_id,'SN',obj.modality);
            [status, commandOut] = system(commandStr);
            if status==0
                fprintf('python conversion was successful\n');
                fprintf(commandOut)
            else 
                fprintf(commandOut)
            end
            
            
%             B=cond_resp;
%             B=obj.get_average(B);
%             word_data=obj.get_value(B,'key','word','type','contain');
%             % flatten the word data 
%             values=cellfun(@(x) x(:,~ismember(x.Properties.VariableNames,{'key','string'})), word_data,'uni',false);
%             % reformat the data to have a cell with averages  
%             func=@(x) cell2mat(transpose(x));
%             func_1= @(x) cellfun(@(var_name) func(x.(var_name)),x.Properties.VariableNames,'uni',false);
%             values_cell=cellfun(@(x) func_1(x), values,'uni',false);
%             % make sure it did what it was supposed to do 
%             t_id=randi(size(values_cell,1));
%             assert(all(all(values_cell{t_id}{1}==cell2mat(values{t_id}.elec_data_dec'))));
%             % combine along word dimension
%             word_values=arrayfun(@(y) cell2mat(cellfun(@(x) x{y}, values_cell,'uni',false)'),[1:size(values_cell{1},2)],'uni',false);
%             word_values_parsed=cellfun(@(x) transpose(mat2cell(x,size(x,1),ones(1,size(x,2)))),word_values,'uni',false);
%             % add additional meta data and save
%             % add strings
%             strings=cellfun(@(x) x(:,ismember(x.Properties.VariableNames,{'string'})), word_data,'uni',false);
%             string_values=cat(2,cellfun(@(x) x.string,strings,'uni',false));
%             string_values=cellflat(string_values);
%             % add keys 
%             keys=cellfun(@(x) x(:,ismember(x.Properties.VariableNames,{'key'})), word_data,'uni',false);
%             key_values=cat(2,cellfun(@(x) x.key,keys,'uni',false));
%             key_values=cellflat(key_values);
%             % add sentence_id
%             s_id_values=cat(2,arrayfun(@(x) x*ones(1,max(size(keys{x}))),1:size(keys,1),'uni',false));
%             s_id_values=cell2mat(s_id_values);
%             s_id_values=mat2cell(s_id_values,1,ones(1,size(s_id_values,2)));
%             % make sure it did what it was supposed to do 
%             overlap_id=find(cell2mat(s_id_values)==t_id);
%             assert(all(all(cell2mat(word_values_parsed{1}(overlap_id)')==cell2mat(values{t_id}.elec_data_dec'))));
%             % combine everything to a table 
%             variable_names=cellflat({'key','string','sentence_id',setdiff(word_data{1}.Properties.VariableNames, {'key','string'},'stable')});
%             %temp_table=cell2table(horzcat({key_values},{string_values},{s_id_values},{D'},word_values),'VariableNames',variable_names);
%             stim_resp_table=cell2table(horzcat(key_values',string_values',s_id_values',[word_values_parsed{:}]),'VariableNames',variable_names);
%             % add stim_resp_table to the object 
%             if not(isprop(obj,'stim_resp_table'))
%                 P = addprop(obj,'stim_resp_table');
%             end 
%             obj.stim_resp_table=stim_resp_table;
%             % save a structure version to use in python 
%             save_struct=struct;
%             save_struct.s_vs_n_struct=table2struct(obj.s_vs_n_sig);
%             % add electrode names 
%             save_struct.s_vs_n_struct.elec_ch_label=obj.elec_ch_label;
%             save_struct.s_vs_n_struct.bip_elec_ch_label=...
%                 arrayfun(@(X) sprintf('%s-%s',obj.bip_ch_label_valid{X,1},obj.bip_ch_label_valid{X,2}),1:size(obj.bip_ch_label_valid,1),'uni',false);
%             save_struct.stim_resp_struct=table2struct(obj.stim_resp_table);
%             save_file=sprintf('%s/%s_%s_%s_stim_resp_struct',obj.save_path,obj.subject_id,obj.experiment,obj.modality)
%             save(save_file,'save_struct');
%             % run python 
%             python_path='/Users/eghbalhosseini/anaconda3/envs/neural_nlp/bin/python';
%             commandStr = sprintf('%s /Users/eghbalhosseini/MyCodes/ecog_SN/construct_stim_resp_data.py %s %s_%s',python_path,obj.subject_id,obj.experiment,obj.modality);
%             [status, commandOut] = system(commandStr);
%             if status==0
%                 fprintf('python conversion was successful\n');
%                 fprintf(commandOut)
%             else 
%                 fprintf(commandOut)
%             end
            
    end
    
    
    function obj=make_extended_stim_resp_mat(obj,varargin)
        p=inputParser();
            addParameter(p, 'words', 'all');
            addParameter(p,'python_path','/Users/eghbalhosseini/miniconda3/envs/neural_align/bin/python')
            addParameter(p,'exec_path','/Users/eghbalhosseini/miniconda3/envs/neural_align/bin/python')
            parse(p, varargin{:});
            ops = p.Results;
            aud_name=obj.events_table.final_audio_filename;
            aud_name=erase(aud_name,'.wav');
            aud_name=erase(aud_name,'_48000');
            aud_transcript=obj.events_table.final_audio_transcript;
            % replace 2 buggs sentences
            cond_key='word_';
            stim_strings=cellfun(@(x) strjoin({x.string{contains(x.key,cond_key)}}),obj.trial_data,'uni',false);
            sent_13=find(contains(aud_name,'013_Eng'));
            aud_transcript{sent_13};
            assert(sum(contains(strsplit(aud_transcript{sent_13}),strsplit(stim_strings{sent_13})))==11);
            aud_transcript{sent_13}=stim_strings{sent_13};
            
            sent_31=find(contains(aud_name,'031_Eng'));
            aud_transcript{sent_31};
            assert(sum(contains(strsplit(aud_transcript{sent_31}),strsplit(stim_strings{sent_31})))==11);
            aud_transcript{sent_31}=stim_strings{sent_31};
            P=1:length(aud_transcript);
            A=[erase(aud_transcript(P),','),stim_strings,cellfun(@(x,y) {contains(x,y)},erase(aud_transcript(P),','),stim_strings)];
            assert(all([A{:,3}]))
            
            
            list_id=obj.events_table.final_list;
            trial_id=obj.events_table.trial;
            B=obj.get_average(obj.trial_data);
            word_data=obj.get_value(B,'key','word','type','contain');
            word_dat_mod=arrayfun(@(x) [cell2table(repmat(aud_name(x),size(word_data{x},1),1),"VariableNames",{'stim_name'}),word_data{x}], 1:length(word_data),'uni',false)';
            word_dat_mod=arrayfun(@(x) [cell2table(repmat({trial_id(x)},size(word_dat_mod{x},1),1),"VariableNames",{'Trial_id'}),word_dat_mod{x}], 1:length(word_dat_mod),'uni',false)';
            word_dat_mod=arrayfun(@(x) [cell2table(repmat({list_id(x)},size(word_dat_mod{x},1),1),"VariableNames",{'list_id'}),word_dat_mod{x}], 1:length(word_dat_mod),'uni',false)';
            
            all_dat=vertcat(word_dat_mod{:});
            stim_id=extract(all_dat.stim_name,digitsPattern);
            stim_id=num2cell(cellfun(@(x) str2num(x),stim_id));
            word_id=extract(all_dat.key,digitsPattern);
            word_id=num2cell(cellfun(@(x) str2num(x),word_id));
            stim_types={'S','N'};
            stim_type=cellfun(@(x) stim_types{isempty(regexp(x,'English'))+1},all_dat.stim_name,'uni',false);
            stim_value=cellfun(@(x) aud_transcript{find(ismember(aud_name,x),1,'first')},all_dat.stim_name,'uni',false);
            stim_resp_table=[cell2table([stim_id,word_id,stim_type,stim_value],"VariableNames",{'stim_id','word_id','stim_type','stim_value'}),all_dat];
            
            if not(isprop(obj,'extended_stim_resp_table'))
                P = addprop(obj,'extended_stim_resp_table');
            end 
            obj.extended_stim_resp_table=stim_resp_table;
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
            save_file=sprintf('%s/%s_%s_%s_extended_stim_resp_struct',obj.save_path,obj.subject_id,'ECOG_SN',obj.modality);
            save(save_file,'save_struct','-v7');
            
            
            python_path=ops.python_path;
            exec_path=ops.exec_path;
             commandStr = sprintf('%s %s %s %s_%s %s',python_path,exec_path,obj.subject_id,'ECOG_SN',obj.modality,obj.save_path);
            [status, commandOut] = system(commandStr);
            if status==0
                fprintf('python conversion was successful\n');
                fprintf(commandOut)
            else 
                fprintf('python conversion was NOT SUCCESSFUL\n');
                fprintf(commandOut)
            end
            
    end     
    
    function plot_trial_timecourse(obj,varargin)
        p=inputParser();
        addParameter(p, 'signal_type', 'envelope_dec');
        parse(p, varargin{:});
        ops = p.Results;
        elec_flag=ops.signal_type;
        analysis_path=strcat(obj.save_path,'analysis/trial_timecourse/');
        save_dir=fullfile(analysis_path,obj.subject_id,elec_flag);
        if ~exist(save_dir)
           mkdir(save_dir);
        end
        % get s_v_n
        if isprop(obj,'s_vs_n_sig')
                s_vs_n=obj.s_vs_n_sig.(elec_flag){1};
        end 
        
        if contains(elec_flag,'bip')
            elec_labels=obj.bip_ch_label_valid;
            elec_labels=arrayfun(@(x) [elec_labels{x,1},'-',elec_labels{x,2}],1:size(elec_labels,1),'uni',false);
        else
            elec_labels=obj.elec_ch_label;
        end 
        
        condtion_flag='S';
        [S_ave_tbl,S_table]=obj.get_ave_cond_trial('words',[1:12],'condition',condtion_flag);
        
        
        condtion_flag='N';
        [N_ave_tbl,N_table]=obj.get_ave_cond_trial('words',[1:12],'condition',condtion_flag);

        S_dat=S_table.(elec_flag){1};
        N_dat=N_table.(elec_flag){1};
        cond_resp=obj.get_cond_resp('S');
        cond_key='fix';
        fix_act=cellfun(@(x) x.(elec_flag){ismember(x.key,cond_key)},cond_resp,'uni',false);
        sizes=(cellfun(@size, fix_act,'uni',false));
        max_second_dim=max(cell2mat(sizes));
        nan_act=cellfun(@(x) nan*ones(x(1),max_second_dim(2)-x(2)),sizes,'uni',false);
        fix_nan_act=cellfun(@(x,y) horzcat(x,y),fix_act,nan_act,'uni',false);
        word_indv_act=[fix_nan_act];

        words_acts=[fix_act];

        for kk=1:12
            cond_key=['word_',num2str(kk)];
            word_act=cellfun(@(x) x.(elec_flag){contains(x.key,cond_key)},cond_resp,'uni',false);
            sizes=(cellfun(@size, word_act,'uni',false));
            max_second_dim=max(cell2mat(sizes));
            nan_act=cellfun(@(x) nan*ones(x(1),max_second_dim(2)-x(2)),sizes,'uni',false);
            word_nan_act=cellfun(@(x,y) horzcat(x,y),word_act,nan_act,'uni',false);
            word_indv_act=[word_indv_act,word_nan_act];
            words_acts=[words_acts,word_act];
        end
        act_index=cumsum(cellfun(@(x) size(x,2), words_acts),2);
        words_act=arrayfun(@(x) horzcat(words_acts{x,:}),1:size(words_acts,1),'uni',false)';
        % add nan
        sizes=(cellfun(@size, words_act,'uni',false));
        max_second_dim=max(cell2mat(sizes));
        nan_acts=cellfun(@(x) nan*ones(x(1),max_second_dim(2)-x(2)),sizes,'uni',false);
        words_nan_act=cellfun(@(x,y) horzcat(x,y),words_act,nan_acts,'uni',false);
        word_tensor=cat(3,words_nan_act{:});
        Sentence_tensor=word_tensor;
        Sentence_word_act=word_indv_act;
        act_index_sent=act_index;
        %% 
        cond_resp=obj.get_cond_resp('N');
        cond_key='fix';
        fix_act=cellfun(@(x) x.(elec_flag){ismember(x.key,cond_key)},cond_resp,'uni',false);
        sizes=(cellfun(@size, fix_act,'uni',false));
        max_second_dim=max(cell2mat(sizes));
        nan_act=cellfun(@(x) nan*ones(x(1),max_second_dim(2)-x(2)),sizes,'uni',false);
        fix_nan_act=cellfun(@(x,y) horzcat(x,y),fix_act,nan_act,'uni',false);
        word_indv_act=[fix_nan_act];
        words_acts=[fix_act];
        for kk=1:12
            cond_key=['word_',num2str(kk)];
            word_act=cellfun(@(x) x.(elec_flag){contains(x.key,cond_key)},cond_resp,'uni',false);
            sizes=(cellfun(@size, word_act,'uni',false));
            max_second_dim=max(cell2mat(sizes));
            nan_act=cellfun(@(x) nan*ones(x(1),max_second_dim(2)-x(2)),sizes,'uni',false);
            word_nan_act=cellfun(@(x,y) horzcat(x,y),word_act,nan_act,'uni',false);
            word_indv_act=[word_indv_act,word_nan_act];
            words_acts=[words_acts,word_act];
        end
        act_index=cumsum(cellfun(@(x) size(x,2), words_acts),2);
        words_act=arrayfun(@(x) horzcat(words_acts{x,:}),1:size(words_acts,1),'uni',false)';
        % add nan
        sizes=(cellfun(@size, words_act,'uni',false));
        max_second_dim=max(cell2mat(sizes));
        nan_acts=cellfun(@(x) nan*ones(x(1),max_second_dim(2)-x(2)),sizes,'uni',false);
        words_nan_act=cellfun(@(x,y) horzcat(x,y),words_act,nan_acts,'uni',false);
        word_tensor=cat(3,words_nan_act{:});
        Nonword_tensor=word_tensor;
        Nonword_word_act=word_indv_act;
        act_index_nonw=act_index;
        %% 
        pbar=ProgressBar(size(Sentence_tensor,1));
        for elec_id=1:size(Sentence_tensor,1)
            %% 
            D=squeeze(Sentence_tensor(elec_id,:,:));
            [~,srt_idx] = sort(sum(isnan(D)));
            D=D(:,srt_idx);
            act_index=act_index_sent(srt_idx,:);
            D_sent=D;
            act_idx_sent=act_index;
            %
            D=squeeze(Nonword_tensor(elec_id,:,:));
            [~,srt_idx] = sort(sum(isnan(D)));
            D=D(:,srt_idx);
            act_index=act_index_nonw(srt_idx,:);
            act_idx_nonw=act_index;
            D_nonw=D;
            
            min_max=prctile([D_sent(:);D_nonw(:)],[2 98]);
            
            f=figure('visible','off');
            
            %clf;
            set(gcf,'position',[31,1,1713,1010]);
            set(gcf,'position',[-2261 73 2162 1091]);
            %ax=axes('position',[.05,.25,.42,.7]);
            ax=axes('position',[.05,.2,.25,.74]);
            hold on
            D=D_sent;
            act_index=act_idx_sent;
            color_map=mat2cell(inferno(2*size(D,2)),ones(1,2*size(D,2)),3);
            
            D_cell=mat2cell(D,size(D,1),ones(1,size(D,2)));
            D_norm_cell=cellfun(@(x) (x-min_max(1))./(min_max(2)-min_max(1)),D_cell,'UniformOutput',false);
            D_norm=cell2mat(D_norm_cell);
            for x=1:size(D_norm,2)
                hold on
                plot(D_norm(:,x)+x,'color',color_map{x},'linewidth',1)
                plot(act_index(x,:),D_norm(act_index(x,:),x)+x,'k','linestyle','none','marker','o','markerfacecolor',[.3,.3,.3],'markersize',2)
                %yline(x,'linewidth',.5,'color',[.5,.5,.5]);
                %plot([0 size(D,1)],[x,x],'linewidth',.5,'color',[.5,.5,.5]);
            end
            set(ax,'ytick',[1:size(D,1)]);
            set(ax,'yticklabel','');
            set(ax,'ylim',[0,size(D,2)+1]);
            ax.XAxis.TickLength=[0.005,0.01];
            ax.YAxis.TickLength=[0.005,0.01];
            ax.Title.String=[elec_labels{elec_id},', ', strrep(elec_flag,'_','-')];
            set(ax,'xlim',[0 size(D,1)]);
            pos=get(ax,'position');
            Newpos=[pos(1) pos(2)-0.1 pos(3) 0.05];
            
            
            ax=axes('position',[.05,.03,.25,.15]);
            for x=1:size(D_norm,2)
                hold on
                plot(D_norm(:,x),'color',[color_map{x},.2],'linewidth',.5)
            end
            s_resp=nanmean(D_norm,2);
            b1=plot(s_resp,'color',[0,0,0],'linewidth',2,'MarkerEdgeColor',[0,0,0],'displayname','S');
            set(ax,'ylim',min_max);
            set(ax,'xlim',[0 size(D_norm,1)]);
            plot([0 size(D_norm,1)],[0,0],'color',[.3,.3,.3],'linewidth',1,'linestyle','--')
            ax=axes('position',[.65,.8,.3,.18]);
            D=cellfun(@(x) x(elec_id,:),Sentence_word_act,'uni',false);
            jumps=cumsum([0,cellfun(@length,D(1,:) )]);
            for jj=1:size(D,2)
                D_word=D(:,jj);
                D_word_mat=cell2mat(D_word);
                for x=1:size(D_word,1)
                    hold on
                    plot(jj*10+[1:length(D_word{x})]+jumps(jj),D_word{x},'color',[color_map{x},.2],'linewidth',.5)
                end
                plot(jj*10+[1:size(D_word_mat,2)]+jumps(jj),nanmean(D_word_mat,1),'color',[0,0,0,1],'linewidth',2)
            end
            set(ax,'ylim',min_max);
            set(ax,'xlim',[0,10*size(D,2)+length([D{1,:}])]);
            plot(get(gca,'xlim'),[0,0],'color',[.3,.3,.3],'linewidth',1,'linestyle','--')
            
            % 
            %ax=axes('position',[.55,.25,.42,.7]);
            ax=axes('position',[.35,.2,.25,.74]);
            D=D_nonw;
            hold on
            color_map=mat2cell(viridis(2*size(D,2)),ones(1,2*size(D,2)),3);
            act_index=act_idx_nonw;
            D_cell=mat2cell(D,size(D,1),ones(1,size(D,2)));
            D_norm_cell=cellfun(@(x) (x-min_max(1))./(min_max(2)-min_max(1)),D_cell,'UniformOutput',false);
            D_norm=cell2mat(D_norm_cell);
            for x=1:size(D_norm,2)
                hold on
                plot(D_norm(:,x)+x,'color',color_map{x},'linewidth',1)
                plot(act_index(x,:),D_norm(act_index(x,:),x)+x,'k','linestyle','none','marker','o','markerfacecolor',[.3,.3,.3],'markersize',2)
                %yline(x,'linewidth',.5,'color',[.5,.5,.5]);
                %plot([0 size(D,1)],[x,x],'linewidth',.5,'color',[.5,.5,.5]);
            end
            set(ax,'ytick',[1:size(D,1)]);
            set(ax,'yticklabel','');
            set(ax,'ylim',[0,size(D,2)+1]);
            ax.XAxis.TickLength=[0.005,0.01];
            ax.YAxis.TickLength=[0.005,0.01];
            set(ax,'xlim',[0 size(D,1)]);
            pos=get(ax,'position');
            Newpos=[pos(1) pos(2)-0.1 pos(3) 0.05];
            
            ax=axes('position',[.35,.03,.25,.15]);
            for x=1:size(D_norm,2)
                hold on
                plot(D_norm(:,x),'color',[color_map{x},.2],'linewidth',.5)
            end
            n_resp=nanmean(D_norm,2);
            b1=plot(n_resp,'color',[0,0,0],'linewidth',2,'MarkerEdgeColor',[0,0,0],'displayname','N');
            set(ax,'ylim',min_max);
            set(ax,'xlim',[0 size(D_norm,1)]);
            plot(get(gca,'xlim'),[0,0],'color',[.3,.3,.3],'linewidth',1,'linestyle','--')
            %
            ax=axes('position',[.65,.6,.3,.18]);
            D=cellfun(@(x) x(elec_id,:),Nonword_word_act,'uni',false);
            jumps=cumsum([0,cellfun(@length,D(1,:) )]);
            for jj=1:size(D,2)
                D_word=D(:,jj);
                D_word_mat=cell2mat(D_word);
                for x=1:size(D_word,1)
                    hold on
                    plot(jj*10+[1:length(D_word{x})]+jumps(jj),D_word{x},'color',[color_map{x},.2],'linewidth',.5)
                end
                plot(jj*10+[1:size(D_word_mat,2)]+jumps(jj),nanmean(D_word_mat,1),'color',[0,0,0,1],'linewidth',2)
            end
            set(ax,'ylim',min_max);
            set(ax,'xlim',[0,10*size(D,2)+length([D{1,:}])]);
            
            plot(get(gca,'xlim'),[0,0],'color',[.3,.3,.3],'linewidth',1,'linestyle','--')

            %
            % plot average over words
            %ax=axes('position',[.05,.05,.42,.18]);
            ax=axes('position',[.65,.05,.3,.18]);
            D=cellfun(@(x) x(elec_id,:),Nonword_word_act,'uni',false);
            nonw_fix_resp=nanmean(cell2mat(D(:,1)),2);
            
            D=cellfun(@(x) x(elec_id,:),Sentence_word_act,'uni',false);
            sent_fix_resp=nanmean(cell2mat(D(:,1)),2);
            
            
            s_electrode_resp=vertcat(sent_fix_resp',squeeze(S_dat(elec_id,:,:))');
            n_electrode_resp=vertcat(nonw_fix_resp',squeeze(N_dat(elec_id,:,:))');
            word_pos=repmat(1:size(s_electrode_resp,1),size(s_electrode_resp,2),1)'-1;
             %append fixation respon 
            b1=plot(mean(word_pos,2)+.1,mean(s_electrode_resp,2),'color',[1,.5,.5],'linewidth',1,'marker','o','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,0],'displayname','S');
            hold on
            b1=plot(mean(word_pos,2)-.1,mean(n_electrode_resp,2),'color',[.5,.5,1],'linewidth',1,'marker','o','MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,1],'displayname','N');
            y=s_electrode_resp;
            x=word_pos;
            bl=errorbar(mean(x,2)+.1,mean(y,2),std(y,[],2)./sqrt(size(y,2)));
            bl.LineStyle='none';
            bl.Color=[1,.5,.5];
            bl.LineWidth=2;
            bl.CapSize=2;
            hAnnotation = get(bl,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
            %
            y=n_electrode_resp;
            x=word_pos;
            bl=errorbar(mean(x,2)-.1,mean(y,2),std(y,[],2)./sqrt(size(y,2)));
            bl.LineStyle='none';
            bl.Color=[.5,.5,1];
            bl.LineWidth=2;
            bl.CapSize=2;
            hAnnotation = get(bl,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
            ax.XAxis.Visible = 'on';
            ax.XTick=0:max(word_pos(:));
            ax.XTickLabel{1}='Fix';
            ax.XLim=[-1,max(word_pos(:))+1];
            %ax.XTick=1:size(electrode_resp,1);
            all_points=[s_electrode_resp(:);n_electrode_resp(:)];
            y_quantile=quantile(all_points,10);
            h=get(ax,'children');
            ax.FontSize=12;
            set(ax,'ydir', 'normal','box','off','ylim',[y_quantile(1),y_quantile(end)]);
            %set(ax,'ydir', 'normal','box','off','ylim',min_max);
            ah =get(ax,'children');
            arrayfun(@(x) set(ah(x),'DisplayName',''),[1:2]);
            ah(3).DisplayName='N';
            ah(4).DisplayName='S';
            set(ax,'children',ah);
            legend(ah(3:4),'Location','northwest','NumColumns',2)
            xlabel('word position');
            ax.YLabel.String='High Gamma (a.u.)';
            ax.XAxis.LineWidth=2;
            ax.YAxis.LineWidth=2;
            ax.Title.String=['S>N, ', num2str(s_vs_n(elec_id))];
            % plot over time not word_id
            % ax=axes('position',[.55,.05,.42,.18]);
            ax=axes('position',[.65,.32,.3,.18]);
            s_resp=nanmean(D_sent,2);
            n_resp=nanmean(D_nonw,2);
            min_max=prctile([s_resp(:);n_resp(:)],[1 99]);
            b1=plot(s_resp,'color',[1,.5,.5],'linewidth',1,'MarkerEdgeColor',[1,0,0],'displayname','S');
            hold on
            b1=plot(n_resp,'color',[.5,.5,1],'linewidth',1,'MarkerEdgeColor',[0,0,1],'displayname','N');
            set(ax,'ydir', 'normal','box','off','ylim',min_max);
            set(findall(gcf,'-property','FontSize'),'FontSize',8);
            set(gcf,'PaperPosition',[.25 .25 8 10])
            set(gcf,'PaperOrientation','landscape');
            fname=sprintf('%s_%s_timecourse.pdf',obj.subject_id,elec_labels{elec_id});            
            print(f, '-fillpage','-dpdf','-opengl', fullfile(save_dir,fname));
            fname=sprintf('%s_%s_timecourse.png',obj.subject_id,elec_labels{elec_id});            
            %print(f,'-dpng','-painters', fullfile(save_dir,fname));
                        
            close(f);
            close all
            pbar.step([],[],[]);
            
            %waitfor(findobj('type','figure','number',1));
            %% 
        end
        
       
    end
    
    function plot_per_shank(obj,varargin)
        p=inputParser();
        addParameter(p, 'signal_type', 'envelope_dec');
        parse(p, varargin{:});
        ops = p.Results;
        elec_flag=ops.signal_type;
        analysis_path=strcat(obj.save_path,'analysis/per_shank/');
        save_dir=fullfile(analysis_path,obj.subject_id,[obj.modality,'_',elec_flag]);
        if ~exist(save_dir)
           mkdir(save_dir);
        end
        % get s_v_n
        if isprop(obj,'s_vs_n_sig')
                s_vs_n=obj.s_vs_n_sig.(elec_flag){1};
        end 
        
        if contains(elec_flag,'bip')
            elec_labels=obj.bip_ch_label_valid;
            elec_labels=arrayfun(@(x) [elec_labels{x,1},'-',elec_labels{x,2}],1:size(elec_labels,1),'uni',false);
        else
            elec_labels=obj.elec_ch_label;
        end 
        
        condtion_flag='S';
        [S_ave_tbl,S_table]=obj.get_ave_cond_trial('words',[1:12],'condition',condtion_flag);
        cond_resp=obj.get_cond_resp(condtion_flag);
        S_dat=S_table.(elec_flag){1};
        cond_key='fix';
        fix_act=cellfun(@(x) x.(elec_flag){ismember(x.key,cond_key)},cond_resp,'uni',false);
        fix_act=cellfun(@(x) nanmean(x,2),fix_act,'uni',false);
        fix_act=cat(2,fix_act{:});
        S_dat=cat(3,fix_act,S_dat);
        
        
        
        condtion_flag='N';
        [N_ave_tbl,N_table]=obj.get_ave_cond_trial('words',[1:12],'condition',condtion_flag);
        cond_resp=obj.get_cond_resp(condtion_flag);
        N_dat=N_table.(elec_flag){1};
        cond_key='fix';
        fix_act=cellfun(@(x) x.(elec_flag){ismember(x.key,cond_key)},cond_resp,'uni',false);
        fix_act=cellfun(@(x) nanmean(x,2),fix_act,'uni',false);
        fix_act=cat(2,fix_act{:});
        N_dat=cat(3,fix_act,N_dat);
        
        
        
        pat = lettersPattern;
        elec_names=vertcat(cellfun(@(x) strsplit(x,'-'),elec_labels,'uni',false )');
        elec_names=cellfun(@(x) x(1),elec_names);
        elec_unique=cellfun(@(x) extract(x,pat) , elec_names);
        
        
        shanks=unique(elec_unique,'stable');
        [shanks,ia,ic] = unique(elec_unique,'stable');
        
        elec_per_shank_counts = accumarray(ic,1);
        
        
        f=figure('visible','off');            
        %clf;
        set(gcf,'position',[-2261 73 2162 1091]);
        subplot_index = reshape(1:length(shanks)*max(elec_per_shank_counts),  length(shanks),max(elec_per_shank_counts)).';
        tiledlayout(max(elec_per_shank_counts),length(shanks),'Padding','tight','TileSpacing','tight');
        %pbar=ProgressBar(length(elec_labels));
        
        for idx=1:length(shanks)
            index_set=subplot_index(:,idx);
            electrodes_in_group=(ic==idx);
            shank_labels=elec_labels(electrodes_in_group);
            for id_col=1:length(shank_labels)
                %ax=subplot(max(elec_per_shank_counts),length(shanks),index_set(id_col));
                ax=nexttile(index_set(id_col));
                elec_id=find(ismember(elec_labels,shank_labels{id_col}));
                assert(length(elec_id)==1);
                s_electrode_resp=squeeze(S_dat(elec_id,:,:))';
                n_electrode_resp=squeeze(N_dat(elec_id,:,:))';
                word_pos=repmat(1:size(s_electrode_resp,1),size(s_electrode_resp,2),1)'-1;
                
                b1=plot(mean(word_pos,2)+.1,mean(s_electrode_resp,2),'color',[1,.5,.5],'linewidth',.5,'marker','o','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,0],'displayname','S','markersize',1);
                hold on
                b1=plot(mean(word_pos,2)-.1,mean(n_electrode_resp,2),'color',[.5,.5,1],'linewidth',.5,'marker','o','MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,1],'displayname','N','markersize',1);
                y=s_electrode_resp;
                x=word_pos;
                bl=errorbar(mean(x,2)+.1,mean(y,2),std(y,[],2)./sqrt(size(y,2)));
                bl.LineStyle='none';
                bl.Color=[1,.5,.5];
                bl.LineWidth=.25;
                bl.CapSize=0;
                hAnnotation = get(bl,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
                %
                y=n_electrode_resp;
                x=word_pos;
                bl=errorbar(mean(x,2)-.1,mean(y,2),std(y,[],2)./sqrt(size(y,2)));
                bl.LineStyle='none';
                bl.Color=[.5,.5,1];
                bl.LineWidth=.25;
                bl.CapSize=0;
                hAnnotation = get(bl,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
                ax.XAxis.Visible = 'off';
                ax.YAxis.Visible = 'off';
                ax.XTick=[];
                ax.XLim=[-1,max(word_pos(:))+1];
                all_points=[s_electrode_resp(:);n_electrode_resp(:)];
                y_quantile=quantile(all_points,6);
                h=get(ax,'children');
                ax.FontSize=4;
                set(ax,'ydir', 'normal','box','off','ylim',[y_quantile(1),y_quantile(end)]);
                ax.Title.String=shank_labels{id_col};
                ax.Title.FontSize=5;
                if s_vs_n(elec_id)==1
                    ax.Title.FontWeight='bold';
                    ax.XAxis.Visible = 'on';
                    ax.YAxis.Visible = 'on';
                    ax.XTick=[];
                    ax.YTick=[];
                    ax.Box='on';
                    ax.LineWidth=1.25;
                    ax.Title.FontSize=5;
                    
                    
                else
                    ax.Title.FontWeight='normal';
                    ax.Title.FontSize=5;
                end
                if idx==1 && id_col==1
                    ax.XAxis.Visible = 'on';
                    ax.XTick=0:2:max(word_pos(:));
                    ax.XTickLabel{1}='Fix';
                    lgd=legend;
                    lgd.FontSize=5;
                    lgd.FontWeight='bold';
                    lgd.NumColumns=1;
                    xx1=lgd.Position;
                    lgd.Position=[xx1(1),xx1(2),.015,.015];
                end 
                %pbar.step([],[],[]);
            end
        end 
        set(gcf,'PaperOrientation','landscape');
        fname=sprintf('%s_%s_%s_per_shank.pdf',obj.subject_id,obj.modality,ops.signal_type);            
        print(f, '-fillpage','-dpdf','-painters', fullfile(save_dir,fname));
        close(f);
            close all
    end 
   
    end 
end 

