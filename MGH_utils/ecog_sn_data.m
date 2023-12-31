
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
            uniq_cond_resp_mod=arrayfun(@(x) [cell2table(repmat(uniq_sent_aud_name(x),size(uniq_cond_resp{x},1),1),"VariableNames",{'sentence_id'}),uniq_cond_resp{x}], 1:length(uniq_cond_resp),'uni',false)';
            
             % create average responses
            B=uniq_cond_resp;
            B=obj.get_average(B);
            word_data=obj.get_value(B,'key','word','type','contain');
            word_dat_mod=arrayfun(@(x) [cell2table(repmat(uniq_sent_aud_name(x),size(word_data{x},1),1),"VariableNames",{'sentence_name'}),word_data{x}], 1:length(word_data),'uni',false)';
            sent_dat=vertcat(word_dat_mod{:});
            sent_id=extract(sent_dat.sentence_name,digitsPattern);
            sent_id=num2cell(cellfun(@(x) str2num(x),sent_id));
            word_id=extract(sent_dat.key,digitsPattern);
            word_id=num2cell(cellfun(@(x) str2num(x),word_id));
            stim_resp_table=[cell2table([sent_id,word_id],"VariableNames",{'sentence_id','word_id'}),sent_dat];
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
    end 
end

