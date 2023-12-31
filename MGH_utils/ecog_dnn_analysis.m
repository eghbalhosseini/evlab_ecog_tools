classdef ecog_dnn_analysis < ecog_analysis
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        
%% method for constructing the object 
        function dnn_obj = ecog_dnn_analysis(ecog_dnn_data,ecog_sn_data)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            dnn_obj@ecog_analysis(ecog_dnn_data,ecog_sn_data);
        end
        
        
        function average_timecourse_across_electrodes(obj,varargin)
            p=inputParser();
            addParameter(p, 'average_brain', false);
            addParameter(p,'elec_mode','bip_elec'); % elec or bip_elec
            addParameter(p,'zscore',false);
            parse(p, varargin{:});
            ops = p.Results;
            % get significant electrodes 
            sig_elec_tbl=obj.ecog_data.s_vs_n_sig;
            sig_elec_loc=cell2mat(sig_elec_tbl.([ops.elec_mode,'_data_dec']));
            % get pials 
            lh_pial=obj.ecog_data.lh_pial;
            rh_pial=obj.ecog_data.rh_pial;
            % get electrode locations 
            switch ops.elec_mode
                case 'bip_elec'
                    elec_loc=obj.ecog_data.bip_ch_pos_anat;
                    elec_label=obj.ecog_data.biop_ch_label_valid;
                    elec_label=cellfun(@(x,y) erase([x,'-',y],'_'),elec_label(:,1),elec_label(:,2),'uni',false);
                    if ops.zscore
                        elec_flag='bip_elec_data_zs_dec';
                    else 
                        elec_flag='bip_elec_data_dec';
                    end 
                case 'elec'
                    elec_loc=obj.ecog_data.elec_ch_pos_anat;
                    elec_label=cellfun(@(x) erase(x,'_'),obj.ecog_data.elec_ch_label,'uni',false);
                    if ops.zscore
                        elec_flag='elec_data_zs_dec';
                    else 
                        elec_flag='elec_data_dec';
                    end 
            end
            % do some checks first 
            assert(size(elec_loc,1)==size(sig_elec_loc,1));
            assert(all(cellfun(@(x,y) strcmp(x,y),elec_loc.label,elec_label)),'labels dont match!');
            sig_elec_tbl=elec_loc(sig_elec_loc,:);
            % 
            % get condition responses 
            condtion_flag='S';
            [S_table]=obj.ecog_data.get_ave_cond_resp('condition',condtion_flag);
            condtion_flag='N';
            [N_table]=obj.ecog_data.get_ave_cond_resp('condition',condtion_flag);
            % 
            func=@(x) cell2mat(permute(x,[3,2,1]));
            S_dat=func(S_table.(elec_flag));
            N_dat=func(N_table.(elec_flag));
            sig_elec=logical(cell2mat(obj.ecog_data.s_vs_n_sig.(elec_flag)));
            
            s_electrode_resp=(nanmean(S_dat(sig_elec,:,:),2));
            s_electrode_resp=reshape(s_electrode_resp,[],size(s_electrode_resp,3));
            % 
            n_electrode_resp=(nanmean(N_dat(sig_elec,:,:),2));
            n_electrode_resp=reshape(n_electrode_resp,[],size(n_electrode_resp,3));
            %% 
             %% plot responses per position ove all electrodes: 
            f=figure(1);
            clf(f);
            f.Units='Inches';
            f.PaperOrientation='landscape';
            f.Position=[9. 1.   11   8.5];
            f.Color=[1,1,1];
            % plot left pial 
            ax=subplot('position',[.1,.6,.3,.3]);
            hold on
            nvertices=size(lh_pial.vert,1);
            color_data = repmat([255,229,204]./255,nvertices,1);
            patch_handle = patch('vertices', lh_pial.vert, 'Faces', lh_pial.face, 'FaceVertexCData', color_data,'FaceLighting','gouraud','SpecularStrength',0,'DiffuseStrength',0.7);
            shading interp;
            xlim(1.8*50*[-1 1]); ylim(1.8*50*[-1 1]-20); zlim(1.8*50*[-1 1]+40);
            view([-90 0]);
            patch_handle.FaceAlpha=.4;
            light_handle = camlight('left','infinite');
            set(light_handle, 'Position', [-1, 1, 0.33]);
            ax.XAxis.Visible='off';
            ax.YAxis.Visible='off';
            ax.ZAxis.Visible='off';
            daspect([1,1,1])
            camzoom(1.3);
            % plot electrodes on the left side 
            L_labels=cellfun(@(x) x(1)=='L',elec_label(sig_elec_loc));
            L_loc=cellfun(@(x) x(1)=='L',sig_elec_tbl.label);
            if any(L_loc)
                left_tbl=sig_elec_tbl(L_loc,:);
                for kk=1:size(left_tbl,1)
                    R=left_tbl(kk,:).R; A=left_tbl(kk,:).A; S=left_tbl(kk,:).S;
                    h=plot3(R,A,S);
                    h.Marker='o';h.MarkerFaceColor=[0,0,0];h.MarkerEdgeColor=[0,0,0];h.MarkerSize=10;
                    plot3(R*[1,1,1],A*[1,1,1],S+[-5,0,5],'LineWidth',2,'color',[1,0,0]);
                    plot3(R+[-5,0,5],A*[1,1,1],S+[1,1,1],'LineWidth',2,'color',[0,1,0]);
                    plot3(R*[1,1,1],A+[-5,0,5],S+[1,1,1],'LineWidth',2,'color',[0,0,1]);
                end 
            end         
            
            % plot right pial 
            ax=subplot('position',[.5,.6,.3,.3]);
            hold on 
            nvertices=size(rh_pial.vert,1);
            color_data = repmat([255,229,204]./255,nvertices,1);
            patch_handle = patch('vertices', rh_pial.vert, 'Faces', rh_pial.face, 'FaceVertexCData', color_data,'FaceLighting','gouraud','SpecularStrength',0,'DiffuseStrength',0.7);
            shading interp;
            xlim(1.8*50*[-1 1]); ylim(1.8*50*[-1 1]-20); zlim(1.8*50*[-1 1]+40);
            view([90 0]);
            patch_handle.FaceAlpha=.4;
            light_handle = camlight('right','infinite');
            ax.XAxis.Visible='off';
            ax.YAxis.Visible='off';
            ax.ZAxis.Visible='off';
            daspect([1,1,1])
            camzoom(1.3);
            % plot electrodes on the right side 
            R_labels=cellfun(@(x) x(1)=='R',elec_label(sig_elec_loc));
            R_loc=cellfun(@(x) x(1)=='R',sig_elec_tbl.label);
            if any(R_loc)
                right_tbl=sig_elec_tbl(R_loc,:);
                for kk=1:size(right_tbl,1)
                    R=right_tbl(kk,:).R; A=right_tbl(kk,:).A; S=right_tbl(kk,:).S;
                    h=plot3(R,A,S);
                    h.Marker='o';h.MarkerFaceColor=[0,0,0];h.MarkerEdgeColor=[0,0,0];h.MarkerSize=10;
                    plot3(R*[1,1,1],A*[1,1,1],S+[-5,0,5],'LineWidth',4,'color',[1,0,0]);
                    plot3(R+[-5,0,5],A*[1,1,1],S+[1,1,1],'LineWidth',4,'color',[0,1,0]);
                    plot3(R*[1,1,1],A+[-5,0,5],S+[1,1,1],'LineWidth',4,'color',[0,0,1]);
                end 
            end 
            
            %%
            ax=subplot('position',[.1,.1,.6,.4])
            s_word_pos=repmat(1:size(s_electrode_resp,2),size(s_electrode_resp,1),1);
            n_word_pos=repmat(1:size(n_electrode_resp,2),size(n_electrode_resp,1),1);
            b1=plot(mean(s_word_pos,1)+.1,mean(s_electrode_resp,1),'color',[1,.5,.5],'linewidth',2,'marker','o','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[1,0,0],'displayname','S');
            hold on
            b1=plot(mean(n_word_pos,1)-.1,mean(n_electrode_resp,1),'color',[.5,.5,1],'linewidth',2,'marker','o','MarkerFaceColor',[0,0,1],'MarkerEdgeColor',[0,0,1],'displayname','N');
            y=s_electrode_resp;
            x=s_word_pos;
            bl=errorbar(mean(x,1)+.1,mean(y,1),std(y,[],1)./sqrt(size(y,1)));
            bl.LineStyle='none';
            bl.Color=[1,.5,.5];
            bl.LineWidth=2;
            bl.CapSize=2;
            hAnnotation = get(bl,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
            %                     %
            y=n_electrode_resp;
            x=n_word_pos;
            bl=errorbar(mean(x,1)-.1,mean(y,1),std(y,[],1)./sqrt(size(y,1)));
            bl.LineStyle='none';
            bl.Color=[.5,.5,1];
            bl.LineWidth=2;
            bl.CapSize=2;
            hAnnotation = get(bl,'Annotation');hLegendEntry = get(hAnnotation,'LegendInformation');set(hLegendEntry,'IconDisplayStyle','off');
            ax.XAxis.Visible = 'on';
            ax.XTick=1:max(s_word_pos(:));
            ax.XLim=[0,max(s_word_pos(:))+1];
            %ax.XTick=1:size(electrode_resp,1);
            all_points=[s_electrode_resp(:);n_electrode_resp(:)];
            y_quantile=quantile(all_points,10);
            h=get(ax,'children');
            ax.FontSize=12;
            set(ax,'ydir', 'normal','box','off');%,[y_quantile(4),y_quantile(7)]);
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
            analysis_path=strcat(obj.ecog_data.save_path,'analysis/all_elec_s_vs_n/');
            if ~exist(strcat(analysis_path,obj.ecog_data.subject_id))
                mkdir(strcat(analysis_path,obj.ecog_data.subject_id));
            end
            set(gcf,'PaperPosition',[.25 .25 8 6]);
            set(gcf,'PaperOrientation','landscape');
            print(f, '-bestfit','-dpdf','-opengl', strcat(analysis_path,obj.ecog_data.subject_id,'/',obj.ecog_data.subject_id,'_',obj.ecog_data.experiment,'_lang_electrode_resp_',obj.ecog_data.modality,'_',ops.elec_mode,'_',ops.elec_mode,'_zs_',num2str(ops.zscore),'.pdf'));
            close(f)
            
            
        end
    end
end

