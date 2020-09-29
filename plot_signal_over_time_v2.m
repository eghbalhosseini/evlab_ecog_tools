function []= plot_signal_over_time_v2(signal,param,session_start,plot_title,save_plots,print_figure,save_path, figure_label, sub_id)
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
    
    x_mean=cellfun(@nanmean,x_cell);
    x_std=cellfun(@nanstd,x_cell);
    std_factor=8;
    outlier=arrayfun(@(x) (((x_cell{x}-x_mean(x))<-std_factor*x_std(x))...
            | ((x_cell{x}-x_mean(x))>std_factor*x_std(x))),1:length(x_mean),'uni',false);

    min_max=nanmean(cell2mat(cellfun(@(y) [min(y),max(y)],x_cell,'UniformOutput',false)));
    x_norm_cell=cellfun(@(x) (x-min_max(1))./(.5*(min_max(2)-min_max(1))),x_cell,'UniformOutput',false);
    %
    x_cell_nan=cellfun(@(x) x*nan,x_norm_cell,'uni',false);
    for i=1:length(x_mean)
        x_cell_nan{i}(outlier{i})=x_norm_cell{i}(outlier{i});
    end 
    
    figure(1);
    col_inf=inferno(floor(.8*size(x_norm_cell,1)));
    col_vir=viridis(floor(.8*size(x_norm_cell,1)));
    colors=[col_vir(1:floor(size(x_norm_cell,1)/2),:);col_inf(1:(floor(size(x_norm_cell,1)/2)+1),:)];

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
        arrayfun(@(x) plot(t_window,x_cell_nan{x}(t_window)+x-1,'color','r','linewidth',2),param.channels)
        set(ax,'ytick',[1:size(x_norm_cell,1)]);
        set(ax,'yticklabel','')
        arrayfun(@(x) text(t_window(1),x,[num2str(x),' '],'Color',colors(x,:),'HorizontalAlignment','right','VerticalAlignment','middle'),param.channels);
        set(ax,'ylim',[0,size(x_norm_cell,1)+1]);
        set(ax,'TickDir','in');
        if any(ismember(t_window,session_start))
        arrayfun(@(x) plot([t_window(x),t_window(x)],[0,size(x_norm_cell,1)+1],'k-','LineWidth',2),find(ismember(t_window,session_start)))
        session_name=num2str(find(ismember(session_start,t_window)));
        for pp=1:length(session_name)
            session_name_temp = session_name(pp);
            arrayfun(@(x) text(t_window(x),size(x_norm_cell,1)+1,['sess: ',session_name_temp],'VerticalAlignment','bottom','fontsize',20),find(ismember(t_window,session_start)))
        end
        end 
        shg;
        title({plot_title,[num2str(kk),' of ', num2str(kk_total)]})  
        pause
    end

    
    if(print_figure && save_plots)
    print_length=5;
    f = figure('visible','off');
        set(0,'units','pixels');
        screen = get(0,'ScreenSize');
        set(f,'color',[.7,.7,.7],'position', screen)
        ax=axes('position',[.02,.02,.95,.95]);
        for kk=1:kk_total
            hold on
            session_end=min([(kk*t_length),size(signal,1)]);
            t_window=((kk-1)*t_length+1):2:(session_end);
            arrayfun(@(x) plot(t_window,x_norm_cell{x}(t_window)+x-1,'color',colors(x,:)),param.channels)
            arrayfun(@(x) plot(t_window,x_cell_nan{x}(t_window)+x-1,'color','r','linewidth',2),param.channels)
            if any(ismember(t_window,session_start))
                arrayfun(@(x) plot([t_window(x),t_window(x)],[0,size(x_norm_cell,1)+1],'k-','LineWidth',2),find(ismember(t_window,session_start)))
                session_name=num2str(find(ismember(session_start,t_window)));
                for pp=1:length(session_name)
                    session_name_temp = session_name(pp);
                    arrayfun(@(x) text(t_window(x),size(x_norm_cell,1)+1,['sess: ',session_name_temp],'VerticalAlignment','bottom','fontsize',20),find(ismember(t_window,session_start)))
                end
                if kk==1 
                        set(ax,'ytick',[1:size(x_norm_cell,1)]);
            set(ax,'yticklabel','')
            arrayfun(@(x) text(t_window(1),x,[num2str(x),' '],'Color',colors(x,:),'HorizontalAlignment','right','VerticalAlignment','middle'),param.channels);
            set(ax,'ylim',[0,size(x_norm_cell,1)+1]);
            set(ax,'TickDir','in');
            title({plot_title})
                end 
            end

        end
        set(ax, 'XLimSpec', 'Tight');
       filename = strcat(save_path, 'channel_plots', filesep, sub_id, figure_label,'.pdf')
       get(f,'Papersize');
       f.Units='Inches';
       set(f,'Papersize',[15,print_length*kk_total]);
       set(f,'PaperOrientation','landscape')
       print(f,'-opengl','-dpdf','-fillpage',filename)
    end 
    %close all;
end
