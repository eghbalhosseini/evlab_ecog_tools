function plot_trigger_with_events(TrigMat1,all_trial_timing,events_table,num_sess)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
figure;
for ses=1:num_sess
    b1=find(events_table.final_list==ses);
    event_b1=events_table(b1,:);
    b1_tri=all_trial_timing(setdiff(b1,[]));
    sig_range=[floor(min([b1_tri{:,1}])):floor(max([b1_tri{:,1}]))];

    %ax=axes('position',[.1,.1,.85,.85]);
    ax=subplot(num_sess,1,ses);
    for k=1:size(TrigMat1,2)
        colormap(flipud(gray()));
        caxis([0,3e4]);
        plot([sig_range-min(sig_range)]/1024,.8*TrigMat1(sig_range,k)+k,'k-','linewidth',2);
        hold on ;
    end
    %set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        ax.YLabel.String=num2str(k);
        ax.YLabel.FontSize=8;
        ax.YLabel.Rotation=0;
        hold on 
        
    hold on 
    
    
    for pp=1:size(b1_tri,1)
        sm_event=[b1_tri{pp,1}-min(sig_range)]/1024;
        corners_x=[min(sm_event),max(sm_event),max(sm_event),min(sm_event)];
        corners_y=[ax.YLim(1),ax.YLim(1),ax.YLim(2),ax.YLim(2)];
        pat=patch(corners_x,corners_y,[0,0,1]);
        pat.FaceAlpha=.1;
        pat.EdgeColor='none';
        text(mean(sm_event),mean(ax.YLim),num2str(pp),'fontsize',12,'horizontalalignment','center')
    end 
    
     event_aud_st=event_b1.trial_onset;
     event_aud_end=event_b1.audio_ended;
    for pp=1:size(b1_tri,1)
        sm_event=[event_aud_st(pp),event_aud_end(pp)];
        corners_x=[min(sm_event),max(sm_event),max(sm_event),min(sm_event)];
        corners_y=[ax.YLim(1),ax.YLim(1),ax.YLim(2),ax.YLim(2)];
        pat=patch(corners_x,corners_y,[0,1,1]);
        pat.FaceAlpha=0;
        pat.EdgeColor='k';
        text(mean(sm_event),mean(ax.YLim)+2,num2str(pp),'fontsize',12,'horizontalalignment','center')
    end 
    end
end 