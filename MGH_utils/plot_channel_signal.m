function TrigMat1=plot_channel_signal(record,hdr,chanel_id)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


figure;
for k=1:size(chanel_id,1)
    ax=axes('position',[.1,.8/16*k+0.05,.8,.8/16]);
    colormap(flipud(gray()));
    caxis([0,3e4]);
    chan_id=find(ismember(hdr.label,chanel_id{k}));
    sig_mic_R=record(chan_id,:);
    plot(sig_mic_R,'k-','linewidth',2);
    %ylim([min(sig_mic_R)*-1-.2,max(sig_mic_R)+.2]);
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    ax.YLabel.String=hdr.label{chan_id};
    ax.YLabel.FontSize=8;
    ax.YLabel.Rotation=0;
    hold on 
   
end
end