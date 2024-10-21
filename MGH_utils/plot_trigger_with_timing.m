function plot_trigger_with_timing(TrigMat1,all_trial_timing,trial_based_frame)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

figure;
ax=axes('position',[.1,.1,.85,.85]);
for k=1:size(TrigMat1,2)
    colormap(flipud(gray()));
    caxis([0,3e4]);
    plot(.8*TrigMat1(:,k)+k,'k-','linewidth',2);
    hold on
end
axis tight
set(gca,'xtick',[]);
set(gca,'ytick',[]);
ax.YLabel.String=num2str(k);
ax.YLabel.FontSize=8;
ax.YLabel.Rotation=0;
hold on

corners_x=[min(trial_based_frame),max(trial_based_frame),max(trial_based_frame),min(trial_based_frame)];
corners_y=[ax.YLim(1),ax.YLim(1),ax.YLim(2),ax.YLim(2)];
pat=patch(corners_x,corners_y,[1,0,0]);
pat.FaceAlpha=.1;

for pp=1:size(all_trial_timing,1)
    sm_event=all_trial_timing{pp,1};
    corners_x=[min(sm_event),max(sm_event),max(sm_event),min(sm_event)];
    corners_y=[ax.YLim(1),ax.YLim(1),ax.YLim(2),ax.YLim(2)];
    pat=patch(corners_x,corners_y,[0,0,1]);
    pat.FaceAlpha=.1;
    pat.EdgeColor='none';
end


end