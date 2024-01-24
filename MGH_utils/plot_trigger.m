function TrigMat1=plot_trigger(TrigMat1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
TrigMat1 = dec2bin(TrigMat1,16); % there are 16 digital channels
TrigMat1 = double(TrigMat1); % convert characters to numbers
TrigMat1 = TrigMat1-48; % 48 -> '0' and 49 -> '1'
TrigMat1 = fliplr(TrigMat1); % bit-16 is the first channel


figure;
for k=1:size(TrigMat1,2)
    ax=axes('position',[.1,.8/size(TrigMat1,2)*k+0.05,.8,.8/size(TrigMat1,2)]);
    colormap(flipud(gray()));
    caxis([0,3e4]);
    plot(TrigMat1(:,k),'k-','linewidth',2);
    ylim([min(TrigMat1(:,k))*-1-.2,max(TrigMat1(:,k))+.2]);
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    ax.YLabel.String=num2str(k);
    ax.YLabel.FontSize=8;
    ax.YLabel.Rotation=0;
    hold on 
   
end
end