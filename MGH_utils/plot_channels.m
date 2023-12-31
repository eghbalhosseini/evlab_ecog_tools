function plot_channels(D,S_rate,color_map,channel_label,valid_channels,varargin)
p=inputParser();
addParameter(p, 'timing', []);
addParameter(p,'events',[]);
parse(p, varargin{:});
ops = p.Results;

if ~exist('color_map','var')
    color_map=mat2cell(inferno(2*size(D,2)),ones(1,2*size(D,2)),3);
end 

if ~exist('valid_channels','var')
    valid_channels=ones(1,size(D,2));
end 



%PLOT_CHANNELS Summary of this function goes here
%   Detailed explanation goes here
down_sample_rate=8;
t_len=60;
t_length=t_len*S_rate/down_sample_rate; % 20 sec * 200 downsample rate
x_cell=mat2cell(D',ones(1,size(D,2)));
x_cell=cellfun(@(x) downsample(x,down_sample_rate),x_cell,'uni',false);
min_max=mean(cell2mat(cellfun(@(y) prctile(y,[5 95]),x_cell,'UniformOutput',false)));
x_norm_cell=cellfun(@(x) (x-min_max(1))./(min_max(2)-min_max(1)),x_cell,'UniformOutput',false);
x_norm_cell=arrayfun(@(x)  x_norm_cell{x}*valid_channels(x),1:size(x_norm_cell,1),'uni',false)';
%
figure(1);
clf;
set(gcf,'position',[31,1,1713,1010]);
ax=axes('position',[.05,.1,.93,.88]);
hold on
time_stamps=([1:size(x_cell{1},2)]/S_rate)*down_sample_rate;
hold on

for x=1:size(x_norm_cell,1)
    hold on 
    plot(time_stamps,x_norm_cell{x}+x,'color',color_map{x},'linewidth',1.5,'tag',sprintf('ch %d, tag %s',x,channel_label{x}))
end 

H=arrayfun(@(x) plot(time_stamps,x_norm_cell{x}+x,'color',color_map{x},'linewidth',1.5,'tag',sprintf('ch %d, tag %s',x,channel_label{x})),[1:size(x_norm_cell,1)]);

if ~isempty(ops.timing)
    timings_ds=cell2mat(cellfun(@(x) floor(x./S_rate),ops.timing,'uni',false));
     for kk=1:length(timings_ds)
        event=timings_ds(kk,:);
        corners_x=[min(event),max(event),max(event),min(event)];
        corners_y=[ax.YLim(1),ax.YLim(1),ax.YLim(2),ax.YLim(2)];
        pat=patch(corners_x,corners_y,[1,.4,.4]);
        pat.FaceAlpha=.1;
    end 
end 

set(gcf,'doublebuffer','on');
set(ax,'ytick',[1:size(x_norm_cell,1)]);
set(ax,'yticklabel','');

set(ax,'ylim',[0,size(x_norm_cell,1)+4]);
ax.XAxis.TickLength=[0.005,0.01];
ax.YAxis.TickLength=[0.005,0.01];
set(ax,'xlim',[0 t_len]);
pos=get(ax,'position');
Newpos=[pos(1) pos(2)-0.1 pos(3) 0.05];
xmax=max(time_stamps);
S=['set(gca,''xlim'',get(gcbo,''value'')+[0 ' num2str(t_len) '])'];

h=uicontrol('style','slider',...
    'units','normalized','position',Newpos,...
    'callback',S,'min',0,'max',xmax-t_len);

datacursormode on
dcm = datacursormode(gcf);
set(dcm,'UpdateFcn',@myupdatefcn)
end

