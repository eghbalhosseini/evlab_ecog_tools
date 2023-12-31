function [ops] = find_epileptic_elec(signal,ops)
%FIND_EPILEPTIC_CH Summary of this function goes here
%   Detailed explanation goes here

    detectionIEDs          = []; % output from Janca et al. script - 1 cell array per segment
    detectionIEDs.settings = '-k1 3.65 -h 60 -dec 200 -dt 0.005 -pt 0.12 -ti 1'; % if you change "-dec 200" here, do not forget to change in selectChannels_Using ... below
    detectionIEDs.segments = [];
    
    detectionIEDs=automaticSpikeDetection_UsingJancaMethod(signal, ops.sr, detectionIEDs.settings);
    
    % plot the detections
    % plot IED events
    t_len=100;
    t_length=t_len*200; % 20 sec * 200 downsample rate
    kk_total=ceil(size(detectionIEDs.envelope,1)/t_length);
    x_cell=mat2cell(detectionIEDs.envelope',ones(1,size(detectionIEDs.envelope,2)));
    min_max=mean(cell2mat(cellfun(@(y) prctile(y,[1 99]),x_cell,'UniformOutput',false)));
    x_norm_cell=cellfun(@(x) (x-min_max(1))./(min_max(2)-min_max(1)),x_cell,'UniformOutput',false);
    
    
    col_inf=inferno(floor(.8*size(x_norm_cell,1)));
    col_vir=viridis(floor(.8*size(x_norm_cell,1)));
    colors=[col_vir(1:floor(size(x_norm_cell,1)/2),:);col_inf(1:(floor(size(x_norm_cell,1)/2)+1),:)];
    figure(1);
    clf;
    set(gcf,'position',[31,1,1713,1010]);
    ax=axes('position',[.05,.1,.93,.88]);
    hold on
    time_stamps=[1:size(detectionIEDs.envelope,1)]/200;
    hold on
    H=arrayfun(@(x) plot(time_stamps,x_norm_cell{x}+x,'color',colors(x,:),'tag',sprintf('ch %d, tag %s',x,ops.elecnames{x})),[1:size(x_norm_cell,1)]);
   
    spike_chan = detectionIEDs.out.chan;
    [spike_chan_sort,sort_idx]=sort(spike_chan);
    spike_times_sort= detectionIEDs.out.pos(sort_idx);
    
    %spike_dur = detectionIEDs.out.dur;
    yval=cell2mat(arrayfun(@(x) x_norm_cell{spike_chan_sort(x)}(floor(spike_times_sort(x)*200))+spike_chan_sort(x),1:length(spike_chan_sort),'uni',false));
    H1=scatter(spike_times_sort,yval,50,'MarkerEdgeColor',[0 0 0],...
        'MarkerFaceColor',[0 .7 .7],...
        'LineWidth',2,'MarkerFaceAlpha',.5);
    
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
    
    waitfor(findobj('type','figure','number',1));
    
    %
    
    % From this automatic assessment, selection of the final pool of channels
    detectionIEDs.tableChanSelection = [];  % info regarding chan selection
    detectionIEDs.threshold          = 6.5; % channels with IEDs higher than threshold are removed
    
    
    
    currIEDs.fs             = 200; % default downsampling during automatic detection (-dec 200)
    currIEDs.discharges.MV  = [];
    currIEDs.numSamples     = 0;
    
    
    currIEDs.discharges.MV = detectionIEDs.discharges.MV;
    currIEDs.numSamples= currIEDs.numSamples + size(detectionIEDs.d_decim, 1);
    
    
    % Compute number of detected spike per channel: [c x 1] where c channels
    numSpikes     = []; numSpikes     = sum(currIEDs.discharges.MV==1, 1);
    totalDuration = []; totalDuration = (currIEDs.numSamples / currIEDs.fs) / 60; % in minutes
    numSpikes_min = []; numSpikes_min = numSpikes / totalDuration;
    numSpikes     = transpose(numSpikes);
    numSpikes_min = transpose(numSpikes_min);
    
    % Select channels with IEDs / minute below threshold - [c x 1] where c channels
    indChanSelected = [];
    indChanSelected = find(numSpikes_min < detectionIEDs.threshold);
    tableChanSelection.numSpikesAll           = numSpikes_min;
    tableChanSelection.indChansSelected       = indChanSelected;
    tableChanSelection.indChansDeselected     = setdiff(ops.elecids,indChanSelected);
    tableChanSelection.nameChansSelected      = transpose(ops.elecnames(indChanSelected));
    tableChanSelection.numSpikesChansSelected = numSpikes_min(indChanSelected);
    
    
    ops.IED_results=tableChanSelection;
    ops.IED_results.threshold=detectionIEDs.threshold;
    %
    
    ops.ecog_channels_IED_deselected=tableChanSelection.indChansDeselected;
    




end

