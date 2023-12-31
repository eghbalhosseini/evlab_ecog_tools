function [outputArg1,outputArg2] = comm_src_ave(signal,ops)
%COMM_SRC_AVE Summary of this function goes here
%   Detailed explanation goes here
%% do common source averaging based on tags (not sure if this should be done in this case)
op_ecog_labels=ops.ecog_channels_labels;
ch_num=cellfun(@str2num,extract(op_ecog_labels,digitsPattern));
ch_tags=extract(op_ecog_labels,lettersPattern);
[ecog_tags,~,~]=unique(ch_tags,'stable');
ecog_loc_for_tag=cellfun(@(x) find(contains(op_ecog_labels,x))', ecog_tags,'uni',false);

% 
signal_ca=nan*zeros(size(signal));
for kk=1:size(ecog_loc_for_tag,1)
    same_shank=ecog_loc_for_tag{kk};
    same_shank_clean=intersect(same_shank,ops.ecog_channels_selected);
    if ~isempty(same_shank_clean)
        com_ave=mean(signal(:,same_shank_clean),2);
        %com_ave_glob=mean(signal_glob(:,same_shank_clean),2);
    else 
        warning('no clean electrode for common source average')
        com_ave=zeros(size(signal,1),1);
        %com_ave_glob=zeros(size(signal,1),1);
    end 
    signal_ca(:,same_shank)=signal(:,same_shank)-com_ave;
    %signal_ca_glob(:,same_shank)=signal_glob(:,same_shank)-com_ave_glob;
end
if ops.do_shankMean
    signal=signal_ca;
    param.preproc_step=[param.preproc_step;'shank_common_source_removal'];
end 
end

