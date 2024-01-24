function [signal_out,ops] = comm_src_ave(signal,ops)
%COMM_SRC_AVE Summary of this function goes here
% do common source averaging based on tags (not sure if this should be done in this case)
op_ecog_labels=ops.ecog_channels_labels;
ch_num=cellfun(@str2num,extract(op_ecog_labels,digitsPattern));
ch_tags=extract(op_ecog_labels,lettersPattern);
[ecog_tags,~,~]=unique(ch_tags,'stable');
ecog_loc_for_tag=cellfun(@(x) find(contains(op_ecog_labels,x))', ecog_tags,'uni',false);

% 
signal_out=nan*zeros(size(signal));
for kk=1:size(ecog_loc_for_tag,1)
    same_shank=ecog_loc_for_tag{kk};
    shank_label=extract(op_ecog_labels(same_shank),lettersPattern);
    assert(size(unique(shank_label),1)==1)
    same_shank_clean=intersect(same_shank,ops.ecog_channels_selected);
    if ~isempty(same_shank_clean)
        com_ave=mean(signal(:,same_shank_clean),2);
        %com_ave_glob=mean(signal_glob(:,same_shank_clean),2);
    else 
        warning('no clean electrode for common source average')
        com_ave=zeros(size(signal,1),1);
        %com_ave_glob=zeros(size(signal,1),1);
    end 
    signal_out(:,same_shank)=signal(:,same_shank)-com_ave;
    %signal_ca_glob(:,same_shank)=signal_glob(:,same_shank)-com_ave_glob;
end
    ops.perfromed_steps=[ops.perfromed_steps;'shank_common_source_removal'];
end

