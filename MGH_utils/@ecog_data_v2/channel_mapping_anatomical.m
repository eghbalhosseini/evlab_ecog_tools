function [chan_mapping,chan_labels]=channel_mapping_anatomical(obj,vera_mat)
    % helper function for mapping channels in datafiles to channels in anatomical files
    % TODO - modify to be less manual

    fprintf(1, '\n> Adding anatomical files to object\n');

    % load channel labels
    chan_labels = obj.elec_ch_label;
    chan_types = obj.elec_ch_type;

    if strcmp(obj.subject,'BJH006')
        % map localized channels to channels in object
        reduced_localized_names = cellfun(@(x) split(x,'-'),vera_mat.electrodeNames,'UniformOutput',false);
        localized_names = cellfun(@(x) strcat(x{1},'_',extract(x{end},digitsPattern)),reduced_localized_names,'UniformOutput',false); 
        localized_names = cellfun(@(x) x{1}, localized_names,'UniformOutput',false);
        chan_mapping = cellfun(@(x) find(strcmp(x,localized_names)),chan_labels,'UniformOutput',false);
        chan_labels = cellfun(@(x) localized_names(x),chan_mapping,'UniformOutput',false);

    elseif strcmp(obj.subject,'SLCH002')
        % map localized channels to channels in object
        reduced_localized_names = cellfun(@(x) split(x,'^'),vera_mat.electrodeNames,'UniformOutput',false);
        localized_names = cellfun(@(x) strcat(x{1},'_',extract(x{2},digitsPattern)),reduced_localized_names,'UniformOutput',false); 
        localized_names = cellfun(@(x) x{1}, localized_names,'UniformOutput',false);
        chan_mapping = cellfun(@(x) find(strcmp(x,localized_names)),chan_labels,'UniformOutput',false);
        chan_labels = cellfun(@(x) localized_names(x),chan_mapping,'UniformOutput',false);

    elseif contains(obj.subject,'BJH') 
        % map localized channels to channels in object
        first_thing = cellfun(@(x) x(1:3),vera_mat.electrodeNames,'UniformOutput',false);
        second_thing = cellfun(@(x) extract(x,digitsPattern),vera_mat.electrodeNames,'UniformOutput',false); 
        second_thing = cellfun(@(x) x{end},second_thing,'UniformOutput',false);
        localized_names = cellfun(@(x,y) strcat(x([1,3]),'_',y),first_thing,second_thing,'UniformOutput',false);
        % localized_names = cellfun(@(x) x{1}, localized_names,'UniformOutput',false);
        if strcmp(obj.subject,'BJH008')
            for i=1:length(localized_names)
                if contains(localized_names{i},'ER')
                    localized_names(i) = {['O' localized_names{i}(3:end)]};
                elseif contains(localized_names{i},'FR')
                    localized_names(i) = {['P' localized_names{i}(3:end)]};
                elseif contains(localized_names{i},'FR')
                    localized_names(i) = {['P' localized_names{i}(3:end)]};
                else
                    localized_names(i) = {localized_names{i}([1,3:end])};
                end
            end
        end
        chan_mapping = cellfun(@(x) find(strcmp(x,localized_names)),chan_labels,'UniformOutput',false);
        chan_labels = cellfun(@(x) localized_names(x),chan_mapping,'UniformOutput',false);

    else % AMC and MCJ subjects
        num_chans = sum(contains(chan_types,'ecog') | strcmp('seeg',chan_types));
        % assert(num_chans==length(vera_mat.electrodeNames),'Mapping may not be correct!');
        chan_mapping = num2cell(1:num_chans)';
        chan_labels = chan_labels(1:num_chans);
        % add empy cells to end of mapping/labels
        empty_to_add = cell(size(obj.elec_data,1)-num_chans,1);
        chan_mapping = [chan_mapping; empty_to_add];
        chan_labels = [chan_labels; empty_to_add];
    end
    
end