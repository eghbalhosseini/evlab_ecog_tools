function [trial_data_epoch, trial_bip_data_epoch] = extract_trial_epochs(obj,epoch_args)
    arguments
        obj ecog_data_v2
        epoch_args.epoch_tw  = [-0.5 3];
        epoch_args.key  = 'fix';
        epoch_args.selectChannels = 1:size(obj.elec_data,1);
    end
    epoch_tw_samples = epoch_args.epoch_tw.*obj.sample_freq;
   
    fprintf(1, '\n> Cutting signal into trial epoch ... \n');
    fprintf(1,'[');
    trial_bip_data_epoch = [];
    trial_data_epoch = [];
    % go through each trial timing table
    for k = 1:size(obj.trial_timing)
        trial_time_tbl = obj.trial_timing{k};
        probe_key=find(ismember(trial_time_tbl.key,epoch_args.key));
        assert(length(probe_key)==1);
        trial_elec_data = arrayfun(@(x) obj.elec_data(epoch_args.selectChannels,trial_time_tbl(x,:).start+epoch_tw_samples(1):trial_time_tbl(x,:).start+epoch_tw_samples(2)),probe_key,'uni',false)';
    
        if ~isempty(obj.bip_elec_data)
            trial_bip_elec_data = arrayfun(@(x) obj.bip_elec_data(:,trial_time_tbl(x,:).start+epoch_tw_samples(1):trial_time_tbl(x,:).start+epoch_tw_samples(2)),probe_key,'uni',false)';
            trial_bip_data_epoch(:,k,:) = cell2mat(trial_bip_elec_data);
        end 
        trial_data_epoch(:,k,:) = cell2mat(trial_elec_data);

        fprintf(1,'.');

    end

end