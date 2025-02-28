function extract_normalization_metrics(obj,epoch_args)

    arguments
        obj ecog_data_v2
        epoch_args.baseTimeRange = [-0.5 0];
        epoch_args.timepad = 0.5;
        epoch_args.key='fix';
    end

    timePad = epoch_args.timepad;
    baseTimeRange = epoch_args.baseTimeRange;
    baseTimeExtract = [baseTimeRange(1)-timePad baseTimeRange(2)+timePad];

    [baseData,baseData_bip] = obj.extract_trial_epochs(epoch_tw = baseTimeExtract,key=epoch_args.key);
    [~,goodtrials] = remove_bad_trials(baseData);
    goodTrialsCommon = extractCommonTrials(goodtrials);
    fprintf(1, '\n>> Extracting normalization metrics for unipolar high gamma envelope \n');
    normFactor = zeros(size(baseData, 1), 2);
    fprintf(1,'[');
    for iChan = 1:size(baseData, 1)
        normFactor(iChan, :) = [mean2(squeeze(baseData(iChan, goodTrialsCommon, :))), std2(squeeze(baseData(iChan, goodTrialsCommon, :)))];
        fprintf(1,'.');
    end
    fprintf(1,'] done\n')
    normMetrics.normFactor = normFactor;

    if(~isempty(baseData_bip))
        
        fprintf(1, '\n>> Extracting normalization metrics for bipolar high gamma envelope \n');
        normFactor = zeros(size(baseData_bip, 1), 2);
        [~,goodtrials] = remove_bad_trials(baseData_bip);
        goodTrialsCommon = extractCommonTrials(goodtrials);
        fprintf(1,'[');
        for iChan = 1:size(baseData_bip, 1)
            normFactor(iChan, :) = [mean2(squeeze(baseData_bip(iChan, goodTrialsCommon, :))), std2(squeeze(baseData_bip(iChan, goodTrialsCommon, :)))];
            fprintf(1,'.');
        end
        fprintf(1,'] done\n')
        normMetrics.normFactor_bip = normFactor;
    end

    obj.stats.normMetrics = normMetrics;
  
    


end