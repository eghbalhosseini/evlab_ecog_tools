function extract_significant_channel(obj,args)
arguments
    obj ecog_data_v2
    args.baseTime = [-0.5 0];
    args.epochTime = [0 0.5];
    args.numPerm = 10000;
    args.p_val = 0.05;
end
    baseTime = args.baseTime;
    epochTime = args.epochTime;

    [baseData,baseData_bip] = obj.extract_trial_epochs(epoch_tw = baseTime);
    [epochData,epochData_bip] = obj.extract_trial_epochs(epoch_tw = epochTime);

    basePower = mean(baseData.^2,3);
    epochPower = mean(epochData.^2,3);

    fprintf(1, '\n>> Extracting significant channels for unipolar high gamma envelope\n');
    fprintf(1,'[');
    [NumTrials, goodtrials] = remove_bad_trials(epochData);
    for iChan = 1:size(basePower,1)
        pSig.pChan(iChan) = permtest(epochPower(iChan,goodtrials(iChan,:)),basePower(iChan,goodtrials(iChan,:)),args.numPerm);
         
        fprintf(1,'.');
    end
    pSig.h_fdr_05 = fdr_bh(pSig.pChan,0.05);
    pSig.h_fdr_01 = fdr_bh(pSig.pChan,0.01);
     fprintf(1,'] done\n');

     if(~isempty(baseData_bip))
        basePower = mean(baseData_bip.^2,3);
        epochPower = mean(epochData_bip.^2,3);
    
        fprintf(1, '\n>> Extracting significant channels for bipolar high gamma envelope\n');
        fprintf(1,'[');
        [NumTrials, goodtrials] = remove_bad_trials(epochData);
        for iChan = 1:size(basePower,1)
            pSig.pChan_bip(iChan) = permtest(epochPower(iChan,goodtrials(iChan,:)),basePower(iChan,goodtrials(iChan,:)),args.numPerm);
            fprintf(1,'.');
        end
         pSig.h_bip_fdr_05 = fdr_bh(pSig.pChan_bip,0.05);
         pSig.h_bip_fdr_01 = fdr_bh(pSig.pChan_bip,0.01);
         fprintf(1,'] done\n');
     end
     obj.stats.sig_hg_channel = pSig;
     
end


function p = permtest(sample1, sample2, numperm)
% permtest - Perform one-sided permutation test to compare the means of two samples.
%
% Syntax: p = permtest(sample1, sample2, numperm)
%
% Inputs:
%   sample1     - First sample data (1 x n1) array
%   sample2     - Second sample data (1 x n2) array
%   numperm     - Number of permutations to perform
%
% Outputs:
%   p           - p-value indicating the significance of the difference between the means
%
% Example:
%   sample1 = [1, 2, 3, 4, 5]; % Example first sample
%   sample2 = [6, 7, 8, 9, 10]; % Example second sample
%   p = permtest(sample1, sample2, 1000); % Perform one-sided permutation test with 1000 permutations

samples = [sample1 sample2]; % Combine the two samples
samplediff = mean(sample1) - mean(sample2); % Calculate the difference between the means of the samples
sampdiffshuff = zeros(1, numperm); % Initialize an array to store shuffled sample differences

for n = 1:numperm
    sampshuff = samples(randperm(length(samples))); % Shuffle the combined samples
    sampdiffshuff(n) = mean(sampshuff(1:length(sampshuff)/2)) - mean(sampshuff(length(sampshuff)/2+1:end)); % Calculate the difference between means for the shuffled samples
end

p = length(find(sampdiffshuff > samplediff)) / numperm; % Calculate the p-value as the proportion of shuffled sample differences greater than the observed difference

end
