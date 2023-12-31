function [outputArg1,outputArg2] = global_ave(inputArg1,inputArg2)
%GLOBAL_AVE Summary of this function goes here
%   Detailed explanation goes here
% %% do global mean removal (not sure if we should do this!)
% if ops.do_globelMean
%     overall_mean=mean(signal(:,param.ecog_channels),2);
%     %overall_mean=mean(signal,2);
%     %figure;plot(overall_mean);
%     signal_glob=signal-repmat(overall_mean,1,size(signal,2));
%     signal=signal_glob;
%     clear signal_glob;
%     param.preproc_step=[param.preproc_step;'global_mean_removal'];
% end

end

