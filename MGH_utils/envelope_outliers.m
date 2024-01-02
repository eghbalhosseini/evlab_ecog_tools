function [envelopes_rjct, ops_out] ...
    = envelope_outliers(envelopes, env_sr, varargin)

% Detects outliers in signal envelopes. The difference between the median
% and the Qth percentile (default: 90%) of the envelope distribution for
% each electrode is measured. Envelopes are considered outliers if they
% fall a certain number of units above this value.
% 
% 2016-08-15 - Created, Sam NH
% 
% 2019-01-21 - Cosmetic changes, Sam NH
%
% 2021-11-29 - update to fit with Evlab ecog pipeline 

p=inputParser();
addParameter(p, 'scale', 0);
addParameter(p, 'threshold', 5);
addParameter(p, 'percentile', 0.9);
addParameter(p,'buffer',20);
addParameter(p,'interpMethod','linear');

parse(p, varargin{:});
ops = p.Results;


% total number of channels
n_channels = size(envelopes,2);

% make scale/threshold parameters electrode specific if not already
% -> scale x channel/electrode
if isvector(ops.scale)
    ops.scales = repmat(ops.scale(:), 1, n_channels);
else
    assert(size(ops.scale) == n_channels);
end
if isvector(ops.threshold)
    ops.threshold = repmat(ops.threshold(:), 1, n_channels);
else
    assert(size(ops.threshold) == n_channels);
end

% normalize by 90% of the distribution
s = quantile(envelopes, ops.percentile);
Z = envelopes ./ repmat(s, size(envelopes,1),1);
clear sd;

%% loop through channels
outliers = Z > ops.threshold;
envelopes_rjct=nan*envelopes;
outlier_prcnt=nan*zeros(size(envelopes,2),1);
outliers_idx={};
pbar=ProgressBar(size(envelopes,2), ...
    'Title', 'Progress');
for q=1:n_channels
    outlie=outliers(:,q);
    envl=envelopes(:,q);
    outlier_prcnt(q)=sum(outlie)./numel(outlie)*1e2;
    envl_rejct=envl;
    rise=find(diff(outlie)==1).';
    fall=find(diff(outlie)==-1).';
    % drop events in the first 50 and last 50 samples 
    len_x = length(outlie);
    keep_idx = (rise > 50) & (rise < (len_x - 50));
    rise= rise(keep_idx);

    keep_idx = (fall > 50) & (fall < (len_x - 50));
    fall= fall(keep_idx);

    try 
        assert(numel(rise)==numel(fall));
        process=true;
    catch err 
        pbar.printMessage(sprintf('unequal rise and fall of outliers, ignoring channel %d \n',q));
        process=false;
    end 


    outlie_idx={};
    if numel(rise)>0 & process
    for r=rise
        f=fall(find(fall>r,1,'first'));
        if f-r>10
           pbar.printMessage(sprintf('rise and fall are too far %d, ignoring event in channel %d \n',f-r,q)); 
           continue
        end 
        interp_win=r:f;
        interp_val=envl(interp_win).';
        lookup_window=max([r-ops.buffer,1]):min([(f+ops.buffer),numel(outlie)]);
        main_sig=envl(lookup_window);
        aux=main_sig;
        fixed_sig=main_sig;
        [C,ia,ib]=intersect(lookup_window,interp_win);
        aux(ia)=[];
        intrp_sig=interp1(setdiff(lookup_window,interp_win,'stable'),aux,interp_win,ops.interpMethod);
        fixed_sig(ia)=intrp_sig;
        envl_rejct(lookup_window)=fixed_sig;
        outlie_idx=[outlie_idx,[interp_win;interp_val]];
    end
    
    end
    envelopes_rjct(:,q)=envl_rejct;
    outliers_idx=[outliers_idx;{outlie_idx}];
    pbar(1,[],[]);
end 
pbar.release()
ops_out=ops;
ops_out.outliers_idx=outliers_idx;
ops_out.outlier_prcnt=outlier_prcnt;
    
end