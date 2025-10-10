% Filename: pLanA_language_analysis.m
% Description: Complete MATLAB code for Bayesian logistic regression
%              and information‐theoretic metrics on pLanA vs. language selectivity.

%% 0. Load Data
% Assumes workspace variables: lana (vector), isLang (logical vector), coords (Nx3 matrix)
% Example: load('results.mat')  % contains lana, isLang, coords

%% 1. Bayesian Logistic Regression via Metropolis–Hastings

% 1.1 Define log-posterior function
logPost = @(b) sum(isLang.*(b(1)+b(2)*lana) - log(1+exp(b(1)+b(2)*lana))) ...
               + log(normpdf(b(1),0,10)) + log(normpdf(b(2),0,5));

% 1.2 MCMC settings
nsamples      = 12000;
burnin        = 2000;
proposalSigma = diag([1,0.5]);
init          = [0;0];

% 1.3 Metropolis–Hastings sampler
% Corrected Metropolis–Hastings sampler call:
samples = mhsample(init, nsamples, 'logpdf', logPost, ...
    'proprnd', @(x) mvnrnd(x, proposalSigma), ...
    'proppdf', @(y,x) mvnpdf(y, x, proposalSigma), ...
    'burnin', burnin);

% Discard burn-in
beta0 = samples(burnin+1:end,1);
beta1 = samples(burnin+1:end,2);

% 1.4 Posterior summaries
mean_beta1 = mean(beta1);
credInt     = prctile(beta1, [2.5,97.5]);
P_pos       = mean(beta1>0);
fprintf('Bayesian slope: mean=%.3f, 95%% CI=[%.3f,%.3f], P(>0)=%.2f\n', ...
        mean_beta1, credInt(1), credInt(2), P_pos);

% 1.5 Posterior predictive curve
xx    = linspace(min(lana), max(lana), 100)';
ppmat = sigmoid(beta0 + beta1.*xx');
meanCurve = mean(ppmat,1);
stdCurve  = std(ppmat,[],1);

figure; hold on;
fill([xx, fliplr(xx)]', [meanCurve+stdCurve, fliplr(meanCurve-stdCurve)]', ...
     [0.8,0.8,1], 'EdgeColor','none');
plot(xx, meanCurve, 'b-', 'LineWidth',1.5);
xlabel('pLanA'); ylabel('P(isLang=1)');
title('Bayesian Posterior Predictive Curve');
hold off;

%% 2. Information‐Theoretic Metrics

% 2.1 Mutual Information Estimation
NBINS = 1000;
edges = linspace(min(lana), max(lana), NBINS+1);
[~,~,xBin] = histcounts(lana, edges);

jointCounts = accumarray([xBin, double(isLang)+1], 1, [NBINS,2]);
pXY = jointCounts / sum(jointCounts(:));
pX  = sum(pXY,2);
pY  = sum(pXY,1);

MI = 0;
for i = 1:NBINS
    for j = 1:2
        if pXY(i,j)>0
            MI = MI + pXY(i,j) * log2(pXY(i,j)/(pX(i)*pY(j)));
        end
    end
end
fprintf('Mutual Information = %.3f bits\n', MI);
%% 2.3 Information-Gain Curve
ths   = linspace(min(lana), max(lana), 100);
H_Y   = -sum(pY.*log2(pY));
infoG = zeros(size(ths));

for k = 1:numel(ths)
    yhat = lana > ths(k);
    H_YgX = 0;
    for v = [0,1]
        mask = (yhat==v);
        if any(mask)
            p_cond = mean(isLang(mask));
            H_YgX = H_YgX + sum(mask)/numel(lana) * ...
                    (-p_cond*log2(p_cond) - (1-p_cond)*log2(1-p_cond));
        end
    end
    infoG(k) = H_Y - H_YgX;
end

figure;
plot(ths, infoG, 'r-', 'LineWidth',1.5);
xlabel('Threshold on pLanA'); ylabel('Information Gain (bits)');
title('Information-Gain vs. Threshold');
%%
% 2.3 Information‐Gain Curve with Balanced Sampling Across Classes

% Precompute thresholds and overall entropy H_Y
ths = linspace(min(lana), max(lana), 50);
pY  = [mean(~isLang), mean(isLang)];  % class marginals
H_Y = -sum(pY .* log2(pY));

% Identify language vs non‐language indices
idx_lang    = find(isLang);
idx_nonlang = find(~isLang);
nSample     = min(numel(idx_lang), numel(idx_nonlang));

% Number of balanced samplings
nSampling = 1000;

% Storage for information‐gain curves
infoG_samps = zeros(nSampling, numel(ths));

% Balanced sampling loop
for s = 1:nSampling
    % Draw balanced sample
    samp_lang    = randsample(idx_lang,    nSample);
    samp_nonlang = randsample(idx_nonlang, nSample);
    samp_idx     = [samp_lang; samp_nonlang];

    % Subsample data
    lana_s   = lana(samp_idx);
    isLang_s = isLang(samp_idx);

    % Compute information‐gain for this sample
    infoG = zeros(1, numel(ths));
    for k = 1:numel(ths)
        yhat = lana_s > ths(k);
        H_YgX = 0;
        for v = [0,1]
            mask = (yhat == v);
            if any(mask)
                p_cond = mean(isLang_s(mask));
                H_YgX = H_YgX + sum(mask)/numel(lana_s) * ...
                    (-p_cond*log2(p_cond) - (1-p_cond)*log2(1-p_cond));
            end
        end
        infoG(k) = H_Y - H_YgX;
    end
    infoG_samps(s, :) = infoG;
end

% Compute mean and standard deviation across samplings
meanInfoG = mean(infoG_samps, 1);
stdInfoG  = std( infoG_samps,  0, 1);

% Plot mean curve with shaded ±1 SD
figure; hold on;
fill([ths, fliplr(ths)], ...
     [meanInfoG+stdInfoG, fliplr(meanInfoG-stdInfoG)], ...
     [1, 0.8, 0.8], 'EdgeColor','none');
plot(ths, meanInfoG, 'r-', 'LineWidth',1.5);
xlabel('Threshold on pLanA');
ylabel('Information Gain (bits)');
title('Information Gain vs. Threshold (Balanced Sampling)');
hold off;


%% Helper: Sigmoid function
function y = sigmoid(x)
    y = 1 ./ (1 + exp(-x));
end
%%

% 1. Fit logistic regression on full dataset
tbl_all  = table(lana, isLang, 'VariableNames', {'pLanA','isLang'});
mdl_all  = fitglm(tbl_all, 'isLang~pLanA', 'Distribution','binomial');
p_pred   = predict(mdl_all, tbl_all);  % predicted P(isLang=1|pLanA)

% 2. Compute marginal entropy H(Y)
pY       = [mean(~isLang), mean(isLang)];
H_Y      = -sum(pY .* log2(pY));

% 3. Define probability‐based thresholds
ths      = linspace(0,1,100);  % since p_pred ∈ [0,1]
infoG    = zeros(size(ths));

% 4. Loop over thresholds to compute info gain
for k = 1:numel(ths)
    % Partition by predicted probability
    groupHigh = (p_pred > ths(k));   % group where model says "likely Lang"
    groupLow  = ~groupHigh;          % group where model says "likely Nonlang"

    H_YgX     = 0;
    for mask = [groupLow, groupHigh]
        for v = [0,1]
            idx = mask & (isLang == v);
            if any(mask)
                p_cond = sum(isLang(mask)==1) / sum(mask);
                H_group = -p_cond*log2(p_cond) - (1-p_cond)*log2(1-p_cond);
            else
                H_group = 0;
            end
        end
        H_YgX = H_YgX + sum(mask)/numel(isLang) * H_group;
    end

    infoG(k) = H_Y - H_YgX;
end

% 5. Plot information‐gain curve
figure;
plot(ths, infoG, 'b-', 'LineWidth',1.5);
xlabel('Threshold on predicted P(isLang=1)');
ylabel('Information Gain (bits)');
title('Information Gain vs. Model‐Predicted Probability Threshold');
%%
% Assumes: lana (Nx1), isLang (Nx1 logical), mdl_all (fitted GLM)

% 1. Compute model predictions and marginals
p_pred = predict(mdl_all, table(lana, isLang, 'VariableNames',{'pLanA','isLang'}));
pY1    = mean(isLang);
pY0    = 1 - pY1;

% 2. Compute pointwise IG for each electrode
IG = zeros(size(lana));
for i = 1:numel(lana)
    if isLang(i)
        IG(i) = log2(p_pred(i) / pY1);
    else
        IG(i) = log2((1 - p_pred(i)) / pY0);
    end
end

% 3. Bin electrodes by raw pLanA values
NBINS = 100;
edges = linspace(min(lana), max(lana), NBINS+1);
[~,~,binIdx] = histcounts(lana, edges);

% 4. Compute average IG per pLanA bin (using the bin center)
binCenters = (edges(1:end-1) + edges(2:end)) / 2;
meanIG     = nan(NBINS,1);
stdIG      = nan(NBINS,1);
for b = 1:NBINS
    members = (binIdx == b);
    if any(members)
        meanIG(b) = mean(IG(members));
        stdIG(b)  = std( IG(members));
    end
end

% 5. Plot mean ±1 SD IG versus pLanA bin centers
figure; hold on;
fill([binCenters, fliplr(binCenters)], ...
     [meanIG'+stdIG', fliplr(meanIG'-stdIG')], ...
     [0.9,0.9,1], 'EdgeColor','none');
plot(binCenters, meanIG, 'b-', 'LineWidth',1.5);
xlabel('Raw pLanA');
ylabel('Average Information Gain (bits)');
title('Average IG vs. pLanA (Model‐Predicted Probabilities)');
hold off;
