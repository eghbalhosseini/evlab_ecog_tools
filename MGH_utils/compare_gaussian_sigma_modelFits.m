% Filename: compare_gaussian_sigmas_modelFits.m
% Compare model fits (true vs. shuffled) for different Gaussian smoothing FWHM

%% 1. Setup and Data Loading
load('results.mat','results');            % contains R,A,S,isSig
V = spm_vol('lanANii.nii');               % lanA probability volume

%% 2. Define Gaussian kernels to test (FWHM in mm)
sigma_list = {[1 1 1], [2 2 2], [4 4 4], [6 6 6], [10 10 10], [15 15 15]};
nSigma     = numel(sigma_list);

%% 3. Assign pLanA for each sigma
pLanA_gauss = nan(height(results), nSigma);
for j = 1:nSigma
    pLanA_gauss(:,j) = assign_pLanA_methods(results, V, 'gaussian', sigma_list{j});
end

%% 4. Permutation testing parameters
isLang      = logical(results.isSig);
coords      = [results.R, results.A, results.S];
idx_lang    = find(isLang);
idx_nonlang = find(~isLang);
nSample     = min(numel(idx_lang), numel(idx_nonlang));
nSampling   = 300;             % adjust for speed
xx          = linspace(0,1,100)';

% Preallocate arrays: [sigma × sampling × xx]
true_curves = nan(nSigma, nSampling, numel(xx));
shuf_curves = nan(nSigma, nSampling, numel(xx));

%% 5. Permutation loop
parfor k = 1:nSampling
    % Balanced subsample
    s_lang = randsample(idx_lang,    nSample);
    s_non  = randsample(idx_nonlang, nSample);
    s_idx  = [s_lang; s_non];
    isBal  = isLang(s_idx);
    coordsBal = coords(s_idx,:);
    
    for j = 1:nSigma
        lana = pLanA_gauss(s_idx, j);
        
        % True fit
        tblT = table(lana, isBal, 'VariableNames', {'pLanA','isLang'});
        mdlT = fitglm(tblT, 'isLang~pLanA', 'Distribution','binomial');
        true_curves(j, k, :) = predict(mdlT, table(xx,'VariableNames',{'pLanA'}));
        
        % Shuffled fit
        lana_spin = generate_spin_permuted_maps(coordsBal, lana, 1);
        tblS = table(lana_spin, isBal, 'VariableNames', {'pLanA','isLang'});
        mdlS = fitglm(tblS, 'isLang~pLanA', 'Distribution','binomial');
        shuf_curves(j, k, :) = predict(mdlS, table(xx,'VariableNames',{'pLanA'}));
    end
end

%% 6. Plot mean ±1 SD curves for each sigma
colors = jet(nSigma);
figure('Color','w','Position',[100 100 800 600]);
hold on;

for j = 1:nSigma
    figure
    hold on;
    tc = squeeze(true_curves(j,:,:));  % [nSampling × xx]
    sc = squeeze(shuf_curves(j,:,:));
    mu_tc = mean(tc,1); sd_tc = std(tc,0,1);
    mu_sc = mean(sc,1); sd_sc = std(sc,0,1);
    
    X  = [xx; flipud(xx)];
    Yt = [mu_tc+sd_tc, fliplr(mu_tc-sd_tc)];
    Ys = [mu_sc+sd_sc, fliplr(mu_sc-sd_sc)];
    
    % % Shuffled ribbon & mean
    patch(X, Ys, colors(j,:), 'FaceAlpha',0.2, 'EdgeColor','none');
    plot(xx, mu_sc, '--', 'Color',colors(j,:), 'LineWidth',1.5);
    
    %True ribbon & mean
    patch(X, Yt, colors(j,:), 'FaceAlpha',0.4, 'EdgeColor','none');
    plot(xx, mu_tc, '-',  'Color',colors(j,:), 'LineWidth',2);
end

xlabel('pLanA');
ylabel('P(isLang=1)');
title('Gaussian Smoothing: True vs. Shuffled Model Fits for Various FWHM');
legendEntries = reshape([...
    strcat(string(cellfun(@(s) sprintf('%dmm',s(1)), sigma_list)),' shuffle'); ...
    strcat(string(cellfun(@(s) sprintf('%dmm',s(1)), sigma_list)),' true')],1,[]);
legend(legendEntries,'Location','NorthWest');
ylim([-0.05 1.05]);
hold off;
