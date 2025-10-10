% Filename: compare_pLanA_methods_modelFits.m

%% 1. Setup and Data Loading
load('results.mat','results');            % must contain R,A,S,isSig
V = spm_vol(lanANii);

methods   = {'trilinear','radial','distweight','gaussian'};
nMethods  = numel(methods);

%% 2. Assign pLanA
results.pLanA_trilin = assign_pLanA_methods(results, V, 'trilinear');
results.pLanA_radial = assign_pLanA_methods(results, V, 'radial', 10);
results.pLanA_distw  = assign_pLanA_methods(results, V, 'distweight');
results.pLanA_gauss  = assign_pLanA_methods(results, V, 'gaussian', [2 2 2]);

%% 3. Prepare for permutation testing
isLang      = logical(results.isSig);
coords      = [results.R, results.A, results.S];
idx_lang    = find(isLang);
idx_nonlang = find(~isLang);
nSample     = min(numel(idx_lang), numel(idx_nonlang));
nSampling   = 500;          
xx          = linspace(0,1,100)';

% Preallocate 3D arrays: [method × sampling × xx]
true_curves = nan(nMethods, nSampling, numel(xx));
shuf_curves = nan(nMethods, nSampling, numel(xx));

%% 4. Permutation testing

b = ProgressBar([], ...
    'IsParallel', true, ...
    'Title', 'Permutation testing' ...
    );

% ALWAYS CALL THE SETUP() METHOD FIRST!!!
b.setup([], [], []);
parfor k = 1:nSampling
    % Balanced subsample
    s_lang = randsample(idx_lang,    nSample);
    s_non  = randsample(idx_nonlang, nSample);
    s_idx  = [s_lang; s_non];
    isBal  = isLang(s_idx);
    coordsBal = coords(s_idx,:);
    
    for m = 1:nMethods
        % Select method-specific pLanA
        switch methods{m}
            case 'trilinear'
                lana = results.pLanA_trilin(s_idx);
            case 'radial'
                lana = results.pLanA_radial(s_idx);
            case 'distweight'
                lana = results.pLanA_distw(s_idx);
            case 'gaussian'
                lana = results.pLanA_gauss(s_idx);
        end
        
        % True logistic fit
        tblT = table(lana, isBal, 'VariableNames', {'pLanA','isLang'});
        mdlT = fitglm(tblT, 'isLang~pLanA', 'Distribution','binomial');
        true_curves(m, k, :) = predict(mdlT, table(xx,'VariableNames',{'pLanA'}));
        
        % Spin‐permutation fit
        spin = generate_spin_permuted_maps(coordsBal, lana, 1);
        tblS = table(spin, isBal, 'VariableNames', {'pLanA','isLang'});
        mdlS = fitglm(tblS, 'isLang~pLanA', 'Distribution','binomial');
        shuf_curves(m, k, :) = predict(mdlS, table(xx,'VariableNames',{'pLanA'}));
    end
    updateParallel()
end
b.release()

%% 5. Plot mean ±1 SD for each method
colors = lines(nMethods);
figure('Color','w','Position',[200 200 800 600]);
hold on;

for m = 1:nMethods
    tc = squeeze(true_curves(m, :, :));  % [nSampling × xx]
    sc = squeeze(shuf_curves(m, :, :));
    mu_tc = mean(tc,1); sd_tc = std(tc,0,1);
    mu_sc = mean(sc,1); sd_sc = std(sc,0,1);
    
    X = [xx; flipud(xx)];
    Yt = [mu_tc+sd_tc, fliplr(mu_tc-sd_tc)];
    Ys = [mu_sc+sd_sc, fliplr(mu_sc-sd_sc)];
    
    %patch(X, Ys, colors(m,:), 'FaceAlpha',0.2, 'EdgeColor','none');
    %patch(X, Yt, colors(m,:), 'FaceAlpha',0.4, 'EdgeColor','none');
   % plot(xx, mu_sc, '--', 'Color',colors(m,:), 'LineWidth',1.5);
    plot(xx, mu_tc, '-',  'Color',colors(m,:), 'LineWidth',2);
end

xlabel('pLanA');
ylabel('P(isLang=1)');
title('True vs. Shuffled Model Fits Across pLanA Methods');
% legendEntries = reshape([strcat(methods,' shuffle'); strcat(methods,' true')],1,[]);
% legend(legendEntries,'Location','NorthWest');
ylim([-0.05 1.05]);
hold off;
