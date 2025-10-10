%% analysis_pipeline.m
% 1. Configuration
base         = '/Volumes/disk/nese/MGH_ECoG_Langloc/';
organizedDir = '/Users/dsuseendar/data/electrode_data';  % output_dir from organize_electrode_data
crunchDir    = fullfile(base, 'crunched');
lanANii      = fullfile(base, 'language_atlas', 'LanA', 'SPM', 'LanA_n806.nii');

% 2. List crunched files
files = dir(fullfile(crunchDir, '*_obj_*.mat'));
files = files(~startsWith({files.name}, '.'));

% 3. Initialize results table
results = table([], [], [], [], [], [], ...
    'VariableNames', {'R','A','S','subject','isSig','pLanA'});

%% 4. Loop: load each file and organized electrode data
% Parallel processing of subject files with test_s_vs_n


nFiles = numel(files);

% Preallocate cell array for results tables
resultsCell = cell(nFiles,1);

parfor iFile = 1:nFiles
    f = files(iFile);
    % Load subject object
    D = load(fullfile(f.folder, f.name));
    fnames = fieldnames(D);
    if isempty(fnames)
        warning('File %s contains no variables.', f.name);
        resultsCell{iFile} = table(); 
        continue;
    end
    od = D.(fnames{1});
    subj = od.subject_id;
    fprintf('Processing %s…\n', subj);
    
    % Run S vs N test
    od.test_s_vs_n();
    isSig_crunched = od.s_vs_n_sig.elec_data_dec{1};
    ch_labels = cellstr(od.elec_ch_label);
    
    % Load combined electrode data
    matFile = fullfile(organizedDir, subj, 'electrodes_combined.mat');
    if ~exist(matFile, 'file')
        warning('Missing organized data for %s, skipping.', subj);
        resultsCell{iFile} = table(); 
        continue;
    end
    C = load(matFile, 'combined');
    combined_labels = cellstr(C.combined.labels);
    coords = C.combined.coords;
    
    % Find common channels
    [chan_common, idx_cr, idx_cb] = intersect(ch_labels, combined_labels, 'stable');
    if isempty(chan_common)
        warning('No channel intersection for %s.', subj);
        resultsCell{iFile} = table(); 
        continue;
    end
    
    % Extract aligned data
    isSig_sub = isSig_crunched(idx_cr);
    coords_sub = coords(idx_cb, :);
    n = numel(idx_cb);
    
    % Build table
    T = table(coords_sub(:,1), coords_sub(:,2), coords_sub(:,3), ...
              repmat({subj}, n, 1), isSig_sub, NaN(n,1), ...
              'VariableNames', {'R','A','S','subject_id','isSig','placeholder'});
    
    resultsCell{iFile} = T;
end

% Concatenate non-empty tables
results = vertcat(resultsCell{:});


%% 5. Interpolate LanA probabilities from NIfTI
V      = spm_vol(lanANii);
[Y,~]  = spm_read_vols(V);
results.pLanA(:) = 0;  % default
for i = 1:height(results)
    xyz = [results.R(i), results.A(i), results.S(i), 1]';
    ijk = round(V.mat \ xyz);
    if all(ijk(1:3) >= 1) && ...
       ijk(1) <= size(Y,1) && ijk(2) <= size(Y,2) && ijk(3) <= size(Y,3)
        results.pLanA(i) = Y(ijk(1), ijk(2), ijk(3));
    end
end

%% 6. Statistical analyses with spin permutation
lana   = results.pLanA_gauss;
isLang = logical(results.isSig);
coords = [results.R, results.A, results.S];  % electrode coordinates



% Identify indices of each class
idx_lang    = find(isLang);
idx_nonlang = find(~isLang);
nSample     = min(numel(idx_lang), numel(idx_nonlang));
nSampling     = 1000;
coef_spin = zeros(nSampling,1);
coef_true = zeros(nSampling,1);

xx = linspace(min(lana_probs), max(lana_probs), 100)';

% Preallocate storage
true_curves  = nan(nSampling, numel(xx));
shuf_curves  = nan(nSampling, numel(xx));

b = ProgressBar([], ...
    'IsParallel', true, ...
    'Title', 'Permutation testing' ...
    );

% ALWAYS CALL THE SETUP() METHOD FIRST!!!
b.setup([], [], []);

parfor k = 1:nSampling
    % Draw a balanced sample for the "true" fit
    s_lang = randsample(idx_lang,    nSample);
    s_non  = randsample(idx_nonlang, nSample);
    s_idx  = [s_lang; s_non];
    
    % Subsample data
    lana_bal   = lana(s_idx);
    isLang_bal = isLang(s_idx);
    coords_bal = coords(s_idx, :);
    
    % 6a. Logistic regression on balanced sample
    tbl_true   = table(lana_bal, isLang_bal, 'VariableNames', {'pLanA','isLang'});
    mdl_true   = fitglm(tbl_true, 'isLang~pLanA', 'Distribution','binomial');
    coef_true(k)  = mdl_true.Coefficients.Estimate('pLanA');
    true_curves(k, :) = predict(mdl_true, table(xx,'VariableNames',{'pLanA'}));

    
    
    % 6b. Spin permutation test on balanced sample
    spin_maps = generate_spin_permuted_maps(coords_bal, lana_bal, 1);   
    lana_spin     = spin_maps;
    mdl_perm      = fitglm(table(lana_spin, isLang_bal, ...
                          'VariableNames', {'pLanA','isLang'}), ...
                          'isLang~pLanA', 'Distribution','binomial');
    coef_spin(k)  = mdl_perm.Coefficients.Estimate('pLanA');
    shuf_curves(k, :) = predict(mdl_perm, table(xx,'VariableNames',{'pLanA'}));
    updateParallel()
end
b.release();
p_coef_spin = mean((coef_spin) >= (coef_true));

% %% 6c. ROC/AUC analysis on balanced sample with spin permutation
% [~,~,~,auc_true] = perfcurve(isLang_bal, lana_bal, true);
% auc_spin = zeros(nSampling,1);
% 
% for k = 1:nSampling
%     [~,~,~,auc_spin(k)] = perfcurve(isLang_bal, spin_maps(:,k), true);
% end
% p_auc_spin = mean(auc_spin >= auc_true);
% 
% %% 6d. Standard permutation test on balanced sample
% coef_perm = zeros(nSampling,1);
% auc_perm  = zeros(nSampling,1);
% 
% for k = 1:nSampling
%     % Shuffle labels within balanced sample
%     y_shuf       = isLang_bal(randperm(nSample));
% 
%     % Logistic slope
%     mdl_p        = fitglm(table(lana_bal, y_shuf, 'VariableNames', {'pLanA','isLang'}), ...
%                           'isLang~pLanA', 'Distribution','binomial');
%     coef_perm(k) = mdl_p.Coefficients.Estimate('pLanA');
% 
%     % ROC/AUC
%     [~,~,~,auc_perm(k)] = perfcurve(y_shuf, lana_bal, true);
% end
% p_coef_perm = mean(abs(coef_perm) >= abs(coef_true));
% p_auc_perm  = mean(auc_perm >= auc_true);

%% 6e. Display results
% fprintf('True slope:        %.4f\n', coef_true);
% fprintf('Standard perm p:   coef=%.4f, auc=%.4f\n', p_coef_perm, p_auc_perm);
% fprintf('Spin perm p:       coef=%.4f, auc=%.4f\n', p_coef_spin, p_auc_spin);
% fprintf('True AUC:          %.4f\n', auc_true);

% %% 7. Package and save analysis output
% analysis.results      = results;
% analysis.balance_idx  = s_idx;
% analysis.mdl_true     = mdl_true;
% analysis.coef_true    = coef_true;
% analysis.coef_perm    = coef_perm;
% analysis.p_coef_perm  = p_coef_perm;
% analysis.coef_spin    = coef_spin;
% analysis.p_coef_spin  = p_coef_spin;
% analysis.auc_true     = auc_true;
% analysis.auc_perm     = auc_perm;
% analysis.p_auc_perm   = p_auc_perm;
% analysis.auc_spin     = auc_spin;
% analysis.p_auc_spin   = p_auc_spin;
% 
% save('analysis_output.mat','analysis');
% fprintf('\nBalanced, spin, and standard permutation analysis complete.\n');
%%

V = spm_vol(lanANii);

% 2. Optional: load pial surface
%surf = load('pial_surf.mat');           % surf.vertices (Nv×3)

% 3. Apply each method
p1 = assign_pLanA_methods(results, V, 'trilinear');
%p2 = assign_pLanA_methods(results, V, 'surface', surf);
p3 = assign_pLanA_methods(results, V, 'radial', 10);
p4 = assign_pLanA_methods(results, V, 'distweight');
p5 = assign_pLanA_methods(results, V, 'gaussian', [2 2 2]);

% 4. Add back to table and save
results.pLanA_trilin   = p1;
%results.pLanA_surface  = p2;
results.pLanA_radial   = p3;
results.pLanA_distw    = p4;
results.pLanA_gauss    = p5;
%save('results_with_pLanA.mat','results');

% Assuming results.isLang is a logical vector of length N
N = height(results);
langIdx = find(results.isSig);

figure; hold on;
plot(results.pLanA_trilin, '.-', 'DisplayName','Trilinear');
plot(results.pLanA_radial,  '.-', 'DisplayName','Radial');
plot(results.pLanA_distw,   '.-', 'DisplayName','DistWeight');
plot(results.pLanA_gauss,   '.-', 'DisplayName','Gaussian');

% Add vertical bars at language-selective indices
yl = ylim;  % get current y-axis limits
for k = 1:numel(langIdx)
    x = langIdx(k);
    plot([x x], yl, 'k--', 'LineWidth', 0.5);
end


xlabel('Electrode index');
ylabel('pLanA');
title('Comparison of pLanA Assignment Methods with Language-Selective Markers');
hold off;
