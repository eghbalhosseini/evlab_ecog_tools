

%% Full s-vs-n Analysis & LanA Overlay for Your Crunched Files
% Modified to work with FreeSurfer surfaces using proper reading functions
% Assumes files named like EM1041_obj_visual.mat, EM1050_obj_auditory.mat, etc.

addpath(genpath('MGH_utils'));
addpath(('/Users/dsuseendar/repos/extra/spm'));

% Add FreeSurfer MATLAB functions to path
if exist(fullfile(getenv('FREESURFER_HOME'), 'matlab'), 'dir')
    addpath(genpath(fullfile(getenv('FREESURFER_HOME'), 'matlab')));
else
    warning('FreeSurfer MATLAB functions not found. Make sure FREESURFER_HOME is set correctly.');
end

%% 1. Paths & File Listing
base        = '/Volumes/disk/nese/MGH_ECoG_Langloc/';
crunchDir   = fullfile(base,'crunched');
elecDir     = fullfile(base,'brainstom_localizations_cat12_fedorenko_nov_2024');
labelDir    = fullfile(base,'brainstorm_localizations/');

% FreeSurfer surfaces - using proper FreeSurfer path structure
fsSubjectsDir = getenv('SUBJECTS_DIR');
if isempty(fsSubjectsDir)
    fsSubjectsDir = fullfile(getenv('FREESURFER_HOME'), 'subjects');
end

% Use fsaverage or specify your template subject
templateSubj = 'fsaverage';
pialLHfile   = fullfile(fsSubjectsDir, templateSubj, 'surf', 'lh.pial');
pialRHfile   = fullfile(fsSubjectsDir, templateSubj, 'surf', 'rh.pial');

% Language atlas
lanANii     = fullfile(base,'language_atlas','LanA','SPM','LanA_n806.nii');

%% 2. List crunched files
files = dir(fullfile(crunchDir,'*_obj_*.mat'));
files = files(~startsWith({files.name},'.'));

%% 3. Initialize results table
results = table([],[],[],[],[],[],'VariableNames', ...
    {'R','A','S','subject','isSig','pLanA'});

%% 4. Loop: load each file, run test_s_vs_n, collect electrodes
for f = files'
    D  = load(fullfile(f.folder,f.name));
    fd = fieldnames(D);
    od = D.(fd{1});           
    
    % Load channel anatomical labels CSV (as before)
    subj = od.subject_id;
    csvFile = fullfile(labelDir,sprintf('sub-%s_channel_anatomical_labels.csv',subj));
    if(isfile(csvFile))
        anatTbl = readtable(csvFile);
        anatTbl.Properties.VariableNames{1} = 'label';
        od.elec_ch_pos_anat   = anatTbl(:, {'label','mni_linear_x','mni_linear_y','mni_linear_z'});
        od.elec_ch_HCP_label  = anatTbl.HCPMMP1_label_1;
        od.elec_ch_HCP_weight = anatTbl.HCPMMP1_weight_1;
    end
    
    % Load channel RAS table from TSV
    tsvPattern = fullfile(elecDir, sprintf('sub-%s*_electrodes*.tsv',subj));
    tbls = dir(tsvPattern); tbls = tbls(~startsWith({tbls.name},'.'));
    if(~isempty(tbls))
        ch_RAS_tbl = readtable(fullfile(tbls(1).folder, tbls(1).name), 'Delimiter', '\t', 'FileType', 'text');
        ch_RAS_tbl.Properties.VariableNames(ismember(ch_RAS_tbl.Properties.VariableNames,'x')) = {'R'};
        ch_RAS_tbl.Properties.VariableNames(ismember(ch_RAS_tbl.Properties.VariableNames,'y')) = {'A'};
        ch_RAS_tbl.Properties.VariableNames(ismember(ch_RAS_tbl.Properties.VariableNames,'z')) = {'S'};
    end
    
    % Use sn.s_vs_n_sig results directly (assuming already computed)
    isSig = od.s_vs_n_sig.elec_data_zs_dec{1};
    
    % Use od.elec_ch_pos_mni if nonempty; otherwise fallback on ch_RAS_tbl
    RAS_cell = od.elec_ch_pos_mni;
    if isempty(RAS_cell)
        % ch_RAS_tbl has numeric R,A,S columns, convert to matrix
        RAS = [ch_RAS_tbl.R, ch_RAS_tbl.A, ch_RAS_tbl.S];
    else
        % convert from cell array of coordinates
        RAS = cell2mat(RAS_cell);
        if size(RAS,2) > 3
            RAS = RAS(:,1:3);
        end
    end
    
    % Append data to results table
    n = size(RAS,1);
    try
        T = table(RAS(:,1), RAS(:,2), RAS(:,3), ...
                  repmat({od.subject_id},n,1), isSig, NaN(n,1), ...
                  'VariableNames', results.Properties.VariableNames);
        results = [results; T];
    catch
        disp('Cannot append table')
        continue
    end
end

%% 5. Load LanA atlas & interpolate probability
V = spm_vol(lanANii);
[Y,~] = spm_read_vols(V);
for i = 1:height(results)
    xyz = [results.R(i),results.A(i),results.S(i),1]';
    ijk = round(V.mat\xyz);
    results.pLanA(i) = Y(ijk(1),ijk(2),ijk(3));
end

%% 6. Load FreeSurfer pial surfaces using read_surf
% Check if surfaces exist
if ~exist(pialLHfile, 'file')
    error('Left hemisphere pial surface not found: %s', pialLHfile);
end
if ~exist(pialRHfile, 'file')
    error('Right hemisphere pial surface not found: %s', pialRHfile);
end

% Read surfaces using FreeSurfer's read_surf function
try
    [lh_vertices, lh_faces] = read_surf(pialLHfile);
    [rh_vertices, rh_faces] = read_surf(pialRHfile);
    
    % Create surface structures
    lh.vertices = lh_vertices;
    lh.faces = lh_faces;
    rh.vertices = rh_vertices;
    rh.faces = rh_faces;
    
    fprintf('Successfully loaded FreeSurfer surfaces:\n');
    fprintf('  LH: %d vertices, %d faces\n', size(lh_vertices,1), size(lh_faces,1));
    fprintf('  RH: %d vertices, %d faces\n', size(rh_vertices,1), size(rh_faces,1));
    
catch ME
    error('Failed to load FreeSurfer surfaces: %s', ME.message);
end

%% 7. Plot hemispheres with proper FreeSurfer surface handling
%figure('Color','w','Position',[100 100 1600 600]);

for h = 1:2
    figure('Color','w');
    ax = gca
    hold(ax,'on');
    
    % Select mesh & view
    if h==1
        mesh = lh;
        view(ax,[-90 0]);
        title(ax,'Left Hemisphere');
        idx = results.R < 0;
    else
        mesh = rh;
        view(ax,[90 0]);
        title(ax,'Right Hemisphere'); 
        idx = results.R > 0;
    end
    
    % Draw cortical surface using trisurf
    % Note: FreeSurfer faces are 0-indexed, MATLAB needs 1-indexed
    faces_matlab = mesh.faces + 1;  % Convert to 1-indexed
    
    trisurf(faces_matlab, ...
            mesh.vertices(:,1), mesh.vertices(:,2), mesh.vertices(:,3), ...
            'FaceColor',[.8 .8 .8], 'EdgeColor','none', 'FaceAlpha',.5, ...
            'Parent',ax);
    
    % Set lighting and appearance
    lighting(ax,'gouraud');
    camlight(ax,'headlight');
    axis(ax,'equal','off');
    
    % Plot electrodes colored by LanA probability
    if sum(idx) > 0
        scatter3(results.R(idx), results.A(idx), results.S(idx), ...
                10, results.pLanA(idx), 'filled', 'Parent',ax);
    end
    
    % Set colormap and color limits
    lmap = lanamap(256);
    colormap(ax,lmap);
    caxis(ax,[0 0.8]);
    
    % Overlay significant electrodes
    sigI = idx & results.isSig;
    if sum(sigI) > 0
        scatter3(results.R(sigI), results.A(sigI), results.S(sigI), ...
                20, 'o', 'MarkerEdgeColor','cyan','LineWidth',1.5, ...
                'Parent',ax);
    end

    % Add shared colorbar
    cb = colorbar('Position',[0.92 0.3 0.02 0.4]);
    cb.Label.String = 'LanA probability';

end

% Add shared colorbar
cb = colorbar('Position',[0.92 0.3 0.02 0.4]);
cb.Label.String = 'LanA probability';

%% 8. Quantify LanA distributions
lana_probs = results.pLanA;
isLang = logical(results.isSig);

%% 8.1. Overlapping Histogram
edges = 0:0.01:0.9;
figure('Position',[100 100 1200 500]);

subplot(1,2,1);
h1 = histogram(lana_probs(~isLang), edges, 'Normalization','probability', ...
               'FaceColor', [0.3 0.5 1], 'FaceAlpha', 0.7); 
hold on;
h2 = histogram(lana_probs(isLang), edges, 'Normalization','probability', ...
               'FaceColor', [1 0.3 0.3], 'FaceAlpha', 0.7);

meanLang = mean(lana_probs(isLang));
meanNonLang = mean(lana_probs(~isLang));
ylims = ylim;
plot([meanLang meanLang], ylims, 'r--', 'LineWidth', 2);
plot([meanNonLang meanNonLang], ylims, 'b--', 'LineWidth', 2);

legend({'non-language', 'language', ...
        sprintf('lang, mean P(x): %.2f',meanLang), ...
        sprintf('non-lang, mean P(x): %.2f',meanNonLang)});
xlabel('LaNA probability');
ylabel('frequency');
title('distribution of LaNA probabilites for language and non-language electrodes');

%% 8.2. Proportion Plot by Bin
edges = 0:0.05:0.8;
subplot(1,2,2);

[countLang,~]    = histcounts(lana_probs(isLang), edges);
[countNonLang,~] = histcounts(lana_probs(~isLang), edges);
countsAll = countLang + countNonLang;
proportionLang = countLang ./ countsAll;
proportionNonLang = countNonLang ./ countsAll;

barWidth = 1;
b = bar(edges(1:end-1)+diff(edges)/2, [proportionLang', proportionNonLang'], ...
        'stacked','BarWidth',barWidth);
b(1).FaceColor = [1 0.3 0.3];
b(2).FaceColor = [0.3 0.5 1];

% Annotate counts
hold on
for i = 1:length(countLang)
    if countsAll(i) > 0
        text(edges(i)+0.02, proportionLang(i)/2, num2str(countLang(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment','middle',...
            'Color',[0 0 0],'FontSize',11);
        text(edges(i)+0.02, proportionLang(i)+proportionNonLang(i)/2, num2str(countNonLang(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment','middle',...
            'Color',[0 0 0],'FontSize',11);
    end
end

xlabel('LaNA proportion');
ylabel('frequency');
legend({'non-lang','lang'});
title('proportion of language and non-language electrodes');
ylim([0 1]);

% Tighten layout
set(gcf,'Color','w');

%% 9. Additional Statistical Analysis (keeping your existing code)
% [Rest of your statistical analysis code remains the same...]

%% Helper function to check FreeSurfer installation
function checkFreeSurferSetup()
    % Check if FreeSurfer is properly set up
    fs_home = getenv('FREESURFER_HOME');
    if isempty(fs_home)
        warning('FREESURFER_HOME environment variable not set');
        return;
    end
    
    % Check if read_surf function is available
    if exist('read_surf', 'file') ~= 2
        warning('read_surf function not found. Make sure FreeSurfer MATLAB functions are in path');
        fprintf('Try adding: addpath(genpath(''%s''))\n', fullfile(fs_home, 'matlab'));
    end
    
    % Check subjects directory
    subjects_dir = getenv('SUBJECTS_DIR');
    if isempty(subjects_dir)
        subjects_dir = fullfile(fs_home, 'subjects');
    end
    
    if ~exist(subjects_dir, 'dir')
        warning('SUBJECTS_DIR not found: %s', subjects_dir);
    end
end
