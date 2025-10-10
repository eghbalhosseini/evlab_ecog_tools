% organize_electrode_data_v4.m
% Organize electrode files (CSV) and anatomical labels (CSV/TSV) for FreeSurfer-compatible visualization

function organize_electrode_data_v3(source_dir, output_dir)
    % Validate inputs
    if nargin < 2
        error('Usage: organize_electrode_data_v4(source_dir, output_dir)');
    end
    
    % Create output directory
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % 1) Gather localization files (CSV only) and ignore hidden
    loc_dir = fullfile(source_dir, 'brainstorm_localizations_cat12_fedorenko_nov_2024');
    loc_csv = dir(fullfile(loc_dir, '*.csv'));
    loc_files = loc_csv(~startsWith({loc_csv.name}, '.'));
    
    % 2) Gather anatomical label files (CSV or TSV) in secondary folder and ignore hidden
    label_dir = fullfile(source_dir, 'brainstorm_localizations');
    lbl_csv = dir(fullfile(label_dir, '*.csv'));
    lbl_tsv = dir(fullfile(label_dir, '*.tsv'));
    label_files = [lbl_csv; lbl_tsv];
    label_files = label_files(~startsWith({label_files.name}, '.'));
    
    fprintf('Found %d localization CSV files and %d label files\n', numel(loc_files), numel(label_files));
    
    % 3) Inspect first few localization files
    fprintf('\nInspecting localization file formats...\n');
    for i = 1:min(3, numel(loc_files))
        fpath = fullfile(loc_files(i).folder, loc_files(i).name);
        fprintf('\nFile: %s\n', loc_files(i).name);
        inspect_table_format(fpath);
    end
    
    % 4) Process each localization file
    processed = {};
    failed = {};
    for i = 1:numel(loc_files)
        fn = loc_files(i).name;
        fpath = fullfile(loc_files(i).folder, fn);
        sid = extract_subject_id(fn);
        if isempty(sid)
            warning('Could not extract subject ID from %s', fn);
            failed{end+1} = fn;
            continue;
        end
        % Find matching label file for subject
        lbl = find_label_file(label_files, sid);
        lbl_path = '';
        if ~isempty(lbl)
            lbl_path = fullfile(lbl.folder, lbl.name);
        else
            fprintf('Warning: No anatomical label file for %s\n', sid);
        end
        % Prepare output folder
        out_sub = fullfile(output_dir, sid);
        if ~exist(out_sub, 'dir')
            mkdir(out_sub);
        end
        % Process
        try
            process_files(sid, fpath, lbl_path, out_sub);
            processed{end+1} = sid;
            fprintf('Processed %s\n', sid);
        catch ME
            warning('Failed %s: %s', sid, ME.message);
            failed{end+1} = sid;
        end
    end
    
    % 5) Create summary
    create_summary_file(output_dir, processed);
    fprintf('\nDone. Success: %d, Failures: %d\n', numel(processed), numel(failed));
end

%% Helper: inspect table format
function inspect_table_format(path)
    try
        opts = detectImportOptions(path);
        fprintf('  Columns: %s\n', strjoin(opts.VariableNames, ', '));
        T = readtable(path, opts);
        disp(T(1:min(3, height(T)), :));
    catch ME
        fprintf('  Error: %s\n', ME.message);
    end
end

%% Helper: process_files
function process_files(subj, loc_path, lbl_path, out_dir)
    % Read localization CSV
    opts = detectImportOptions(loc_path);
    opts.VariableNamingRule = 'preserve';
    Tloc = readtable(loc_path, opts);
    % Standardize coordinates
    Tcoord = standardize_table(Tloc, subj);
    writetable(Tcoord, fullfile(out_dir, 'electrodes_coords.tsv'), 'Delimiter', '\t');
    
    % Read anatomical labels if available
    if ~isempty(lbl_path) && isfile(lbl_path)
        opts2 = detectImportOptions(lbl_path);
        opts2.VariableNamingRule = 'preserve';
        Tlbl = readtable(lbl_path, opts2);
        writetable(Tlbl, fullfile(out_dir, 'electrodes_labels.csv'));
    else
        Tlbl = table();
    end
    
    % Combine and save
    combined = combine_data(Tcoord, Tlbl);
    save(fullfile(out_dir, 'electrodes_combined.mat'), 'combined', '-v7.3');
end

%% Helper: standardize_table (localization)
function coords_table = standardize_table(tbl, subject_id)
    cols = tbl.Properties.VariableNames;
    xcol = find_col(cols,{'x','coord_x','mni_x'});
    ycol = find_col(cols,{'y','coord_y','mni_y'});
    zcol = find_col(cols,{'z','coord_z','mni_z'});
    if isempty(xcol)||isempty(ycol)||isempty(zcol)
        error('Missing coordinate columns for %s', subject_id);
    end
    coords_table = table();
    coords_table.x = tbl.(xcol);
    coords_table.y = tbl.(ycol);
    coords_table.z = tbl.(zcol);
    namecol = find_col(cols,{'name','label','channel'});
    if ~isempty(namecol)
        coords_table.name = tbl.(namecol);
    else
        coords_table.name = arrayfun(@(i)sprintf('%s_%d',subject_id,i),1:height(coords_table),'Uniform',false)';
    end
end

%% Helper: find column
function name = find_col(cols,cands)
    name = '';
    for c=cands
        idx = find(strcmpi(cols,c{1}),1);
        if isempty(idx)
            idx = find(contains(lower(cols),lower(c{1})),1);
        end
        if ~isempty(idx)
            name = cols{idx};
            return;
        end
    end
end

%% Helper: extract subject ID
function sid = extract_subject_id(fname)
    tok = regexp(fname,'EM\d+','match');
    sid = tok{1};
end

%% Helper: find label file
function f = find_label_file(list,sid)
    f = [];
    for k=1:numel(list)
        if contains(list(k).name,sid)
            f = list(k); return;
        end
    end
end

%% Helper: combine_data
function combined = combine_data(Tc,Tl)
    combined.coords = [Tc.x Tc.y Tc.z];
    combined.names = Tc.name;
    if istable(Tl) && height(Tl)==height(Tc)
        if ismember('mni_linear_x',Tl.Properties.VariableNames)
            combined.mni_coords = [Tl.mni_linear_x Tl.mni_linear_y Tl.mni_linear_z];
        end
        if ismember('label',Tl.Properties.VariableNames)
            combined.labels = Tl.label;
        end
        if ismember('ictal_activity',Tl.Properties.VariableNames)
            combined.ictal_activity = Tl.ictal_activity;
        end
        if ismember('interictal_activity',Tl.Properties.VariableNames)
            combined.interictal_activity = Tl.interictal_activity;
        end
    end
end

%% Helper: create summary
function create_summary_file(output_dir, processed)
    fid = fopen(fullfile(output_dir,'subjects_summary.txt'),'w');
    fprintf(fid,'Processed subjects (%d):\n',numel(processed));
    for i=1:numel(processed)
        fprintf(fid,'  %d. %s\n',i,processed{i});
    end
    fclose(fid);
end
