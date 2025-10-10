%% Group-Level ECoG Analysis Using Existing Crunched Data
% Updated to handle macOS metadata files and proper TSV file loading

%% Helper function to filter valid MAT files
function valid_files = filter_valid_mat_files(file_list)
    valid_files = [];
    for i = 1:length(file_list)
        filename = file_list(i).name;
        % Exclude files that start with ._ (macOS metadata)
        % Exclude files that start with . (hidden files)
        % Only include .mat files
        if ~startsWith(filename, '.') && ~startsWith(filename, '._') && endsWith(filename, '.mat')
            % Check if file size is reasonable (> 1KB to avoid empty files)
            if file_list(i).bytes > 1024
                valid_files = [valid_files; file_list(i)];
            end
        end
    end
end

%% 1. Setup paths and parameters
base_path = '/Volumes/disk/nese/MGH_ECoG_Langloc/'; % Your actual path
crunched_path = [base_path 'crunched/'];
electrode_path = [base_path 'brainstom_localizations_cat12_fedorenko_nov_2024/'];
anatomical_path = [base_path 'language_atlas/'];

% Subject list from your folder structure
subjects = {'EM1036', 'EM1041', 'EM1042', 'EM1049', 'EM1050', 'EM1051', ...
           'EM1054', 'EM1057', 'EM1059', 'EM1065', 'EM1070', 'EM1071', ...
           'EM1078', 'EM1080', 'EM1083', 'EM1086', 'EM1091', 'EM1094', ...
           'EM1096', 'EM1100', 'EM1101', 'EM1102', 'EM1103', 'EM1118', ...
           'EM1122', 'EM1123', 'EM1125', 'EM1126', 'EM1127', 'EM1168', ...
           'EM1169', 'EM1175', 'EM1177'};

%% 2. Load all subjects' data and electrode coordinates
all_subjects_data = {};
all_electrode_coords = table();
all_electrode_labels = {};
all_anatomical_labels = {};
successful_loads = 0;
failed_loads = {};

fprintf('Loading data for %d subjects...\n', length(subjects));

for i = 1:length(subjects)
    subj = subjects{i};
    fprintf('Processing subject %s (%d/%d)\n', subj, i, length(subjects));
    
    try
        % Get all MAT files for this subject and filter them
        all_files = dir([crunched_path '*' subj '*.mat']);
        crunched_files = filter_valid_mat_files(all_files);
        
        if ~isempty(crunched_files)
            % Try to load the first valid file
            file_loaded = false;
            for f = 1:length(crunched_files)
                try
                    fprintf('  Attempting to load: %s\n', crunched_files(f).name);
                    data = load([crunched_path crunched_files(f).name]);
                    all_subjects_data{i} = data;
                    successful_loads = successful_loads + 1;
                    file_loaded = true;
                    fprintf('  Successfully loaded: %s\n', crunched_files(f).name);
                    break; % Exit loop once we successfully load a file
                catch load_error
                    fprintf('  Warning: Could not load %s: %s\n', ...
                            crunched_files(f).name, load_error.message);
                    continue; % Try next file
                end
            end
            
            if ~file_loaded
                fprintf('  No valid files found for subject %s\n', subj);
                failed_loads{end+1} = subj;
                continue;
            end
        else
            fprintf('  No MAT files found for subject %s\n', subj);
            failed_loads{end+1} = subj;
            continue;
        end
        
        % Load electrode coordinates - FIXED: Use dir() to find TSV files
        fprintf('  Looking for electrode files for %s...\n', subj);
        
        % Try multiple possible electrode file patterns
        possible_patterns = {
            [electrode_path 'sub-' subj '_ses-implantation-01_space-MNI-cat12_electrodes.tsv'],
            [electrode_path 'sub-' subj '_ses-implantation-02_space-MNI-cat12_electrodes.tsv'],
            [electrode_path 'sub-' subj '_*electrodes.tsv']
        };
        
        elec_file_found = false;
        
        for p = 1:length(possible_patterns)
            if contains(possible_patterns{p}, '*')
                % Use dir() for wildcard patterns
                elec_files = dir(possible_patterns{p});
                if ~isempty(elec_files)
                    % Filter out hidden files
                    valid_elec_files = elec_files(~startsWith({elec_files.name}, '.'));
                    if ~isempty(valid_elec_files)
                        elec_file = fullfile(valid_elec_files(1).folder, valid_elec_files(1).name);
                        fprintf('    Found electrode file: %s\n', valid_elec_files(1).name);
                        elec_file_found = true;
                        break;
                    end
                end
            else
                % Direct file check
                if exist(possible_patterns{p}, 'file')
                    elec_file = possible_patterns{p};
                    fprintf('    Found electrode file: %s\n', elec_file);
                    elec_file_found = true;
                    break;
                end
            end
        end
        
        if elec_file_found
            try
                % Load the TSV file
                fprintf('    Loading electrode coordinates from: %s\n', elec_file);
                
                % Try different methods to read the TSV file
                try
                    ch_RAS_tbl = readtable(elec_file, 'FileType', 'text', 'Delimiter', '\t');
                catch
                    % Alternative loading method
                    ch_RAS_tbl = readtable(elec_file, 'Delimiter', '\t');
                end
                
                % Check if we have the expected columns
                if any(ismember({'x', 'y', 'z'}, ch_RAS_tbl.Properties.VariableNames))
                    % Add subject identifier
                    ch_RAS_tbl.subject_id = repmat({subj}, height(ch_RAS_tbl), 1);
                    ch_RAS_tbl.subject_idx = repmat(i, height(ch_RAS_tbl), 1);
                    
                    % Rename coordinates to standard format
                    if ismember('x', ch_RAS_tbl.Properties.VariableNames)
                        ch_RAS_tbl.Properties.VariableNames{'x'} = 'R';
                    end
                    if ismember('y', ch_RAS_tbl.Properties.VariableNames)
                        ch_RAS_tbl.Properties.VariableNames{'y'} = 'A';
                    end
                    if ismember('z', ch_RAS_tbl.Properties.VariableNames)
                        ch_RAS_tbl.Properties.VariableNames{'z'} = 'S';
                    end
                    
                    % Concatenate to main table
                    if isempty(all_electrode_coords)
                        all_electrode_coords = ch_RAS_tbl;
                    else
                        % Ensure compatible table structures
                        common_vars = intersect(all_electrode_coords.Properties.VariableNames, ...
                                              ch_RAS_tbl.Properties.VariableNames);
                        all_electrode_coords = [all_electrode_coords(:, common_vars); 
                                              ch_RAS_tbl(:, common_vars)];
                    end
                    
                    fprintf('    Loaded %d electrodes\n', height(ch_RAS_tbl));
                else
                    fprintf('    Warning: TSV file does not contain expected columns (x,y,z)\n');
                    fprintf('    Available columns: %s\n', strjoin(ch_RAS_tbl.Properties.VariableNames, ', '));
                end
                
            catch coord_error
                fprintf('    Warning: Could not load electrode coordinates for %s: %s\n', ...
                        subj, coord_error.message);
            end
        else
            fprintf('    Electrode coordinate file not found for %s\n', subj);
            
            % Debug: List available files
            debug_files = dir([electrode_path 'sub-' subj '*']);
            if ~isempty(debug_files)
                fprintf('    Available files for %s:\n', subj);
                for df = 1:length(debug_files)
                    fprintf('      %s\n', debug_files(df).name);
                end
            else
                fprintf('    No files found matching pattern: sub-%s*\n', subj);
            end
        end
        
        % Load anatomical labels if available
        anat_file = [electrode_path 'sub-' subj '_channel_anatomical_labels.csv'];
        if exist(anat_file, 'file')
            try
                anat_labels = readtable(anat_file);
                all_anatomical_labels{i} = anat_labels;
                fprintf('    Loaded anatomical labels\n');
            catch anat_error
                fprintf('    Warning: Could not load anatomical labels for %s: %s\n', ...
                        subj, anat_error.message);
            end
        end
        
    catch subj_error
        fprintf('  Error processing subject %s: %s\n', subj, subj_error.message);
        failed_loads{end+1} = subj;
        continue;
    end
end

%% 3. Summary of data loading
fprintf('\n=== DATA LOADING SUMMARY ===\n');
fprintf('Successfully loaded: %d subjects\n', successful_loads);
fprintf('Failed to load: %d subjects\n', length(failed_loads));
if ~isempty(failed_loads)
    fprintf('Failed subjects: %s\n', strjoin(failed_loads, ', '));
end
fprintf('Total electrodes loaded: %d\n', height(all_electrode_coords));
if height(all_electrode_coords) > 0
    fprintf('Unique subjects with coordinates: %d\n', length(unique(all_electrode_coords.subject_id)));
    
    % Display electrode coordinate ranges for sanity check
    if ismember('R', all_electrode_coords.Properties.VariableNames)
        fprintf('Electrode coordinate ranges:\n');
        fprintf('  R: %.1f to %.1f mm\n', min(all_electrode_coords.R), max(all_electrode_coords.R));
        fprintf('  A: %.1f to %.1f mm\n', min(all_electrode_coords.A), max(all_electrode_coords.A));
        fprintf('  S: %.1f to %.1f mm\n', min(all_electrode_coords.S), max(all_electrode_coords.S));
    end
end

%% 4. Continue with analysis only if we have data
if successful_loads > 0
    fprintf('\nProceeding with analysis...\n');
    
    % Extract language-responsive electrodes across subjects
    fprintf('Extracting language-responsive electrodes...\n');
    
    language_electrodes = [];
    language_responses = [];
    all_sig_electrodes = table();
    
    for i = 1:length(all_subjects_data)
        if ~isempty(all_subjects_data{i})
            try
                data = all_subjects_data{i};
                
                % Debug: Display structure of loaded data
                fprintf('  Subject %s data structure:\n', subjects{i});
                fields = fieldnames(data);
                fprintf('    Top-level fields: %s\n', strjoin(fields, ', '));
                
                % Extract significant electrodes (adapt field names to your data structure)
                if isfield(data, 'obj') % If data is wrapped in obj structure
                    ecog_obj = data.obj;
                    fprintf('    Using data.obj\n');
                elseif isfield(data, 'ecog_obj')
                    ecog_obj = data.ecog_obj;
                    fprintf('    Using data.ecog_obj\n');
                else
                    % Find the main data structure
                    fields = fieldnames(data);
                    % Look for likely field names
                    likely_fields = fields(contains(fields, {'obj', 'data', 'ecog', 'results'}));
                    if ~isempty(likely_fields)
                        ecog_obj = data.(likely_fields{1});
                        fprintf('    Using data.%s\n', likely_fields{1});
                    else
                        ecog_obj = data.(fields{1});
                        fprintf('    Using data.%s (first field)\n', fields{1});
                    end
                end
                
                % Display structure of ecog_obj
                if isstruct(ecog_obj)
                    obj_fields = fieldnames(ecog_obj);
                    fprintf('    Object fields: %s\n', strjoin(obj_fields, ', '));
                end
                
                % Extract significance results (adapt field names as needed)
                sig_found = false;
                if isfield(ecog_obj, 's_vs_n_sig')
                    sig_results = ecog_obj.s_vs_n_sig;
                    fprintf('    Found s_vs_n_sig field\n');
                    
                    if isfield(sig_results, 'elec_data_dec') && ~isempty(sig_results.elec_data_dec)
                        sig_idx = find(sig_results.elec_data_dec{1});
                        fprintf('    Found %d significant electrodes\n', length(sig_idx));
                        
                        % Store significant electrodes with subject info
                        if ~isempty(sig_idx) && height(all_electrode_coords) > 0
                            subj_coords = all_electrode_coords(all_electrode_coords.subject_idx == i, :);
                            if height(subj_coords) >= max(sig_idx)
                                sig_coords = subj_coords(sig_idx, :);
                                
                                if isempty(all_sig_electrodes)
                                    all_sig_electrodes = sig_coords;
                                else
                                    all_sig_electrodes = [all_sig_electrodes; sig_coords];
                                end
                                
                                fprintf('    Subject %s: %d significant electrodes added\n', ...
                                        subjects{i}, length(sig_idx));
                                sig_found = true;
                            else
                                fprintf('    Warning: Electrode index mismatch for subject %s\n', subjects{i});
                            end
                        end
                    end
                elseif isfield(ecog_obj, 'stats') && isfield(ecog_obj.stats, 'sig_channels')
                    % Alternative field structure
                    sig_results = ecog_obj.stats.sig_channels;
                    fprintf('    Found stats.sig_channels field\n');
                    % Add processing for this structure...
                end
                
                if ~sig_found
                    fprintf('    No significance data found for subject %s\n', subjects{i});
                end
                
            catch extract_error
                fprintf('    Warning: Could not extract significance data for subject %d: %s\n', ...
                        i, extract_error.message);
            end
        end
    end
    
    fprintf('\nTotal significant electrodes found: %d\n', height(all_sig_electrodes));
    
    % Continue with visualization if we have significant electrodes
    if height(all_sig_electrodes) > 0 && height(all_electrode_coords) > 0
        fprintf('Generating visualizations...\n');
        
        % Simple visualization example
        figure('Position', [100 100 1200 400]);
        
        subplot(1,3,1);
        scatter3(all_electrode_coords.R, all_electrode_coords.A, all_electrode_coords.S, ...
                 20, [0.7 0.7 0.7], 'filled');
        hold on;
        scatter3(all_sig_electrodes.R, all_sig_electrodes.A, all_sig_electrodes.S, ...
                 60, 'r', 'filled');
        title('All Electrodes vs Language-Responsive');
        xlabel('R'); ylabel('A'); zlabel('S');
        axis equal; grid on;
        
        subplot(1,3,2);
        left_idx = all_sig_electrodes.R < 0;
        if any(left_idx)
            scatter3(all_sig_electrodes.R(left_idx), all_sig_electrodes.A(left_idx), ...
                     all_sig_electrodes.S(left_idx), 60, 'r', 'filled');
        end
        title('Left Hemisphere Language Electrodes');
        xlabel('R'); ylabel('A'); zlabel('S');
        axis equal; grid on; view([-90 0]);
        
        subplot(1,3,3);
        right_idx = all_sig_electrodes.R > 0;
        if any(right_idx)
            scatter3(all_sig_electrodes.R(right_idx), all_sig_electrodes.A(right_idx), ...
                     all_sig_electrodes.S(right_idx), 60, 'r', 'filled');
        end
        title('Right Hemisphere Language Electrodes');
        xlabel('R'); ylabel('A'); zlabel('S');
        axis equal; grid on; view([90 0]);
        
    elseif height(all_electrode_coords) > 0
        fprintf('Plotting all electrodes (no significance data available)...\n');
        
        figure('Position', [100 100 800 600]);
        scatter3(all_electrode_coords.R, all_electrode_coords.A, all_electrode_coords.S, ...
                 40, 'b', 'filled');
        title('All Loaded Electrodes');
        xlabel('R (mm)'); ylabel('A (mm)'); zlabel('S (mm)');
        axis equal; grid on;
        
    else
        fprintf('No electrode coordinates available for visualization\n');
    end
    
else
    fprintf('No data loaded successfully. Please check your file paths and data structure.\n');
end

%% 5. Save results if successful
if exist('all_sig_electrodes', 'var') && height(all_sig_electrodes) > 0
    save('group_analysis_results.mat', 'all_subjects_data', 'all_electrode_coords', ...
         'all_sig_electrodes', 'subjects', 'successful_loads', 'failed_loads');
    fprintf('\nResults saved to group_analysis_results.mat\n');
elseif height(all_electrode_coords) > 0
    save('group_electrode_coords.mat', 'all_electrode_coords', 'subjects', ...
         'successful_loads', 'failed_loads');
    fprintf('\nElectrode coordinates saved to group_electrode_coords.mat\n');
end
