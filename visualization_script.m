%% Simulated Electrode Analysis Script
% This script demonstrates the IntracranialElectrodeVisualizer with simulated data

%% 1. Initialize the Visualizer
fprintf('Initializing IntracranialElectrodeVisualizer...\n');

viz = IntracranialElectrodeVisualizer( ...
    'freesurferHome', '/Applications/freesurfer/8.1.0', ...
    'subjectsDir',    '/Applications/freesurfer/8.1.0/subjects', ...
    'reconDir',       '/Volumes/disk/nese/MGH_ECoG_Langloc', ...
    'elecDataDir',    '/Users/dsuseendar/data/electrode_data' ...
);

%% 2. Load fsaverage Surfaces
fprintf('Loading fsaverage surfaces...\n');
viz.loadSubject('fsaverage', ...
    'surfaceTypes', {'pial','inflated','white'}, ...
    'hemispheres',  {'lh','rh'} ...
);

%% 3. Load Atlas
fprintf('Loading aparc atlas...\n');
viz.loadAtlas('fsaverage', 'aparc');

%% 4. Load Simulated Electrode Data
fprintf('Loading simulated electrode data...\n');

% Read the simulated electrode table
elecTable = readtable('simulated_electrodes.csv');

% Load into visualizer
viz.loadElectrodes({}, 'table', elecTable);

fprintf('Loaded %d electrodes from %d subjects\n', height(elecTable), ...
        length(unique(elecTable.subject)));

%% 5. Display electrode statistics
fprintf('\nElectrode Statistics:\n');
fprintf('Coordinate ranges - X: [%.1f, %.1f], Y: [%.1f, %.1f], Z: [%.1f, %.1f]\n', ...
        min(elecTable.x), max(elecTable.x), ...
        min(elecTable.y), max(elecTable.y), ...
        min(elecTable.z), max(elecTable.z));
fprintf('pLanA range: [%.3f, %.3f], Mean: %.3f\n', ...
        min(elecTable.pLanA), max(elecTable.pLanA), mean(elecTable.pLanA));
fprintf('High language electrodes (pLanA > 0.7): %d\n', sum(elecTable.pLanA > 0.7));

%% 6. Plot Pial Surface with All Electrodes
fprintf('\nPlotting pial surface with electrodes...\n');
viz.plotSurface('fsaverage', 'pial', 'both');
viz.plotElectrodes();  % Default uniform coloring

%% 7. Plot Electrodes Colored by Language Probability
fprintf('Plotting electrodes colored by pLanA...\n');
viz.plotSurface('fsaverage', 'inflated', 'both');
pLanA = elecTable.pLanA;
viz.plotElectrodes(pLanA, 0.2, 'both');

%% 8. Overlay Selected Anatomical Regions
fprintf('Overlaying anatomical regions...\n');
try
    % Note: Region names need to match FreeSurfer aparc labels exactly
    viz.plotAnnotation('both', 'aparc', {'superiortemporal', 'middletemporal'});
catch ME
    fprintf('Warning: Could not overlay regions - %s\n', ME.message);
    fprintf('This may be due to region name mismatch with FreeSurfer labels\n');
end

%% 9. Generate Time Series Data and Coordinates
fprintf('Generating simulated time series data...\n');

% Extract coordinates
coords = [elecTable.x, elecTable.y, elecTable.z];
n_electrodes = size(coords, 1);
n_timepoints = 200;

% Generate realistic neural time series
% Base oscillatory activity with noise
t = linspace(0, 2, n_timepoints);  % 2 seconds
base_freq = 10;  % 10 Hz base oscillation

timeSeries = zeros(n_electrodes, n_timepoints);
for i = 1:n_electrodes
    % Base oscillation with some randomness in frequency and phase
    freq_var = base_freq + randn() * 2;
    phase_var = rand() * 2 * pi;
    
    % Generate base signal
    signal = sin(2*pi*freq_var*t + phase_var);
    
    % Add language-related activity for high pLanA electrodes
    if elecTable.pLanA(i) > 0.7
        % Add language-specific gamma activity (30-80 Hz)
        gamma_freq = 30 + rand() * 50;
        gamma_activity = 0.3 * sin(2*pi*gamma_freq*t + rand()*2*pi);
        signal = signal + gamma_activity;
        
        % Add event-related response (simulated language task response)
        event_time = 0.5;  % Event at 0.5 seconds
        event_response = 0.5 * exp(-((t - event_time) / 0.2).^2);
        signal = signal + event_response;
    end
    
    % Add noise
    noise = 0.2 * randn(size(t));
    timeSeries(i, :) = signal + noise;
end

fprintf('Generated time series: %d electrodes x %d timepoints\n', ...
        size(timeSeries, 1), size(timeSeries, 2));

%% 10. Plot Electrode Density Across Multiple Subjects
fprintf('Computing and plotting electrode density...\n');
subjectList = {'EM1036','EM1041','EM1042','EM1049','EM1050'};

try
    viz.plotElectrodeDensity(subjectList, 8, 'fsaverage', ...
        'surfaceType', 'inflated', ...
        'threshold', 0.25, ...
        'showElectrodes', false, ...
        'colormap', flipud(hot) ...
    );
catch ME
    fprintf('Warning: Could not plot density - %s\n', ME.message);
    fprintf('This may require organized electrode data structure\n');
end

%% 11. Analysis Summary
fprintf('\n=== Analysis Summary ===\n');
fprintf('Total electrodes analyzed: %d\n', n_electrodes);
fprintf('Subjects: %s\n', strjoin(unique(elecTable.subject), ', '));
fprintf('High language probability electrodes (pLanA > 0.7): %d (%.1f%%)\n', ...
        sum(pLanA > 0.7), 100*sum(pLanA > 0.7)/length(pLanA));

% Language electrode locations
high_lang_idx = pLanA > 0.7;
high_lang_regions = elecTable.anatomical_label(high_lang_idx);
[unique_regions, ~, idx] = unique(high_lang_regions);
region_counts = accumarray(idx, 1);

fprintf('\nHigh language electrodes by region:\n');
for i = 1:length(unique_regions)
    fprintf('  %s: %d electrodes\n', unique_regions{i}, region_counts(i));
end

%% 12. Optional: Animate Time Series (commented out - requires implementation)
fprintf('Animating time series...\n');
try
    viz.animateTimeSeries(coords, timeSeries, ...
        'frameRate', 15, ...
        'filename', 'language_timeseries.mp4' ...
    );
catch ME
    fprintf('Note: Animation feature not yet implemented - %s\n', ME.message);
end

%% 13. Save Results
fprintf('\nSaving analysis results...\n');

% Save electrode data with time series
results = struct();
results.elecTable = elecTable;
results.coords = coords;
results.timeSeries = timeSeries;
results.high_language_electrodes = find(pLanA > 0.7);
results.analysis_timestamp = datetime('now');

save('electrode_analysis_results.mat', 'results');
fprintf('Results saved to electrode_analysis_results.mat\n');

% Optional: Save figure
try
    viz.saveVisualization('language_electrode_mapping.png', ...
        'format', 'png', ...
        'resolution', 300 ...
    );
    fprintf('Visualization saved to language_electrode_mapping.png\n');
catch ME
    fprintf('Note: Could not save visualization - %s\n', ME.message);
end

fprintf('\n=== Analysis Complete ===\n');

%%

% Set your paths
source_dir = '/Volumes/disk/nese/MGH_ECoG_Langloc';
output_dir = '/Users/dsuseendar/data/electrode_data';

% Run organization
organize_electrode_data(source_dir, output_dir);