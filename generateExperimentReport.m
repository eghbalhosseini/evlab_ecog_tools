function generateExperimentReport(obj, reportName)
     % Create the crunched folder if it doesn't exist
    crunchedFolder = fullfile(obj.crunched_file_path);
    if ~exist(crunchedFolder, 'dir')
        mkdir(crunchedFolder);
    end

    % Create a new PDF file in the crunched folder
    pdfFileName = fullfile(crunchedFolder, [reportName '.pdf']);
    if exist(pdfFileName, 'file')
        delete(pdfFileName);
    end

    % Start the PDF document with proper title page
    import mlreportgen.report.*
    import mlreportgen.dom.*
    
    % Create report object with title page
    rpt = Report(pdfFileName, 'pdf');
    
    % Add title page
    tp = TitlePage();
    tp.Title = 'Experiment Results Report';
    tp.Subtitle = 'Neural Data Analysis';
    tp.Author = 'Kumar Duraivel';
    tp.PubDate = datestr(now, 'dd-mmm-yyyy');
    add(rpt, tp);

    % General Experiment Information
    add(rpt, Heading1('General Experiment Information'));
    
    infoTable = {
        'Subject Name', obj.subject;
        'Experiment Name', obj.experiment;
        'Total Trials Completed', num2str(size(obj.events_table,1));
        'Total Electrodes Implanted', num2str(length(obj.elec_ch_label));
        'Total SEEG Electrodes', num2str(sum(contains(obj.elec_ch_type, 'seeg')));
        'Unipolar Contacts', num2str(length(obj.elec_ch_valid));
        'Bipolar Contacts', num2str(length(obj.bip_ch_label));
         'Significant Unipolar Electrodes (High Gamma)', num2str(sum(cellfun(@(x) any(x.h_sig_05), obj.stats.time_series.pSigChan)));
         'Significant Bipolar Electrodes (High Gamma)', num2str(sum(cellfun(@(x) any(x.h_sig_05), obj.stats.time_series.pSigChan_bip)));
    };
    
    tbl = Table(infoTable);
    tbl.Style = {Border('solid'), ColSep('solid'), RowSep('solid')};
    tbl.TableEntriesStyle = {HAlign('left')};
    add(rpt, tbl);

   % Data Extraction
    acc = [obj.events_table.accuracy];
    rt = [obj.events_table.RT];
    cond = [obj.condition];

    % Performance Metrics
    add(rpt, Heading2('Performance Metrics'));

    % Overall Accuracy
    accuracy_percentage = mean(acc)*100;
    add(rpt, Paragraph(['Overall accuracy: ' num2str(accuracy_percentage, '%.1f') '%']));

    % Reaction Time Analysis
    rt_correct = rt(acc == 1);
    cond_correct = cond(acc == 1);

    switch obj.experiment

        case {'LangLocVisual','LangLoc'}
           
            mean_rt_all = mean(rt_correct, 'omitnan');
            mean_rt_sentence = mean(rt_correct(strcmp(cond_correct, 'sentence')), 'omitnan');
            mean_rt_nonword = mean(rt_correct(strcmp(cond_correct, 'nonword')), 'omitnan');
        
            % Results Table
            tableContent = {
                'Condition', 'Mean RT (ms)';
                'Overall', sprintf('%6.1f ± %.1f', mean_rt_all*1000, std(rt_correct)*1000);
                'Sentence (S)', sprintf('%6.1f ± %.1f', mean_rt_sentence*1000, std(rt_correct(strcmp(cond_correct, 'sentence')))*1000);
                'Nonword (N)', sprintf('%6.1f ± %.1f', mean_rt_nonword*1000, std(rt_correct(strcmp(cond_correct, 'nonword')))*1000)
            };
            tbl = Table(tableContent);
            tbl.Style = {Border('solid'), ColSep('solid'), RowSep('solid')};
            tbl.TableEntriesStyle = {HAlign('center')};
            add(rpt, tbl);
        
            % Reaction Time Distribution Plot
            add(rpt, Heading2('Reaction Time Distributions'));
            debugMode = false; % Set to true for debugging
            if debugMode
                f = figure('Visible', 'on', 'Position', [100 100 800 600]);
            else
                f = figure('Visible', 'off', 'Position', [100 100 800 600]);
            end
            subplot(1,2,1)
            histogram(rt_correct(strcmp(cond_correct, 'sentence'))*1000, 'BinWidth', 50)
            title('Sentence Condition RT Distribution')
            xlabel('Reaction Time (ms)')
            ylabel('Frequency')
        
            subplot(1,2,2)
            histogram(rt_correct(strcmp(cond_correct, 'nonword'))*1000, 'BinWidth', 50)
            title('Nonword Condition RT Distribution')
            xlabel('Reaction Time (ms)')
            ylabel('Frequency')
        
            % Add the figure to the PDF
            add(rpt, Figure(f));
            close(f);

            % High Gamma Plot Generation (All trials)
            conds.A = find(strcmp(obj.condition, 'sentence'));
            conds.B = find(strcmp(obj.condition, 'nonword'));
            add(rpt, Heading2('High Gamma Plots (All Trials)'));
            high_gamma_plot(obj, conds, rpt);

            

            
        case {'LangLocAudio','LangLocAudio-2'}

            mean_rt_all = mean(rt_correct, 'omitnan');
            mean_rt_sentence = mean(rt_correct(strcmp(cond_correct, 'sentence')), 'omitnan');
            mean_rt_nonword = mean(rt_correct(strcmp(cond_correct, 'nonword')), 'omitnan');
        
            % Results Table
            tableContent = {
                'Condition', 'Mean RT (ms)';
                'Overall', sprintf('%6.1f ± %.1f', mean_rt_all*1000, std(rt_correct)*1000);
                'Sentence (S)', sprintf('%6.1f ± %.1f', mean_rt_sentence*1000, std(rt_correct(strcmp(cond_correct, 'sentence')))*1000);
                'Nonword (N)', sprintf('%6.1f ± %.1f', mean_rt_nonword*1000, std(rt_correct(strcmp(cond_correct, 'nonword')))*1000)
            };
            tbl = Table(tableContent);
            tbl.Style = {Border('solid'), ColSep('solid'), RowSep('solid')};
            tbl.TableEntriesStyle = {HAlign('center')};
            add(rpt, tbl);
        
            % Reaction Time Distribution Plot
            add(rpt, Heading2('Reaction Time Distributions'));
            debugMode = false; % Set to true for debugging
            if debugMode
                f = figure('Visible', 'on', 'Position', [100 100 800 600]);
            else
                f = figure('Visible', 'off', 'Position', [100 100 800 600]);
            end
            subplot(1,2,1)
            histogram(rt_correct(strcmp(cond_correct, 'sentence'))*1000, 'BinWidth', 50)
            title('Sentence Condition RT Distribution')
            xlabel('Reaction Time (ms)')
            ylabel('Frequency')
        
            subplot(1,2,2)
            histogram(rt_correct(strcmp(cond_correct, 'nonword'))*1000, 'BinWidth', 50)
            title('Nonword Condition RT Distribution')
            xlabel('Reaction Time (ms)')
            ylabel('Frequency')
        
            % Add the figure to the PDF
            add(rpt, Figure(f));
            close(f);
        
            % Audio Duration Histogram
            add(rpt, Heading2('Audio Duration Distribution'));
            f = figure('Visible', 'off', 'Position', [100 100 800 600]);
            audioDur=obj.events_table.audio_ended_natus-obj.events_table.audio_onset_natus;
            sentence_durations = audioDur(strcmp(obj.condition, 'sentence'));
            nonword_durations = audioDur(strcmp(obj.condition, 'nonword'));
            size(sentence_durations)
            subplot(1,2,1)
            histogram(sentence_durations, 'BinWidth', 0.1)
            title('Sentence Audio Duration Distribution')
            xlabel('Duration (s)')
            ylabel('Frequency')
            
        
            subplot(1,2,2)
            histogram(nonword_durations, 'BinWidth', 0.1)
            title('Nonword Audio Duration Distribution')
            xlabel('Duration (s)')
            ylabel('Frequency')
            
        
            figReporter = Figure(f);
            figReporter.Scaling = 'none';
            figReporter.Snapshot.ScaleToFit = true;
            add(rpt, figReporter);
            close(f);
        
            % High Gamma Plot Generation (All trials)
            conds.A = find(strcmp(obj.condition, 'sentence'));
            conds.B = find(strcmp(obj.condition, 'nonword'));
            add(rpt, Heading2('High Gamma Plots (All Trials)'));
            high_gamma_plot(obj, conds, rpt);
        
        
           
    end
     % Add summary report of significant channels
            add(rpt, Heading1('Summary of Significant Channels'));
            
            % Unipolar channels
            add(rpt, Heading2('Unipolar Channels with Significant Time Clusters'));
            sigUnipolarChannels = find(cellfun(@(x) any(x.h_sig_05), obj.stats.time_series.pSigChan));
            if ~isempty(sigUnipolarChannels)
                unipolarList = cell(length(sigUnipolarChannels), 1);
                for i = 1:length(sigUnipolarChannels)
                    unipolarList{i} = obj.elec_ch_label{sigUnipolarChannels(i)};
                end
                add(rpt, UnorderedList(unipolarList));
            else
                add(rpt, Paragraph('No significant unipolar channels found.'));
            end
            
            % Bipolar channels
            add(rpt, Heading2('Bipolar Channels with Significant Time Clusters'));
            sigBipolarChannels = find(cellfun(@(x) any(x.h_sig_05), obj.stats.time_series.pSigChan_bip));
            if ~isempty(sigBipolarChannels)
                bipolarList = cell(length(sigBipolarChannels), 1);
                for i = 1:length(sigBipolarChannels)
                    bipolarList{i} = obj.bip_ch_label{sigBipolarChannels(i)};
                end
                add(rpt, UnorderedList(bipolarList));
            else
                add(rpt, Paragraph('No significant bipolar channels found.'));
            end

    % Save the updated object in the crunched folder
    saveUpdatedObject(obj);

    % Close the PDF document
    
    
        close(rpt);
   
end

function high_gamma_plot(obj, conds, rpt)
    import mlreportgen.report.*
    import mlreportgen.dom.*
    
    % Define LangLoc settings
    langLocSettings = {
        'All Trials', conds;
        'Accurate Trials', struct('A', conds.A(obj.events_table.accuracy(conds.A) == 1), ...
                                  'B', conds.B(obj.events_table.accuracy(conds.B) == 1));
    };
    
    % Add audio duration condition for LangLocAudio
    % if strcmp(obj.experiment, 'LangLocAudio')
    %     audioDur = obj.events_table.audio_ended_natus - obj.events_table.audio_onset_natus;
    %     langLocSettings(length(langLocSettings)+1,:) = {'Trials_less_than_5_seconds', ...
    %         struct('A', conds.A(audioDur(conds.A) <= 5), ...
    %                'B', conds.B(audioDur(conds.B) <= 5))};
    %      langLocSettings(length(langLocSettings)+1,:) = {'Accurate_Trials_less_than_5_seconds', ...
    %         struct('A', conds.A(audioDur(conds.A) <= 5 & obj.events_table.accuracy(conds.A) == 1), ...
    %                'B', conds.B(audioDur(conds.B) <= 5 & obj.events_table.accuracy(conds.B) == 1))};
    % end

    % Data epoching
    epochTimeRange = [-0.5 6];
    [epochData, epochData_bip] = obj.extract_trial_epochs('epoch_tw', epochTimeRange, 'selectChannels', obj.elec_ch_clean);
    
    for iTrial = 1:size(epochData, 2)
        audioDurSample = obj.trial_timing{iTrial,:}.end(13)-obj.trial_timing{iTrial,:}.start(1)-epochTimeRange(1)*obj.sample_freq;
        mask = false(size(epochData, 3), 1);
        mask(audioDurSample+1:end) = true;
        
        epochData(:, iTrial, mask) = NaN;
        epochData_bip(:, iTrial, mask) = NaN;
    end

    % Process each LangLoc setting
    for i = 1:size(langLocSettings, 1)
        settingName = langLocSettings{i, 1};
        settingConds = langLocSettings{i, 2};
        
        add(rpt, Heading2(['High Gamma Plots: ' settingName]));
        
        % Process unipolar data
        add(rpt, Heading3('Unipolar Data'));
        obj = process_and_plot_data(obj, epochData, settingConds, 'unipolar', obj.elec_ch_label(obj.elec_ch_valid), epochTimeRange, rpt, settingName);
        
        % Process bipolar data
        if ~isempty(epochData_bip)
            add(rpt, Heading3('Bipolar Data'));
            obj = process_and_plot_data(obj, epochData_bip, settingConds, 'bipolar', obj.bip_ch_label, epochTimeRange, rpt, settingName);
        end
    end
end

function obj = process_and_plot_data(obj, data, conds, data_type, chanLab, epochTimeRange, rpt, settingName)
    data2plot = data;

    import mlreportgen.report.*
    import mlreportgen.dom.*
    
    % Time series plotting setup
    aTrials = conds.A;
    bTrials = conds.B;
    
    duration = size(data2plot, 3);
    x = linspace(epochTimeRange(1), epochTimeRange(2), duration);
    aColor = [0.8500 0.3250 0.0980]; % Red
    bColor = [0 0.4470 0.7410]; % Blue

    numChan = 10; % Number of channels to plot per figure
    totChanBlock = ceil(size(data2plot,1) / numChan);

    % Set figure size to match standard PDF page size
    pdfPageWidth = 10;  % inches
    pdfPageHeight = 12;  % inches
    figureWidth = pdfPageWidth * 100;  % Convert to pixels
    figureHeight = pdfPageHeight * 100;  % Convert to pixels
    
    % Initialize pSig cell array
    pSig = cell(1, size(data2plot,1));
    
    parfor iChan = 1:size(data2plot,1)
        aTrialData = squeeze(data2plot(iChan, aTrials, :));
        bTrialData = squeeze(data2plot(iChan, bTrials, :));

        pSig{iChan} = timePermCluster(aTrialData,bTrialData,numTail=1);
    end

    % Save pSig results in obj.stats
    if strcmp(data_type, 'unipolar')
        obj.stats.time_series.(['pSigChan_contrast_' strrep(settingName, ' ', '_')]) = pSig;
    else
        obj.stats.time_series.(['pSigChan_bip_contrast_' strrep(settingName, ' ', '_')]) = pSig;
    end

    for iF = 0:totChanBlock-1
        % Create a figure with the size of a standard PDF page
        f = figure('Visible', 'off', 'Position', [100, 100, figureWidth, figureHeight], 'Renderer', 'painters');
        
        % Define the grid layout based on the number of channels per figure
        numRows = 5; % 5 rows for 10 subplots (2 columns x 5 rows)
        numCols = 2; % 2 columns for 10 subplots

        for iChan = 1:min(numChan, size(data2plot,1) - iF*numChan)
            iChan2 = iChan + iF*numChan;
            
            % Calculate row and column for the current subplot
            row = floor((iChan - 1) / numCols) + 1;
            col = mod(iChan - 1, numCols) + 1;
            
            % Create subplot
            subplot(numRows, numCols, iChan);
            hold on;

            % Process both trial types
            for i = 1:2
                if i == 1
                    trialData = squeeze(data2plot(iChan2, aTrials, :));
                    trialColor = aColor;
                else
                    trialData = squeeze(data2plot(iChan2, bTrials, :));
                    trialColor = bColor;
                end
                trialMean = nanmean(trialData, 1);
                trialSEM = nanstd(trialData, 0, 1) / sqrt(size(trialData, 1));
                
                % patch([x fliplr(x)], [trialMean+trialSEM fliplr(trialMean-trialSEM)], ...
                %       trialColor, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                % plot(x, trialMean, 'Color', trialColor, 'LineWidth', 1);

                % Find valid (non-NaN) indices
                validIndices = ~isnan(trialMean) & ~isnan(trialSEM);
                
                % Extract valid data
                trialMean_cleaned = trialMean(validIndices);
                trialSEM_cleaned = trialSEM(validIndices);
                x_cleaned = x(validIndices);
                
                % Create the patch
                patch([x_cleaned fliplr(x_cleaned)], [trialMean_cleaned + trialSEM_cleaned fliplr(trialMean_cleaned - trialSEM_cleaned)], ...
                      trialColor, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                % plot the mean
                plot(x_cleaned, trialMean_cleaned, 'Color', trialColor, 'LineWidth', 1);

            end
            % Plot significance markers and horizontal line
            if strcmp(data_type, 'unipolar')
                sigTime = x(obj.stats.time_series.pSigChan{iChan2}.h_sig_05 == 1);
            else
                sigTime = x(obj.stats.time_series.pSigChan_bip{iChan2}.h_sig_05 == 1);
            end
            
            sigTimeCond = x(pSig{iChan2}.h_sig_05 == 1);
            scatter(sigTime, -0.5 * ones(size(sigTime)), 10, 'k', 'filled'); % significant HG activations
            scatter(sigTimeCond, -0.75 * ones(size(sigTimeCond)), 10, 'r', 'filled'); % significant langloc activations
            plot([epochTimeRange(1) epochTimeRange(2)], [0 0], 'k', 'LineWidth', 0.25);
            
            % Add X-label only for the last row
            if row == numRows
                xlabel('Time (s)');
            else
                xlabel('');
            end
            
            % Add Y-label only for the first column
            if col == 1
                ylabel('Z-score');
            else
                ylabel('');
            end
            
            % Add title with minimum font size
            title([chanLab{iChan2}]);
            
            % Adjust y-axis limits
            ylim([-1 3]);
            xlim(epochTimeRange);
            
            hold off;
        end
        
      % Create a new axes for the legend
        legendAxes = axes('Position', [0.1, 0.05, 0.8, 0.05], 'Visible', 'off');

        % Create dummy scatter plots for the legend
        hold(legendAxes, 'on');
        scatter(legendAxes, NaN, NaN, 10, 'k', 'filled');
        scatter(legendAxes, NaN, NaN, 10, 'r', 'filled');
        plot(legendAxes, NaN, NaN, 'Color', aColor, 'LineWidth', 2);
        plot(legendAxes, NaN, NaN, 'Color', bColor, 'LineWidth', 2);

        % Create the legend
        legend(legendAxes, {'Significant HG activations', 'Significant langloc activations', 'Sentence', 'Nonword'}, ...
               'Orientation', 'horizontal', 'Location', 'south');

        % Add footer with time-permutation cluster statistics details
        footerText = sprintf('Time-Permutation Cluster Test: nPerm=1000; pThresh=%.3f; numTail=1.', 0.05);
        annotation('textbox', [0.1, 0.01, 0.8, 0.03], 'String', footerText,...
                   'EdgeColor', 'none', 'HorizontalAlignment', 'center');

        figReporter0 = Figure(f);
        figReporter0.Scaling = 'none';
        figReporter0.Snapshot.ScaleToFit = true;
        add(rpt, figReporter0);
    
        close(f);
    end

end

function saveUpdatedObject(obj)
    % Create the crunched folder if it doesn't exist
    crunchedFolder = fullfile(obj.crunched_file_path);
    if ~exist(crunchedFolder, 'dir')
        mkdir(crunchedFolder);
    end

    % Generate the filename
    filename = fullfile(crunchedFolder, [obj.subject '_' obj.experiment '_crunched_HG_ZScore.mat']);

    % Save the object
    save(filename, 'obj', '-v7.3');
    
    fprintf('Updated object saved as: %s\n', filename);
end

