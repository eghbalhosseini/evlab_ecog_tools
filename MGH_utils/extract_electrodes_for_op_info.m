function [electrode_labels,electrode_type] = extract_electrodes_for_op_info(cellArray)
%EXTRACT_ELECTRODES_FOR_OP_INFO Summary of this function goes here
%   Detailed explanation goes here
    stringParts = {};
    numberParts = {};
    electrode_labels={};
    electrode_type = {};
    % Loop through each element in the cell array
    for i = 1:length(cellArray)
        % Check if the current cell contains a string
        if ischar(cellArray{i}) || isstring(cellArray{i})
            % Use regular expression to find matches
            % The pattern is 'L' or 'R' followed by any letters and then some digits
            [matches, tokens] = regexp(cellArray{i}, '([LR][A-Za-z]*)(\d+)', 'match', 'tokens');

            % If a match is found, extract the string and number parts
            if ~isempty(matches)
                stringParts{end+1} = tokens{1}{1}; % String part
                numberParts{end+1} = str2double(tokens{1}{2}); % Number part, converted to a double
            end
        end
    end

    for j = 1:length(stringParts)
        electrode_labels{j, 1} = [stringParts{j}, '_', num2str(numberParts{j})];
    end

    % Create electrode_type array with 'stereo_eeg' for each entry
    electrode_type = repmat({'stereo_eeg'}, size(electrode_labels));

end

