function numericData = extractNodalDataBlock(logFile, stepNumber, outputFile)
    % extractNodalDataBlock extracts the nodal coordinates block from a log file
    %
    % It searches for:
    %   1. A line "Step = <stepNumber>"
    %   2. Two lines below, the line "Data = nodal coordinates"
    %   3. Then collects all subsequent lines until the next empty line
    %
    % Inputs:
    %   logFile     - path to .log file
    %   stepNumber  - numeric step to search for (e.g. 9)
    %   outputFile  - (optional) path to save the numeric data (as text)
    %
    % Output:
    %   numericData - numeric matrix of parsed values
    %
    % Example:
    %   coords = extractNodalDataBlock('simulation.log', 9, 'step9_coords.txt');

    if nargin < 3
        outputFile = '';
    end

    % --- Read all lines ---
    fid = fopen(logFile, 'r');
    if fid == -1
        error('Could not open file %s for reading.', logFile);
    end
    lines = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    fclose(fid);
    lines = lines{1};
    nLines = numel(lines);

    % --- Find "Step = N" line ---
    stepString = sprintf('Step = %d', stepNumber);
    stepIdx = find(strcmp(strtrim(lines), stepString));

    if isempty(stepIdx)
        warning('No line with "%s" found in file.', stepString);
        numericData = [];
        return;
    end

    numericData = []; % initialize

    for i = 1:length(stepIdx)
        checkLine = stepIdx(i) + 2; % two lines below

        if checkLine <= nLines && strcmp(strtrim(lines{checkLine}), 'Data = nodal coordinates')
            % Start collecting after that line
            startIdx = checkLine + 1;
            stopIdx = startIdx;
            while stopIdx <= nLines && ~isempty(strtrim(lines{stopIdx}))
                stopIdx = stopIdx + 1;
            end

            % Extract the text lines
            dataLines = lines(startIdx:stopIdx-1);

            % Convert each line to numeric
            try
                rowData = cellfun(@(line) sscanf(line, '%f')', dataLines, 'UniformOutput', false);
                rowData = vertcat(rowData{:});
            catch
                warning('Could not fully parse numeric data for Step = %d.', stepNumber);
                rowData = [];
            end

            numericData = [numericData; rowData]; %#ok<AGROW>
        end
    end

    % --- Optionally write to file ---
    if ~isempty(outputFile) && ~isempty(numericData)
        fid = fopen(outputFile, 'w');
        if fid == -1
            error('Could not open %s for writing.', outputFile);
        end
        fmt = [repmat('%g ', 1, size(numericData, 2)) '\n'];
        fprintf(fid, fmt, numericData');
        fclose(fid);
        fprintf('✅ Data block written to "%s"\n', outputFile);
    end

    if isempty(numericData)
        fprintf('⚠️  No nodal coordinate block found for Step = %d.\n', stepNumber);
    else
        fprintf('✅ Extracted nodal data for Step = %d (%d lines, %d columns)\n', ...
                stepNumber, size(numericData,1), size(numericData,2));
    end
end

%function dataBlock = extractNodalDataBlock(logFile, stepNumber, outputFile)
%    % extractNodalDataBlock extracts the nodal coordinates data block from a log file.
%    %
%    % It looks for a line "Step = <stepNumber>", checks two lines below for
%    % "Data = nodal coordinates", and if found, copies all following lines
%    % until an empty line.
%    %
%    % Inputs:
%    %   logFile     - path to the .log file
%    %   stepNumber  - integer step number to search for
%    %   outputFile  - optional output file path to save extracted block
%    %
%    % Output:
%    %   dataBlock   - cell array of strings (the extracted lines)
%    %
%    % Example:
%    %   dataBlock = extractNodalDataBlock('run.log', 9, 'step9_data.txt');
%
%    if nargin < 3
%        outputFile = '';
%    end
%
%    % Read all lines
%    fid = fopen(logFile, 'r');
%    if fid == -1
%        error('Could not open file %s for reading.', logFile);
%    end
%    lines = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
%    fclose(fid);
%    lines = lines{1};
%    nLines = numel(lines);
%
%    % Construct the step string
%    stepString = sprintf('Step = %d', stepNumber);
%
%    % Find all matching step lines
%    stepIdx = find(strcmp(strtrim(lines), stepString));
%
%    if isempty(stepIdx)
%        warning('No line with "%s" found in file.', stepString);
%        dataBlock = {};
%        return;
%    end
%
%    dataBlock = {};  % initialize
%
%    for i = 1:length(stepIdx)
%        checkLine = stepIdx(i) + 2;  % two lines below
%
%        if checkLine <= nLines && strcmp(strtrim(lines{checkLine}), 'Data = nodal coordinates')
%            % Start collecting from next line
%            startIdx = checkLine + 1;
%
%            % Find next empty line
%            stopIdx = startIdx;
%            while stopIdx <= nLines && ~isempty(strtrim(lines{stopIdx}))
%                stopIdx = stopIdx + 1;
%            end
%
%            % Extract block
%            block = lines(startIdx:stopIdx-1);
%            dataBlock = [dataBlock; block]; %#ok<AGROW>
%        end
%    end
%
%    % Write to file if requested
%    if ~isempty(outputFile)
%        fid = fopen(outputFile, 'w');
%        if fid == -1
%            error('Could not open %s for writing.', outputFile);
%        end
%        fprintf(fid, '%s\n', dataBlock{:});
%        fclose(fid);
%        fprintf('Data block written to "%s"\n', outputFile);
%    end
%
%    % Print info
%    fprintf('Found %d matching data block(s) for Step = %d.\n', numel(dataBlock) > 0, stepNumber);
%end
%
