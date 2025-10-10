function combinedArray = concatenateFieldValues(cellArray, fieldName)
    % concatenateFieldValues - Concatenate all elements from a specified field 
    %                          in a cell array of structs into one array
    %
    % Usage:
    %   combinedArray = concatenateFieldValues(cellArray, fieldName)
    %
    % Inputs:
    %   cellArray : Cell array where each cell contains a struct
    %   fieldName : Name of the field to extract and concatenate (string or char)
    %
    % Output:
    %   combinedArray : Array containing concatenated values from all structs

    % Validate inputs
    if ~iscell(cellArray)
        error('Input must be a cell array.');
    end

    if nargin < 2 || (~ischar(fieldName) && ~isstring(fieldName))
        error('You must specify a field name as a string or character vector.');
    end

    combinedArray = [];
    for i = 1:numel(cellArray)
        s = cellArray{i};
        if isstruct(s) && isfield(s, fieldName)
            vals = s.(fieldName);
            combinedArray = [combinedArray, vals(:)']; % Make sure vals is a row vector for concatenation
        else
            warning('Cell %d does not contain a struct with field "%s". Skipping.', i, fieldName);
        end
    end
end

