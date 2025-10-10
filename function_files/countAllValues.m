function totalCount = countAllValues(cellArray, fieldName)
    % countAllValues - Count total number of elements in a specified field 
    %                  across a cell array of structs
    %
    % Usage:
    %   totalCount = countAllValues(cellArray, fieldName)
    %
    % Inputs:
    %   cellArray : Cell array where each cell contains a struct
    %   fieldName : Name of the field to count elements from (string or char)
    %
    % Output:
    %   totalCount : Total number of elements across all specified fields

    % Validate inputs
    if ~iscell(cellArray)
        error('Input must be a cell array.');
    end

    if nargin < 2 || ~ischar(fieldName) && ~isstring(fieldName)
        error('You must specify a field name as a string or character vector.');
    end

    totalCount = 0;
    for i = 1:numel(cellArray)
        s = cellArray{i};
        if isstruct(s) && isfield(s, fieldName)
            totalCount = totalCount + numel(s.(fieldName));
        else
            warning('Cell %d does not contain a struct with field "%s". Skipping.', i, fieldName);
        end
    end
end

