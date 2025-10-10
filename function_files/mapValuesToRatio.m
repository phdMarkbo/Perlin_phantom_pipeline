function mappedValues = mapValuesToRatio(cellArray, maxValue, minRatio)
    % mapValuesToRatio maps struct values from [1, maxValue] to [1, minRatio]
    %
    % Inputs:
    %   cellArray - 1xN cell array with structs having field 'value'
    %   maxValue - scalar, max expected value
    %   minRatio - scalar in [0,1], ratio that maxValue maps to
    %
    % Output:
    %   mappedValues - 1xN array with mapped ratios
    
    N = numel(cellArray);
    mappedValues = zeros(1, N);
    
    for i = 1:N
        v = cellArray{i}.avgValue;
        
        % Clamp value between 1 and maxValue
        v = max(1, min(v, maxValue));
        
        % Linear mapping to [1, minRatio]
        mappedValues(i) = 1 - ((v - 1) / (maxValue - 1)) * (1 - minRatio);
    end
end

