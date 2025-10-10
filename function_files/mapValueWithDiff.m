function mappedValues = mapValuesWithDiff(cellArray, diffValue)
    % mapValuesWithDiff maps values with step diffValue per unit increase,
    % starting from ratio 1 at value=1.
    %
    % Inputs:
    %   cellArray - cell array with structs having field 'value'
    %   diffValue - scalar, decrement per unit increase in value
    %
    % Output:
    %   mappedValues - mapped ratios between 0 and 1 (clamped)
    
    N = numel(cellArray);
    mappedValues = zeros(1, N);
    
    for i = 1:N
        v = cellArray{i}.avgValue;
        ratio = 1 - (v - 1) * diffValue;
        
        % Clamp between 0 and 1 (optional)
        mappedValues(i) = max(0, min(1, ratio));
    end
end

