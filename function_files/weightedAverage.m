function results = weightedAverage(val1, val2, ratios)
    % weightedAverage computes weighted averages of two values over a list of ratios
    %
    % Inputs:
    %   val1   - scalar, first value (weighted by ratio)
    %   val2   - scalar, second value (weighted by 1 - ratio)
    %   ratios - vector of ratios (each between 0 and 1)
    %
    % Output:
    %   results - vector of same size as ratios, containing weighted averages
    
    % Validate inputs
    if ~isscalar(val1) || ~isscalar(val2)
        error('val1 and val2 must be scalars.');
    end
    if ~isvector(ratios)
        error('ratios must be a vector.');
    end

    % Clamp ratios between 0 and 1 (optional)
    ratios = max(0, min(1, ratios));
    
    % Vectorized computation
    results = val1 .* ratios + val2 .* (1 - ratios);
end

