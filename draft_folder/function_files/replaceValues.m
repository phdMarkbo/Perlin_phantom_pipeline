function M_out = replaceValues(M_in, pairs)
% replaceValues - replaces specified values in a matrix
%
% Usage:
%   M_out = replaceValues(M_in, pairs)
%
% Inputs:
%   M_in  - input matrix of any size
%   pairs - N×2 matrix where:
%           pairs(i,1) = value to replace
%           pairs(i,2) = replacement value
%
% Output:
%   M_out - matrix with replacements applied

    M_out = M_in; % start with input matrix unchanged

    for i = 1:size(pairs,1)
        oldVal = pairs(i,1);
        newVal = pairs(i,2);
        M_out(M_out == oldVal) = newVal;
    end
end

