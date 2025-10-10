function appendValuesToFile_chunked(filename, values, chunkSize)
    % appendValuesToFile_chunked appends numeric values to a file efficiently in chunks.
    %
    % Inputs:
    %   filename - path to text file
    %   values - numeric vector
    %   chunkSize - how many values to write per batch (default 10000)

    if nargin < 3
        chunkSize = 10000;
    end

    N = numel(values);
    fid = fopen(filename, 'a');
    if fid == -1
        error('Could not open file %s for appending.', filename);
    end

    fprintf(fid, '\n--- New Data Block (%d values) ---\n', N);

    % Stream the data in chunks
    for i = 1:chunkSize:N
        idxEnd = min(i + chunkSize - 1, N);
        chunk = values(i:idxEnd);
        fprintf(fid, '%.6f\n', chunk);
    end

    fclose(fid);
end

