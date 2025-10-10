function appendFormattedValuesToFile(filename, values, chunkSize)
    % appendFormattedValuesToFile appends multi-line formatted text blocks
    % for each numeric value in 'values' to the given file.
    %
    % Inputs:
    %   filename  - path to text file
    %   values    - numeric vector
    %   chunkSize - (optional) number of values per write batch (default: 5000)
    %
    % Each block written looks like:
    %   this value is {value}. that is quite a number. lorem ipsum
    %   just padding text here
    %   padding even more text but i also want the value printed on this line so here it is: {value}

    if nargin < 3
        chunkSize = 5000; % default: write in chunks of 5k blocks
    end

    N = numel(values);
    fid = fopen(filename, 'a');
    if fid == -1
        error('Could not open file %s for appending.', filename);
    end

	ii = 0;
    for i = 1:chunkSize:N
        idxEnd = min(i + chunkSize - 1, N);
        chunk = values(i:idxEnd);

        % Prebuild formatted text for the chunk
        block = '';
        for j = 1:numel(chunk)
		ii = ii +1;
            val = chunk(j);
            block = [block, sprintf([ ...
                '<material id="%i" name="element_tissue_%i" type="neo-Hookean">\n' ...
                '	<density>1</density><E>%.6f</E><v>0.49</v>\n' ...
                '</material>\n'], ...
                ii,ii, val)];
        end

        % Write the whole chunk to file
        fwrite(fid, block);
    end

    fclose(fid);
end

