function appendFileToFeb(febFile, textFile)
    % appendFileToFeb appends the contents of a text file to a FEBio .feb file.
    %
    % Inputs:
    %   febFile  - path to .feb file (string or char)
    %   textFile - path to text file whose contents you want to append
    
    % Read the text file
    fidIn = fopen(textFile, 'r');
    if fidIn == -1
        error('Could not open text file %s for reading.', textFile);
    end
    textData = fread(fidIn, '*char')';  % read entire text as char
    fclose(fidIn);

    % Open the FEB file for appending
    fidOut = fopen(febFile, 'a');
    if fidOut == -1
        error('Could not open FEB file %s for appending.', febFile);
    end

    % Append with a separating newline
    fprintf(fidOut, '\n%s\n', textData);
    fclose(fidOut);

    fprintf('Appended contents of "%s" to "%s".\n', textFile, febFile);
end

