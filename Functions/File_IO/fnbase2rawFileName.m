function [ rawFilename, error ] = fnbase2rawFileName( fnbase, ch, ext )
    % Returns raw file name for both 2 and 3 digit channel number file name
    % formats
    % error = 1 if file not found
    rawFilenameStruct = dir([fnbase,'*',num2str(ch),ext]) ;
    if (size(rawFilenameStruct,1) == 0 )
        error = 1 ; % File doesn't exist
        rawFilename = '' ;
    else
        error = 0 ;
        rawFilename = rawFilenameStruct.name ;
    end
end

