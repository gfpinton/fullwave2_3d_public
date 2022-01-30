function [filesize] = sizeOfFile(fname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: MAY 3, 2020
% Determine the size of a file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fid msg] = fopen(fname);
if(fid~=-1)
    fseek(fid, 0, 'eof');
    filesize = ftell(fid);
    fclose(fid);
else
    filesize=0;
end
