%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Description: ask the user to select a movie file (*.tif).  Then a 
%   folder is created with the movie name and the movie is moved inside
%   this folder.  Then it is loaded.
%
% Input:
%   tmpdir - the initial directory to look for the movie file
%
% Output:
%   mov - movie matrix from tif file
%   numFrames - number of frames in movie
%   frameRate - frame rate of movie
%   dirpath - directory the movie is in
%   fname - name of the movie
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mov,numFrames,dirpath,frameDuration,fname] = LoadMovie(filepath,filename)
warning off;

[tmp,fname] = fileparts(filename);

% create folder for movie if it does not exist
tmpdir = filepath(find(filepath(1:end-1)==filesep,1,'last')+1:end-1);
if ~strcmp(tmpdir,fname)
    dirpath = fullfile(filepath,fname);
    if exist(dirpath,'dir') ~= 7
        mkdir(dirpath);
    end
else
    dirpath = filepath;
end


% load the movie frames
status = 1;
i = 1;
while status == 1
    try
        mov(:,:,i) = imread(fullfile(filepath,filename),i);       
        i = i + 1;
    catch
        status = 0;
    end
end
numFrames = i - 1;


% get the frame rate
try
    movFile = mmreader(fullfile(filepath,filename));
    frameDuration = movFile.FrameRate;
catch
    disp('no frame rate found... using 0.032 s');
    frameDuration= 0.032;
end



