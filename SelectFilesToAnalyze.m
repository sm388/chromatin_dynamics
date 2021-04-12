function  [filename,directory,numFiles] = SelectFilesToAnalyze(tmpDirectory)

if nargin < 1
    tmpDirectory = fileparts(pwd);
end

filename = {};
directory = {};
numFiles = 0;
status = 1;
while status == 1
    clc;
    for i = 1:numFiles
        disp([int2str(i) '. ' filename{i}]);
    end
    disp('__________________________________');
    disp('1. add another file');
    disp('2. remove a file');
    disp('3. Done');
    reply = input('Please select from above:');
    if isnumeric(reply) & reply > 0 & reply < 4
        switch(reply)
            case 1        
                [files,dirpath] = uigetfile(fullfile(tmpDirectory,'*.tif'),'Select mat files to analyze','Multiselect','on');
                if ~iscell(files)
                    numFiles = numFiles + 1;
                    filename{numFiles} = files;
                    directory{numFiles} = dirpath;
                else
                    for i = 1:length(files)
                        numFiles = numFiles + 1;
                        filename{numFiles} = files{i};
                        directory{numFiles} = dirpath;
                    end
                end
                tmpDirectory = dirpath;
            case 2
                reply2 = input('Select file to delete: ');
                if isnumeric(reply2) & reply2 > 0 & reply2 < numFiles+1
                    if numFiles > 0
                        numFiles = 0;
                        tmpFiles = {};
                        tmpdir = {};
                        for i = 1:length(filename)
                            if reply2 ~= i
                                numFiles = numFiles + 1;
                                tmpFiles{numFiles} = filename{i};
                                tmpdir{numFiles} = directory{i};
                            end
                        end
                        filename = tmpFiles;
                        directory = tmpdir;
                    else
                        clc;
                        disp('No files to delete... idiot');
                    end
                else
                    clc;
                    disp('Invalid response... idiot');
                end
            case 3
                status = 0;
            otherwise
                clc;
                disp('Invalid response... idiot');
        end
    end
end