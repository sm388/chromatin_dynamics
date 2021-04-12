clear all;
clc;
close all;

%% Select Movies

% load movie file
[files,filepath,numFiles] = SelectFilesToAnalyze;

%%
% over79len=[];
% over100len=[];
% over300len=[];
fullAnalysis = 1; %

for z = 1:numFiles
    close all;
    
    % load movie
    disp(files{z});
    if fullAnalysis == 1
        [mov numFrames dirpath frameDuration fname] = LoadMovie(filepath{z},files{z});
    else
        [tmp fname] = fileparts(files{z});
        dirpath = fullfile(filepath{z},fname);
        frameDuration = 0.058;
    end
    
    % load gaussian fit positions
    try
        pixelPos = dlmread(fullfile(dirpath,'pixelPositions.txt'));
        status = 1;
    catch
        disp('Pixel position does not exist');
        status = 0;
    end
    
    if status == 1
        pixelPos(isnan(pixelPos(:,end)),:) = [];
        if size(pixelPos,2) == 7
            pixelPos = pixelPos(:,[1:5 7]); %excluded column 6 (snr)
        end
        % [ x y intensity error frameIndex SNR]
        
        
        % Filter positions!!!!!!!!!!
        
        
        % Track particles
        minFrames = 50;
        max_dsp = 2;%10
        trackParam.mem = 0;  %changed from 10 to 2 on 11/8/17    % number of frames that the fluorophore can disappear %0
        trackParam.good = minFrames;     % minimum number of frames that is considered a good track
        trackParam.dim = 2;
        trackParam.quiet = 1;
%    trackParam.pixels2um = .1605; %pixel size for 100x  - 0.1605 um/pixel %coolSNAP: 0.0645 % for 51.36um = 0.10
 trackParam.pixels2um = .065; %mochrie lab movies
%        trackParam.pixels2um = .1; %10/18/19 x16
        set_parameters = [minFrames max_dsp trackParam.mem];
        try
            %[pixelTracks, posTracks] = FindTracks(filteredPixelPos,max_dsp,trackParam);
            [pixelTracks, posTracks] = FindTracks(pixelPos,max_dsp,trackParam); %
            trackStatus = 1;
        catch
            trackStatus = 0;
        end
        
        if trackStatus == 1
            % save tracks
            dlmwrite(fullfile(dirpath,'PixelTracks.txt'),pixelTracks,'delimiter','\t','newline','pc');
            dlmwrite(fullfile(dirpath,'PosTracks.txt'),posTracks,'delimiter','\t','newline','pc');
            dlmwrite(fullfile(dirpath,'analysis_parameters.txt'),set_parameters,'delimiter','\t','newline','pc');
        else
            disp(['Errors tracking in ' fname]);
        end
    end
% %     tracklength = [];
%     for i=1:max(posTracks(:,8))
%         tracklength = [tracklength;length(find(posTracks(:,8)==i))];
%     end
    
%     over79 = find(tracklength>=79);
%     over300 = find(tracklength>=300);
%     over100 = find(tracklength>=100);
%     
%     over79len = [over79len;length(over79)];
%     over300len = [over300len;length(over300)];
%     over100len = [over100len;length(over100)];
%     
    
end

%% use this to plot tracks from main_analysis. Uses pixelTracks variable from main_Analysis
frameIndex=5; 
img = mov(:,:,frameIndex);
figure; imagesc(img); colormap gray;hold on;
m=max(posTracks(:,8));

%use pixel positions, not um positions, when plotting 
for i=1:m
    if i < m
    group=find(posTracks(:,8)==i);
    x=pixelTracks(group(1):group(end),1);
    y=pixelTracks(group(1):group(end),2);
    plot(x,y);
    hold on;
    else
    group=find(posTracks(:,8)==i);
    x=pixelTracks(group(1):group(end),1);
    y=pixelTracks(group(1):group(end),2);
    plot(x,y);
    hold off;
    end
end
    
%% find lengths of tracks 

tracklength = [];
for i=1:max(posTracks(:,8))
tracklength = [tracklength;length(find(posTracks(:,8)==i))];
end

over79 = find(tracklength>=79);
over300 = find(tracklength>=300);
over100 = find(tracklength>=100);
%% only plot tracks that are over 79 steps 

frameIndex=5; 
img = mov(:,:,frameIndex);
figure; imagesc(img); colormap gray;hold on;
m=max(posTracks(:,8));

%use pixel positions, not um positions, when plotting 
for i=1:length(over79)
    if i < length(over79)
    group=find(posTracks(:,8)==over79(i));
    x=pixelTracks(group(1):group(end),1);
    y=pixelTracks(group(1):group(end),2);
    plot(x,y);
    hold on;
    else
    group=find(posTracks(:,8)==over79(i));
    x=pixelTracks(group(1):group(end),1);
    y=pixelTracks(group(1):group(end),2);
    plot(x,y);
    hold off;
    end
end
















