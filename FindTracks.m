function [pixelTracks,posTracks] = FindTracks(pixelPos,max_dsp,trackParam)

trackPositions = track(pixelPos(:,[1 2 6]),max_dsp,trackParam); %[1 2 7]
% Format: trackPositions = [x y frameIndex trackIndex]
    
% match positions to get intensity & error of positions
pixelTracks = [];
for i = 1:size(trackPositions,1)
    index = find(pixelPos(:,1)==trackPositions(i,1) & pixelPos(:,2) == trackPositions(i,2),1,'first');
    pixelTracks = [pixelTracks; pixelPos(index,1:6) trackPositions(i,3:4)];
end
% Format: posTracks = [x y intensity error frameIndex ratio SNR trackIndex]
% NO
  
numTracks = pixelTracks(end,end);
id = pixelTracks(:,end);

pixels2um = trackParam.pixels2um;
posTracks = pixelTracks;
posTracks(:,1:2) = posTracks(:,1:2)*pixels2um;

disp([num2str(numTracks) ' tracked particles']);

