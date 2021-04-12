%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script prompts the user to select multiple .tif files movie for single
% molecule particle tracking analysis.  All of the results are saved within
% their respective folder where the selected movie is located. This one is 
% fully automated.  The user will only have to input parameters at the
% beginning.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CHECK GAUSSIAN FITTING ALGORITHM%%%%%%
clear all;
clc;
close all;

%% Select Movies

%addpath('./bioformats');

% load tif movie file
[files,filepath,numFiles] = SelectFilesToAnalyze;

%% Extract All Positions
tic
particleSiz = 4; %3 original, 5 good
brightThreshold =5000; %200 original, 5000 good
%don't include particles where the fwhm of the fitted gaussian is >6.5
%pixels
gaussWidth = 5;

showplot = 0;
showFit = 0; %0
showcnt = 0; %
set_parameters = [particleSiz brightThreshold gaussWidth];
for z = 1:numFiles
    tic
    close all
    disp(files{z})
    [mov,numFrames,dirpath,frameRate,fname] = LoadMovie(filepath{z},files{z});
    fits = []; %array containing number of fitted centroids for each frame
    % get flurophore positions
    pixelPos = [];
    intens_width =[];
   for frameIndex = 1:numFrames 
       disp(frameIndex)
        %for frameIndex = 1:1
        img = mov(:,:,frameIndex);
        
        % bandpass filter image
        imgb = bpass(double(img),1,particleSiz);

        % find peaks
        pk = pkfnd(double(imgb),brightThreshold,particleSiz);
        
        % if peaks are found, find centroids
        if ~isempty(pk)
            
            % calculate centroid moments
            cnt = cntrd(double(imgb),pk,particleSiz+1);
        
            
             % show movie with initial centroids
%             if showcnt == 1
%                 %img = mov(:,:,frameIndex); %
%                 figure(1); cla; hold on;
%                 imagesc(img); colormap gray;
%                 scatter(cnt(:,1),cnt(:,2),200,'gs');
%                 %scatter(cnt(badIndex,1),cnt(badIndex,2),200,'rs');
%                 axis tight
%                 title(['centroids ' num2str(frameIndex) ' out of ' num2str(numFrames)]);
%             end
            
            % 2D Gaussian fit around each centroid
            scale = 5; %5
            p0(3) = ceil(particleSiz/2);
            p0(4) = max(img(:));
            p0(5) = min(img(:));
            centroid = []; badIndex = []; error = []; sigma = []; fwhm_list=[];
            intensity = []; SNR = []; ratio = []; 
            
            fwhm_listG=[];
            
            for i = 1:size(cnt,1) %loop through every centroid
            %for i = 1:5 
                 try 
                    p0(1:2) = cnt(i,1:2); %get x y positions of centroid 
                    %[parameters,img,e1] = Gaussian2DFiterror(img,cnt(i,1:2),p0,scale,0);
                    %[parameters,img] = Gaussian2DFit(img,cnt(i,1:2),p0,scale,0);
                    
                    %perform gaussian fit around centroid
                    [parameters,imgbb] = Gaussian2DFit(img,cnt(i,1:2),p0,scale,0);
                    % standard deviation of fit
                    meanSigma = parameters(3); 
                    %FWHM
                    fwhm=meanSigma*2.355;
                    fwhm_list=[fwhm_list;fwhm];
                    %check that width of gauss fit is less than specified
                    %threshold value gaussWidth
                    if fwhm<gaussWidth
                        centroid = [centroid; parameters(1:2)]; % position of fit
                        N = parameters(4)*meanSigma*sqrt(2*pi); % number of photons
                        sigma = [sigma; meanSigma];
                        intensity = [intensity; N];
                        error = [error; meanSigma/sqrt(N)];
                        SNR = [SNR; parameters(4)/parameters(5)];
                        fwhm_listG=[fwhm_listG;fwhm];
                        
%                         fits1 = [frameIndex length(centroid)];
                    else
%                         fits1 = [frameIndex 0];
                    end
                   
                    
                 

                    % show 2D Gaussian Fit
                    if showFit == 1
                        disp(['SNR = ' num2str(SNR(i)) '; Error = ' num2str(error(i))]);
                        [sizey sizex] = size(img);
                        [X,Y] = meshgrid(1:sizey,1:sizex);
                        figure(1); cla; 
                        subplot(2,2,1); imagesc(img); colormap gray;
                        subplot(2,2,3); surf(img); 
                        subplot(2,2,2);
                        contour3(img,'--');
                        hold on;
                        contour3(Gaussian2D(parameters,X,Y),'-');
                    end
                    
               

                catch
                    badIndex = [badIndex; i];
                    disp('Error fitting Gaussian');
                end

            end
            fits1 = [frameIndex length(centroid)];
            fits = vertcat(fits,fits1); %add frame number v. # of centroids fitted to fits array
            % store fit info
            pixelPos = [pixelPos; centroid intensity sigma error SNR ones(size(centroid,1),1)*frameIndex];
            intens_width = [intens_width; intensity fwhm_listG ones(size(intensity,1),1)*frameIndex];
            
            
            % show movie with fitted centroids
            if showplot == 1
                %img = mov(:,:,frameIndex); %
                figure(2); cla; hold on;
                %figure; cla; hold on;
                imagesc(img); colormap gray;
                scatter(centroid(:,1),centroid(:,2),200,'gs');
                %scatter(cnt(badIndex,1),cnt(badIndex,2),200,'rs');
                axis tight
                title([num2str(frameIndex) ' out of ' num2str(numFrames)]);
            end
           
        end        
    end
    
    % save centroid positions
    dlmwrite(fullfile(dirpath,'pixelPositions.txt'),pixelPos,'delimiter','\t','newline','pc');
    dlmwrite(fullfile(dirpath,'parameters.txt'),set_parameters,'delimiter','\t','newline','pc');
    dlmwrite(fullfile(dirpath,'fits.txt'),fits,'delimiter','\t','newline','pc');
    
       % save intens_width as mat file 
    save(strcat(dirpath,'\intens_width.mat'),'intens_width')
    
    
    toc
end

    %% Extract one frame 5/15. Many edits December 2017 for gaussian fitting modifications
    
    particleSiz = 3; %3 original, 5 good
    brightThreshold = 5000; %200 original, 5000 good\
    %don't include particles where the fwhm of the fitted gaussian is >6.5
    %pixels
    gaussWidth = 5;
    frame_test=500;
    
    showplot = 1;
    showplotBad = 1;
    showFit = 0; %0
    showcnt = 0; %
    set_parameters = [particleSiz brightThreshold gaussWidth];
    for z = 1:numFiles
        tic
        close all
        disp(files{z})
        [mov,numFrames,dirpath,frameRate,fname] = LoadMovie(filepath{z},files{z});
        fits = []; %array containing number of fitted centroids for each frame
        % get flurophore positions
        pixelPos = [];
        for frameIndex = frame_test:frame_test
            img = mov(:,:,frameIndex);
            
            % bandpass filter image
            imgb = bpass(double(img),1,particleSiz);
            
            % find peaks
            pk = pkfnd(double(imgb),brightThreshold,particleSiz);
            
            % if peaks are found, find centroids
            if ~isempty(pk)
                
                % calculate centroid moments
                cnt = cntrd(double(imgb),pk,particleSiz+1);
                
                
                % show movie with initial centroids
                %             if showcnt == 1
                %                 %img = mov(:,:,frameIndex); %
                %                 figure(1); cla; hold on;
                %                 imagesc(img); colormap gray;
                %                 scatter(cnt(:,1),cnt(:,2),200,'gs');
                %                 %scatter(cnt(badIndex,1),cnt(badIndex,2),200,'rs');
                %                 axis tight
                %                 title(['centroids ' num2str(frameIndex) ' out of ' num2str(numFrames)]);
                %             end
                
                % 2D Gaussian fit around each centroid
                scale = 5; %5
                p0(3) = ceil(particleSiz/2);
                p0(4) = max(img(:));
                p0(5) = min(img(:));
                
%                 %include gaussian linear offset guess term: assume no
%                 %offset first
%                 p0(6)=0;
%                 p0(7)=0;
                
                intens_width =[];
                centroid = []; badIndex = []; error = []; sigma = []; fwhm_list=[];fwhm_listG=[];
                intensity = []; SNR = []; ratio = [];
                b_lin=[];d_lin=[];b_linB=[];d_linB=[];
                % pixels around centroid
                fit_img={}; 
                count=1; %keep track of number of fits with proper widths

                %fit within window
                window=[];
                
                %for bigger fits
                 centroidB = []; errorB = []; sigmaB = []; fwhm_listB=[];
                intensityB = []; SNRB = []; ratioB = [];
                % pixels around centroid
                fit_imgB={}; 
                countB=1; %keep track of number of fits with proper widths
                %fit within window
                windowB=[];
                %list of fitted centroids within gw
                cnt_orig=[];
                %list of fitted centroids bigger than gw
                cnt_origB=[];
                
                for i = 1:size(cnt,1) %loop through every centroid
                   % for i = 1:5
                    try
                        p0(1:2) = cnt(i,1:2); %get x y positions of centroid
                        %[parameters,img,e1] = Gaussian2DFiterror(img,cnt(i,1:2),p0,scale,0);
                        %[parameters,img] = Gaussian2DFit(img,cnt(i,1:2),p0,scale,0);
                        
                        %preform gaussian fit around centroid
                        %[parameters,imgbb] = Gaussian2DFit(img,cnt(i,1:2),p0,scale,0);
                        % try gaussian fit with linear term
                        [parameters,imgbb] = Gaussian2DFit(img,cnt(i,1:2),p0,scale,0);
                        % standard deviation of fit
                        meanSigma = parameters(3);
                        %FWHM
                        fwhm=meanSigma*2.355;
                        %check that width of gauss fit is less than specified
                        %threshold value gaussWidth
                        fwhm_list=[fwhm_list;fwhm];
                        
                        
                        if fwhm<gaussWidth
                            fit_img{count}=imgbb;
                            centroid = [centroid; parameters(1:2)]; % position of fit
                            N = parameters(4)*meanSigma*sqrt(2*pi); % number of photons
                            sigma = [sigma; meanSigma];
                            intensity = [intensity; N];
                            error = [error; meanSigma/sqrt(N)];
                            SNR = [SNR; parameters(4)/parameters(5)];
                            window=[window;parameters(6:7)];
                            %list of centroids withing gw
                            cnt_orig = [cnt_orig;cnt(i,1:2)];
%                             fits1 = [frameIndex length(centroid)];
                            fwhm_listG=[fwhm_listG;fwhm];
                            b_lin=[b_lin;parameters(6)];
                            d_lin=[d_lin;parameters(7)];    
                            count=count+1;
                            %look at peaks that were wider than GW
                        else
                            fit_imgB{countB}=imgbb;
                            centroidB = [centroidB; parameters(1:2)]; % position of fit
                            N = parameters(4)*meanSigma*sqrt(2*pi); % number of photons
                            sigmaB = [sigmaB; meanSigma];
                            intensityB = [intensityB; N];
                            errorB = [errorB; meanSigma/sqrt(N)];
                            SNRB = [SNRB; parameters(4)/parameters(5)];
                            windowB=[windowB;parameters(6:7)];
                            cnt_origB = [cnt_origB;cnt(i,1:2)];
                            countB=countB+1;
                            fwhm_listB=[fwhm_listB;fwhm];
                             b_linB=[b_linB;parameters(6)];
                            d_linB=[d_linB;parameters(7)];  
%                             fits1 = [frameIndex 0];
                        end
                        
                        
                        
                        %                     N = parameters(4)*meanSigma*sqrt(2*pi); % number of photons
                        %                     sigma = [sigma; meanSigma];
                        %                     intensity = [intensity; N];
                        %                     error = [error; meanSigma/sqrt(N)];
                        %                     SNR = [SNR; parameters(4)/parameters(5)];
                        %
                        %                     fits1 = [frameIndex length(centroid)];
                        
                        % show 2D Gaussian Fit
                        if showFit == 1
                            disp(['SNR = ' num2str(SNR(i)) '; Error = ' num2str(error(i))]);
                            [sizey sizex] = size(img);
                            [X,Y] = meshgrid(1:sizey,1:sizex);
                            figure(1); cla;
                            subplot(2,2,1); imagesc(img); colormap gray;
                            subplot(2,2,3); surf(img);
                            subplot(2,2,2);
                            contour3(img,'--');
                            hold on;
                            contour3(Gaussian2D(parameters,X,Y),'-');
                        end
                        
                        
                        
                    catch
                        badIndex = [badIndex; i];
                        disp('Error fitting Gaussian');
                    end
                    
                end
                fits1 = [frameIndex length(centroid)];
                fits = vertcat(fits,fits1); %add frame number v. # of centroids fitted to fits array
                % store fit info
                pixelPos = [pixelPos; centroid intensity sigma error SNR ones(size(centroid,1),1)*frameIndex];
                %store intensity and width of spot
                intens_width = [intens_width; intensity fwhm_listG];
                
                % show movie with fitted GOOD centroids
                if showplot == 1 && ~isempty(centroid)
                    %img = mov(:,:,frameIndex); %
                    figure(2); cla; hold on;
                    %figure; cla; hold on;
                    imagesc(img); colormap gray;
                    scatter(centroid(:,1),centroid(:,2),200,'gs');
                    %scatter(cnt(badIndex,1),cnt(badIndex,2),200,'rs');
                    axis tight
                    title([num2str(frameIndex) ' out of ' num2str(numFrames)]);
                end
                 % show movie with fitted BAD centroids if there are bad
                 % centroids
                if showplotBad == 1 && ~isempty(centroidB)
                    %img = mov(:,:,frameIndex); %
                    figure; cla; hold on;
                    %figure; cla; hold on;
                    imagesc(img); colormap gray;
                    scatter(centroidB(:,1),centroidB(:,2),200,'rs');
                    if ~isempty(centroid)
                   hold on;
                    scatter(centroid(:,1),centroid(:,2),200,'gs');
                    end
                    axis tight
                    title([num2str(frameIndex) ' out of ' num2str(numFrames)]);
                end
            end
        end
        
        % save centroid positions
        dlmwrite(fullfile(dirpath,'pixelPositions.txt'),pixelPos,'delimiter','\t','newline','pc');
        dlmwrite(fullfile(dirpath,'parameters.txt'),set_parameters,'delimiter','\t','newline','pc');
        dlmwrite(fullfile(dirpath,'fits.txt'),fits,'delimiter','\t','newline','pc');
        toc
    end
    
    % save intens_width as mat file 
    save(strcat(dirpath,'\intens_width.mat'),'intens_width')
    
     %% Extract Positions - one frame
%     
%     particleSiz = 4; %3 original, 5 good
%     brightThreshold = 200; %200 original, 5000 good
%     
%     showplot = 1;
%     showFit = 0; %0
%     showcnt = 0; %
%     set_parameters = [particleSiz brightThreshold];
%     for z = 1:numFiles
%         tic
%         close all
%         disp(files{z})
%         [mov,numFrames,dirpath,frameRate,fname] = LoadMovie(filepath{z},files{z});
%         fits = []; %array containing number of fitted centroids for each frame
%         % get flurophore positions
%         pixelPos = [];
%         frameIndex = 1000;
%         %for frameIndex = 1:numFrames
%         %for frameIndex = 1:1
%         if frameIndex == 1
%             img = mov(:,:,frameIndex);
%             
%             % bandpass filter image
%             imgb = bpass(double(img),1,particleSiz);
%             
%             % find peaks
%             pk = pkfnd(double(imgb),brightThreshold,particleSiz);
%             
%             % if peaks are found, find centroids
%             if ~isempty(pk)
%                 
%                 % calculate centroid moments
%                 cnt = cntrd(double(imgb),pk,particleSiz+1);
%                 
%                 
%                 % show movie with initial centroids
%                 %             if showcnt == 1
%                 %                 %img = mov(:,:,frameIndex); %
%                 %                 figure(1); cla; hold on;
%                 %                 imagesc(img); colormap gray;
%                 %                 scatter(cnt(:,1),cnt(:,2),200,'gs');
%                 %                 %scatter(cnt(badIndex,1),cnt(badIndex,2),200,'rs');
%                 %                 axis tight
%                 %                 title(['centroids ' num2str(frameIndex) ' out of ' num2str(numFrames)]);
%                 %             end
%                 
%                 % 2D Gaussian fit around each centroid
%                 scale = 5; %5
%                 p0(3) = ceil(particleSiz/2);
%                 p0(4) = max(img(:));
%                 p0(5) = min(img(:));
%                 centroid = []; badIndex = []; error = []; sigma = [];
%                 intensity = []; SNR = []; ratio = [];
%                 
%                 
%                 
%                 for i = 1:size(cnt,1) %loop through every centroid
%                     %for i = 1:5 %loop through every centroid
%                     try
%                         p0(1:2) = cnt(i,1:2); %get x y positions of centroid
%                         %[parameters,img,e1] = Gaussian2DFiterror(img,cnt(i,1:2),p0,scale,0);
%                         %[parameters,img] = Gaussian2DFit(img,cnt(i,1:2),p0,scale,0);
%                         [parameters,imgbb] = Gaussian2DFitlin(img,cnt(i,1:2),p0,scale,0);
%                         centroid = [centroid; parameters(1:2)]; % position of fit
%                         meanSigma = parameters(3); % standard deviation of fit
%                         N = parameters(4)*meanSigma*sqrt(2*pi); % number of photons
%                         sigma = [sigma; meanSigma];
%                         intensity = [intensity; N];
%                         error = [error; meanSigma/sqrt(N)];
%                         SNR = [SNR; parameters(4)/parameters(5)];
%                         
%                         fits1 = [frameIndex length(centroid)];
%                         
%                         % show 2D Gaussian Fit
%                         if showFit == 1
%                             disp(['SNR = ' num2str(SNR(i)) '; Error = ' num2str(error(i))]);
%                             [sizey sizex] = size(img);
%                             [X,Y] = meshgrid(1:sizey,1:sizex);
%                             figure(1); cla;
%                             subplot(2,2,1); imagesc(img); colormap gray;
%                             subplot(2,2,3); surf(img);
%                             subplot(2,2,2);
%                             contour3(img,'--');
%                             hold on;
%                             contour3(Gaussian2D(parameters,X,Y),'-');
%                         end
%                         
%                         
%                         
%                     catch
%                         badIndex = [badIndex; i];
%                         disp('Error fitting Gaussian');
%                     end
%                     
%                 end
%                 
%                 fits = vertcat(fits,fits1); %add frame number v. # of centroids fitted to fits array
%                 % store fit info
%                 pixelPos = [pixelPos; centroid intensity sigma error SNR ones(size(centroid,1),1)*frameIndex];
%                 
%                 % show movie with fitted centroids
%                 if showplot == 1
%                     %img = mov(:,:,frameIndex); %
%                     figure(2); cla; hold on;
%                     %figure; cla; hold on;
%                     imagesc(img); colormap gray;
%                     scatter(centroid(:,1),centroid(:,2),200,'gs');
%                     %scatter(cnt(badIndex,1),cnt(badIndex,2),200,'rs');
%                     axis tight
%                     title([num2str(frameIndex) ' out of ' num2str(numFrames)]);
%                 end
%                 
%                 % show specific frames
%                 %             if frameIndex == 1 || 1000
%                 %                 %img = mov(:,:,frameIndex); %
%                 %                 figure; cla; hold on;
%                 %                 imagesc(img); colormap gray;
%                 %                 scatter(centroid(:,1),centroid(:,2),200,'gs');
%                 %                 %scatter(cnt(badIndex,1),cnt(badIndex,2),200,'rs');
%                 %                 axis tight
%                 %                 title([num2str(frameIndex) ' out of ' num2str(numFrames)]);
%                 %             end
%             end
%         end
%         
%         % save centroid positions
%         dlmwrite(fullfile(dirpath,'pixelPositions.txt'),pixelPos,'delimiter','\t','newline','pc');
%         dlmwrite(fullfile(dirpath,'parameters.txt'),set_parameters,'delimiter','\t','newline','pc');
%         dlmwrite(fullfile(dirpath,'fits.txt'),fits,'delimiter','\t','newline','pc');
%         toc
%     end