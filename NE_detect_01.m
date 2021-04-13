
% Script to try for spots/objects detection for NE detection for 2D NE-labeled
% images. Also expects to have another fluo channel with labled nuclear
% protein (diffuse/spotty)
%

%% INTIALIZATION: specify general parameters, constants etc

dx=0.1605; % pixel size
%dz=0.2; % z-stack step

%=========================================================================

%% Provide file information

path='C:\TMP\2020_03_04_MKSP_3140';
filename='3140_no_MBC_02_R3D.dv';
file_type='dv'; % or 'dv'
channel_1=1; % what channel to use from the multichannel image
channel_2=2; % what channel to use from the multichannel image

% path and name to save figures and Matlab variables
save_fig=0;
path_save='C:\_Users_\Ivan\NE dynamics\Nucleus drift estimate';
file_save='3140_no_MBC_Test';

file_suffix=filename(end-8:end-7);


%=========================================================================
%%  code INITIALIZATION
% to take of difefrences between Mac and PC in path representation
if ismac
    bf_slash='/';
elseif ispc
    bf_slash='\';
else
    disp('Platform not supported')
end


%=========================================================================

%% LOAD image data

% load image stack
switch file_type
    case 'tif'
        [Fluo_1, image_info]=load_tiff_stack([path, bf_slash, filename]);
        % n_frames=image_info(3);
        % im_size=image_info(1:2);
    case 'dv'
        % loading using bfread: filename, load 1-st series, all timepoints (or
        % specify the range), n-th channel...
        Fluo_1 = bfread([path, bf_slash, filename], 1, 'TimePoints', 'all', 'Channel', channel_1);
        Fluo_2 = bfread([path,  bf_slash, filename], 1, 'TimePoints', 'all', 'Channel', channel_2);
        
        if iscell(Fluo_1)
            Fluo_1=dv_movie_2_matrix(Fluo_1);
            Fluo_2=dv_movie_2_matrix(Fluo_2);
        end
        Fluo_1=uint16(Fluo_1);
        Fluo_2=uint16(Fluo_2);
        %disp ('not working yet, under costruction')
end
% get image dimensions
im_size_3=size(Fluo_1);
n_frames=im_size_3(3);
im_size=im_size_3(1:2);


%=========================================================================
%% PREPARE SPECIAL IMAGE for nuclei detection

% combine Fluo_1 diffuse nuclear signal with Fluo_2 nucleus edge signal

f1=max(Fluo_1(:))-min(Fluo_1(:));
f2=max(Fluo_2(:))-min(Fluo_2(:));

FluoX=uint16(mean(Fluo_1,3))+uint16(round(f2/f1))*uint16(mean(Fluo_2,3));


%=========================================================================

%% NUCLEI DETECTION using irregular spot detection

disp('Starting spot detection...')

param.peakRadius=20; % expected typical radius (or a bit more) of the nuceleus
param.intensityRatioThreshold=1.15;
param.shellThickness=1;
param.edgeDist=30; % want to be beyongd the nucleus
param.centerDist=2;
param.fitRadius=10;

% create fake 3D image since we using 3D algorithmm, will be replaced with
% 2D one

tic
Fluo0=uint16(zeros(im_size(1),im_size(2),3));
Fluo0(:,:,2)=FluoX(:,:);

% find approximate spots(i.e. nuclei) centers
[spotData, spotDetection_params]= findIrregularSpots3_only(Fluo0,'peakRadius',param.peakRadius,...
    'intensityRatioThreshold',param.intensityRatioThreshold,...
    'shellThickness',param.shellThickness,...
    'edgeDist',param.edgeDist,...
    'centerDist',param.centerDist, ...
    'fitRadius',param.fitRadius);

toc

disp('DONE with the spot detection!!!')


%=========================================================================
%% SHOW RESULTS of NUCLEI (i.e. spot) detection

% setup the image scaling
im_min=min(FluoX(:));
im_max=max(FluoX(:));
im_range=im_max-im_min;
im_0=im_min+0.0*im_range;
im_1=im_min+0.5*im_range;

% show results as max_projection
image=FluoX;
figure;
imshow(image,[im_0,im_1],'InitialMagnification',200,'Border','tight');
hold on

for ii=1:length(spotData)
    pos=spotData{ii}.spotPosition;
    plot(pos(1),pos(2),'+g','LineWidth',2)
    text(pos(1)+5,pos(2)+5,num2str(ii),'Color','g')
    
end
%plot(potentialXYZ(:,1),potentialXYZ(:,2),'or')

if save_fig
    savefig(gcf,[path_save,file_save,file_suffix,' SpotsDetection'],'compact')
end


%=========================================================================
%% RING DETECTION 1st frame

Rfit=20; % sets radius for the image square to be fitted
show_fit=1;
sc=[0.1,0.8]; % low and upper scale for image intenisty relative to min and max values, only for image display
frame=1;
% spot=1;
Rnuc_guess=8;

phi=0:2*pi/100:2*pi;

Fluo_2_ROI=zeros(im_size+2*Rfit);
Fluo_2_ROI(Rfit+1:end-Rfit, Rfit+1:end-Rfit)=Fluo_2(:, :, frame);

fit_results=zeros(length(spotData),9);
coeff_w_int=zeros(length(spotData),18);

for spot=1:length(spotData)
    oneSpotData=spotData{spot};
    spot_XYZ=oneSpotData.spotPosition;
    spot_RC=round([spot_XYZ(2),spot_XYZ(1)]);
    
    %  get 1st frame image
    % assuming that Fluo_2 is one with NE marker...
    image=Fluo_2_ROI(spot_RC(1)-Rfit+Rfit:spot_RC(1)+Rfit+Rfit, spot_RC(2)-Rfit+Rfit:spot_RC(2)+Rfit+Rfit, frame);
    im_size_0=size(image);
    roi=[0,0; 0,im_size_0(2); im_size_0(1),im_size_0(2); im_size_0(1),0];
    
    
    
    % prepare x and y and z values for fitting
    [x_values,y_values,px_values,roi_indx] = extractCellPixels_ROI(image,roi);
    % define fit function:
    model = fittype('a*circle_w_PSF(x, y, sigma, xc, yc, R)+c',...
        'dependent',{'z'},'independent',{'x','y'},...
        'coefficients',{'a','sigma','xc','yc','R','c'});
    % intial estimates for coefficients:
    coeff_0= [max(px_values(:))-min(px_values(:)), 2, im_size_0(1)/2, im_size_0(1)/2,  Rnuc_guess, min(px_values(:))];
    % set bounadries for coefficeints:
    lb=[0, 1, 0.2*im_size_0(1), 0.2*im_size_0(2), 2, 0];
    ub=[Inf,im_size_0(1)/2, 0.8*im_size_0(1), 0.8*im_size_0(2),min(im_size)/2,max(px_values(:))];
    % let's try to fit
    [fit1,gof1,fit_out] = fit([x_values,y_values],px_values,model,'Startpoint',coeff_0,'Lower',lb,'Upper',ub);
    % get coeff values, conf intervals ans smoe goodness of fit measures
    coeff=coeffvalues(fit1);
    a=coeff(1);
    sigma=coeff(2);
    xc=coeff(3);
    yc=coeff(4);
    R=coeff(5);
    c=coeff(6);
    coeff_int=confint(fit1);
    delta_coeff_int=diff(coeff_int,1);
    av_err_sum=sum(delta_coeff_int./coeff)/6;
    adjRsq=gof1.adjrsquare;
    rmse=gof1.adjrsquare;
    
    % collest results:     1       2         3  4  5  6    7            8       9
    %                                   a   sigma    x  y  R  c AvErr  RMSE aR2
    fit_results(spot,:)=[a, sigma, xc+spot_RC(2)-(Rfit+1),  yc+spot_RC(1)-(Rfit+1), R, c, av_err_sum, rmse, adjRsq];
    coeff_w_int(frame,1:3:18)=coeff;
    coeff_w_int(frame,2:3:18)=coeff_int(1,:);
    coeff_w_int(frame,3:3:18)=coeff_int(2,:);
    
    if show_fit
        if spot==1
            figure; hold off
        end
        hold off
        imshow(image,[min(image(:)),sc(1)*min(image(:))+sc(2)*(max(image(:))-min(image(:)))],'InitialMagnification',1600,'Border','tight');
        set(gcf, 'Name', ['spot #', num2str(spot), ' R = ',num2str(R)]);
        hold on
        plot((xc-R*cos(phi)), yc-R*sin(phi),':c','LineWidth',2)
        plot(xc, yc,'+c','LineWidth',1)
        text(xc-R, yc-R-2, ['AvRelEr =', num2str(av_err_sum)],'Color','b')
        text(xc+R/2, yc-R-2, ['RMSE=', num2str(rmse)],'Color','b')
        text(xc-R, yc+R+2, ['AdjR2 =', num2str(adjRsq)],'Color','b')
        pause (0.1)
    end
end

n_nuc=size(fit_results,1);

disp('Done with NE ring fitting.')

%=========================================================================
%% SHOWING RING DETECTION RESULTS for one frame
% showing all resulst for the frame
figure
image=Fluo_2(:, :, frame);
imshow(image,[min(image(:)),sc(1)*min(image(:))+sc(2)*(max(image(:))-min(image(:)))],'InitialMagnification',1600,'Border','tight');
set(gcf, 'Name', ['spot #',num2str(spot), ' R = ',num2str(R)]);
hold on
for spot=1:length(spotData)
    xc_yc_R=fit_results(spot,3:5);
    plot((xc_yc_R(1)+xc_yc_R(3)*cos(phi)), xc_yc_R(2)+xc_yc_R(3)*sin(phi),'-g','LineWidth',1)
    plot(xc, yc,'+g','LineWidth',1)
    text(xc_yc_R(1)+xc_yc_R(3), xc_yc_R(2)+xc_yc_R(3), num2str(spot),'Color','g')
end

%=========================================================================
%% FILTER NUCLEI RINGs
% Nuclei will be filtered to keep only well-fitted nucleis, based on total relative error of 
% fitting parameters. Also remove "overlapping" nuclei (usually coming from nuclei undegoing mitosis)

use_filter=1;

max_RelEr=0.1;

if use_filter
    fit_results_0=fit_results;
    
    ind_0=true(n_nuc,1);
    % Total Rel Error filter
    ind_1=fit_results(:,7)<max_RelEr;
    % overlapping nuclei filter
    x=fit_results(:, 3);
    y=fit_results(:, 4);
    r=fit_results(:, 5);
    
    XX=repmat(x', n_nuc, 1);
    YY=repmat(y, 1, n_nuc);
    RRi=repmat(r', n_nuc, 1);
    RRj=repmat(r, 1, n_nuc);
    
    RR2=(XX-XX').^2+(YY-YY').^2;
    
    ind_2= (RR2>(RRi+RRj).^2);
    ind_2=ind_2+eye(n_nuc);
    [ind_2i, ind_2j]=find(ind_2==0);
    
    ind_2=ind_0;
    ind_2(ind_2i)=false;
    
    % combine all filters in one
    ind=ind_1 & ind_2;
    
    fit_results=[(1:n_nuc)', fit_results];
    
    fit_results_F=fit_results(ind,:);
    
    % keep track of filter pars
    filter_pars.max_RelEr=max_RelEr;
    filter_pars.overlap=true;
end

spotData_0=spotData;

spotData=spotData_0(ind);

%  ========================================================
%% RING DETECTION for entire TIME LAPSE

show_fit=0;

% we will use the same params as before

fit_results_all={};
fit_results=[];

frame=1;
fit_results_all{frame}=fit_results_F(:, 2:end);

n_nuc_0=n_nuc;

tic

for frame=2:n_frames
    
    disp(['Running NE ring fitting for frame=', num2str(frame),'...'])
    
    n_nuc=length(spotData);
    
    Fluo_2_ROI=zeros(im_size+2*Rfit);
    Fluo_2_ROI(Rfit+1:end-Rfit, Rfit+1:end-Rfit)=Fluo_2(:, :, frame);
    
    fit_results=zeros(length(spotData),9);
    %coeff_w_int=zeros(length(spotData),18);
    
    for spot=1:length(spotData)
       % parfor spot=1:length(spotData)
        oneSpotData=spotData{spot};
        spot_XYZ=oneSpotData.spotPosition;
        spot_RC=round([spot_XYZ(2),spot_XYZ(1)]);
        
        % get image of the "spot"
        image=Fluo_2_ROI(spot_RC(1)-Rfit+Rfit:spot_RC(1)+Rfit+Rfit, spot_RC(2)-Rfit+Rfit:spot_RC(2)+Rfit+Rfit);
        im_size_0=size(image);
        roi=[0,0; 0,im_size_0(2); im_size_0(1),im_size_0(2); im_size_0(1),0];
        
        % prepare x and y and z values for fitting
        [x_values,y_values,px_values,roi_indx] = extractCellPixels_ROI(image,roi);
        % define fit function:
        model = fittype('a*circle_w_PSF(x, y, sigma, xc, yc, R)+c',...
            'dependent',{'z'},'independent',{'x','y'},...
            'coefficients',{'a','sigma','xc','yc','R','c'});
        % intial estimates for coefficients:
        coeff_0= [max(px_values(:))-min(px_values(:)), 2, im_size_0(1)/2, im_size_0(1)/2,  Rnuc_guess, min(px_values(:))];
        % set bounadries for coefficeints:
        lb=[0, 1, 0.2*im_size_0(1), 0.2*im_size_0(2), 2, 0];
        ub=[Inf,im_size_0(1)/2, 0.8*im_size_0(1), 0.8*im_size_0(2),min(im_size)/2,max(px_values(:))];
        % let's try to fit
        [fit1,gof1,fit_out] = fit([x_values,y_values],px_values,model,'Startpoint',coeff_0,'Lower',lb,'Upper',ub);
        % get coeff values, conf intervals ans smoe goodness of fit measures
        coeff=coeffvalues(fit1);
        a=coeff(1);
        sigma=coeff(2);
        xc=coeff(3);
        yc=coeff(4);
        R=coeff(5);
        c=coeff(6);
        coeff_int=confint(fit1);
        delta_coeff_int=diff(coeff_int,1);
        av_err_sum=sum(delta_coeff_int./coeff)/6;
        adjRsq=gof1.adjrsquare;
        rmse=gof1.adjrsquare;
        
        % collest results:     1       2         3  4  5  6    7            8       9
        %                                   a   sigma    x  y  R  c AvErr  RMSE aR2
        fit_results(spot,:)=[a, sigma, xc+spot_RC(2)-(Rfit+1),  yc+spot_RC(1)-(Rfit+1), R, c, av_err_sum, rmse, adjRsq];
        coeff_w_int_1(spot, :)=coeff;
        coeff_w_int_2(spot, :)=coeff_int(1,:);
        coeff_w_int_3(spot, :)=coeff_int(2,:);
        
        %         if show_fit
        %             if spot==1
        %                 figure; hold off
        %             end
        %             hold off
        %             imshow(image,[min(image(:)),sc(1)*min(image(:))+sc(2)*(max(image(:))-min(image(:)))],'InitialMagnification',1600,'Border','tight');
        %             set(gcf, 'Name', ['frame=', num2str(frame), 'spot #', num2str(spot), ' R = ',num2str(R)]);
        %             hold on
        %             plot((xc-R*cos(phi)), yc-R*sin(phi),':c','LineWidth',2)
        %             plot(xc, yc,'+c','LineWidth',1)
        %             text(xc-R, yc-R-2, ['AvRelEr =', num2str(av_err_sum)],'Color','b')
        %             text(xc+R/2, yc-R-2, ['RMSE=', num2str(rmse)],'Color','b')
        %             text(xc-R, yc+R+2, ['AdjR2 =', num2str(adjRsq)],'Color','b')
        %             pause(0.1)
        %         end
    end
    
    % do the same quality filtering as before
    
    if use_filter
        fit_results_0=fit_results;
        
        ind_0=true(n_nuc,1);
        % Total Rel Error filter
        ind_1=fit_results(:,7)<max_RelEr;
        % overlapping nuclei filter
        x=fit_results(:, 3);
        y=fit_results(:, 4);
        r=fit_results(:, 5);
        
        XX=repmat(x', n_nuc, 1);
        YY=repmat(y, 1, n_nuc);
        RRi=repmat(r', n_nuc, 1);
        RRj=repmat(r, 1, n_nuc);
        
        RR2=(XX-XX').^2+(YY-YY').^2;
        
        ind_2= (RR2>(RRi+RRj).^2);
        ind_2=ind_2+eye(n_nuc);
        [ind_2i, ind_2j]=find(ind_2==0);
        
        ind_2=ind_0;
        ind_2(ind_2i)=false;
        
        % combine all filters in one
        ind=ind_1 & ind_2;
        
        fit_results=[(1:n_nuc)', fit_results];
        
        fit_results_F=fit_results(ind,:);
        
        coeff_w_int_1=coeff_w_int_1(ind, :);
        coeff_w_int_2=coeff_w_int_2(ind, :);
        coeff_w_int_3=coeff_w_int_3(ind, :);
    
        spotData=spotData(ind);
    end
    
    n_spots=length(spotData);
    coeff_w_int(end+1:end+n_spots, 1:3:18)=coeff_w_int_1;
    coeff_w_int(end+1:end+n_spots, 2:3:18)=coeff_w_int_2;
    coeff_w_int(end+1:end+n_spots, 3:3:18)=coeff_w_int_3;
    
    fit_results_all{frame}=fit_results_F;
    fit_results_all_no_filter{frame}=fit_results_0;
    coeff_w_int_all{frame}=coeff_w_int;
    
    
    disp(['Done with frame=', num2str(2),'; found ', num2str(length(spotData)), ' good NE fits'])
    toc
    
end

% keep track of filter pars
    filter_pars.max_RelEr=max_RelEr;
    filter_pars.overlap=true;
    
disp('DONE with NE detection in all frames')



%=========================================================================
%% SAVE SPOTS and detection params
% % %% ???? save the most important variables
% % % save
save([path_save,file_save,file_suffix,'.mat'],'path', 'filename', 'dx', 'spotDetection_params', 'filter_pars', 'fit_results_all_no_filter', 'fit_results_all','spotData','spotData_0')
disp(['Data saved into ',[path_save,file_save,file_suffix,'.mat']])
