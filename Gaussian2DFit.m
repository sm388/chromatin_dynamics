function [parameters2, mxy] = Gaussian2DFit(img2,centroid,initpar,scale,showplot)

cx = centroid(1); 
cy = centroid(2);

% fit 2D gaussian
mxy = double(img2(round(cy)-scale:round(cy)+scale,round(cx)-scale:round(cx)+scale));%for scale=5: 11*11 matrix of pixel values around centroid position
% max_y = size(img2,1);
% max_x = size(img2,2);
% mxy = double(img2(max(1,round(cy)-scale):min(max_y,round(cy)+scale),max(1,round(cx)-scale):min(max_x,round(cx)+scale)));
[sizey,sizex] = size(mxy);%11*11
[X,Y] = meshgrid(1:sizey,1:sizex);
%options = optimset('TolX',1e-8,'Display','Off','MaxIter',1e4,'MaxFunEvals',1e3,'Jacobian','on');
options = optimset('TolX',1e-8,'Display','Off','MaxIter',1e4,'MaxFunEvals',1e3);

%original
lowpar = double([scale-.5,scale-.5,0,0,0]);
highpar = double([scale+2.5,scale+2.5,8,65535,65535]);


% lowpar = double([scale-3,scale-3,0,0,0]);
% highpar = double([scale+3,scale+3,8,65535,65535]);

%ORIGINAL CODE
% Gauss2D = inline('(Gaussian2D(p,x,y) - z)','p','x','y','z');
% [parameters,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = ...
%         lsqnonlin(Gauss2D,initpar,lowpar,highpar,options,X,Y,mxy);

%Gauss2D = 
% [parameters,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = ...
%         lsqnonlin(@(p) Gauss2D(p,X,Y,mxy),initpar,lowpar,highpar,options);
%     
%     [parameters,fval,exitflag,output,lambda,grad,hessian] = ...
%         fmincon(@(p) Gauss2D(p,X,Y,mxy),initpar,[],[],[],[],lowpar,highpar,[],options);
%     err=sqrt(diag(inv(hessian)));
    [parameters,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = ...
        lsqnonlin(@(p) Gauss2D(p,X,Y,mxy),initpar,lowpar,highpar,options);
    
if showplot == 1 
    figure(234234);clf
    subplot(2,1,1); surf(mxy); 
    subplot(2,1,2);
    contour3(mxy,'--');
    colormap cool
    hold on;
    contour3(Gaussian2D(parameters,X,Y),'-');
end

%include position within mxy
parameters2(6)=parameters(1);
parameters2(7)=parameters(2);


%orig
parameters2(1) = parameters(1) + cx - scale - 1;
parameters2(2) = parameters(2) + cy - scale - 1;

% mlb edit 12/1/17
% x
% parameters2(1) = parameters(2) + cx - scale - 1;
% y
% parameters2(2) = parameters(1) + cy - scale - 1;


parameters2(3) = parameters(3);
parameters2(4) = parameters(4);
parameters2(5) = parameters(5);

end