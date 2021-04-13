
function I_xy=circle_w_PSF(x, y, sigma, xc, yc, R)
%--------------------------------------------------------------------------
% returns I(x,y) - light intensity  at (x, y) point for the image of the
% circular (ring) source of radius R with the center at (xc, yc) , given
% a point-spread-function with the width sigma. I.e., calulates convolution
% of the ring with intensity 1 and point-spread function (PSF). PSF is approximated as a Gaussian. 
%--------------------------------------------------------------------------
%
%**********INPUT********:
 % (x, y)  - coordinates of the image
% sigma - width of PSF 
% (xc, yc) - coordinates of the ring center
% R - radius of the ring
%*********OUTPUT********:
% I_xy - Intensity at the (x,y) point of the image 
%
%@author:  Ivan Surovtsev
%@date:    2018.12.14
        
        circle_w_PSF = @ (phi,x,y,xc,yc,R,sigma) exp(-((x-xc-R*cos(phi)).^2+(y-yc-R*sin(phi)).^2)/(2*sigma^2));
        
        I_xy=integral(@(phi) circle_w_PSF(phi,x,y,xc,yc,R,sigma),0,2*pi,'ArrayValued',true);
     
end