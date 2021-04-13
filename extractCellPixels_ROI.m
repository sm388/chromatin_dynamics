function [x_values,y_values,px_values,roi_indx] = extractCellPixels_ROI(image,roi)
%--------------------------------------------------------------------------
% returns vectors of X,Y and pixel values for the image ROIs 
%--------------------------------------------------------------------------
%
%**********INPUT********:
% image = image 
% roi = ROI in Matlab sense, nx2 matrix with x znd y coordinates of n vertices of closed ROI 
%
%*********OUTPUT********:
% x_values: x-coordiantes of the pixels within ROI in Image.
% y_values: y-coordiantes of the pixels within ROI in Image.
% px_values: values of the pixels within ROI in Image.
% roi: coordinates (in image coordinates) of ROI
% roi_indx: linear indices of the pixels within ROI in Image.
%
%@author:  Ivan Surovtsev
%@date:    March 31, 2017

im_size=size(image);
roi=double(roi);

xx=repmat(1:im_size(2),im_size(1),1);
 yy=repmat((1:im_size(1))',1,im_size(2));

% make a mask over ROI
msk = poly2mask(roi(:,1),roi(:,2),im_size(1),im_size(2));


% getting indices and values
lin_indx=1:length(image(:));
 roi_indx=lin_indx(msk==1);
px_values=double(image(msk==1));
 x_values=double(xx(msk==1));
 y_values=double(yy(msk==1));

end