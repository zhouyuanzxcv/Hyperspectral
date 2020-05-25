function pixelList = EF_reshape(imageData)

% Inputs:
%   imageData - hyperspectral data cube, M by N by d
%
% Outputs:
%   pixelList - linear data, d by (M*N)

imageData=double(imageData);

pixelList = reshape(shiftdim(imageData(:,:,:),2),size(imageData,3),size(imageData,1)*size(imageData,2));