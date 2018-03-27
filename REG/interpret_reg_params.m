function [I1,wl1,rgb1,g,H1,extra,s,G,S] = interpret_reg_params(...
    I,wl,bbl,rgb,U2,V2,rho,s2,T2,sigma2)
%INTERPRET_REG_PARAMS remove bad bands and warp the hyperspectral image,
%rigidly transform the color image, and estimate the PSF and SRF.
% Input:
%   I,wl,bbl - hyperspectral image, wavelength, and bad band list
%   rgb - color image
%   U2,V2,rho,s2,T2,sigma2 - registration parameters
% Output:
%   I1 - bad-band-removed and warped hyperspectral image
%   wl1 - bad-band-removed wavelength
%   rgb1 - registered hi-resolution color image
%   g - R by 1 PSF
%   H1 - (d+1) by b SRF
%   extra,s - scale (sx,sy) and margin (extra_x,extra_y)
%   G - PSF matrix
%   S - selection matrix
%
%% get bad-band-removed and warped hyperspectral image
options = [];
options.isRigid = 0;
options.useTranslationField = 1;
options.U = U2;
options.V = V2;
options.rho = rho;

I(:,:,~bbl) = [];
I1 = transform(I,eye(3),[1,1],1e-3,size(I,2),size(I,1),options);
wl1 = wl(bbl);

%% get registered color image
[rows,cols,B] = size(I1);
s = s2;
T = T2;
sigma = sigma2;

s1 = ceil(s(1)); % scale_x
s2 = ceil(s(2)); % scale_y

options.isRigid = 1;
[rgb_small,rgb1,g,extra] = transform(double(rgb),T,s,sigma,cols,rows,options);

%% calculate SRF
% No bad bands available in the visible range
[I_sel,bands_sel] = select_relevant_bands(I1,wl1,'color');
wl_sel = wl1(bands_sel)';
S = zeros(B,length(wl_sel));
for i = 1:size(S,2)
    S(bands_sel(i),i) = 1;
end

Y1 = reshape_hsi(I_sel);
Y1 = [ones(size(Y1,1),1),Y1];

X = reshape_hsi(rgb_small);
options = [];
options.lambda = 1e-3;

H1 = solve_for_H1(X, Y1, options);

%% calculate G 
s = [s1,s2];
Cs = construct_C(I1,rgb1,s,extra);
C = cat(2,Cs{:});
g = g(:);
G = g2G(C,g,rows*cols);

end