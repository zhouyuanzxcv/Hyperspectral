function [tform, degree, s, t] = reg_hyper_fft(I1, I, wl, transformation, windowing)
%REG_HYPER_FFT Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
    transformation = 'similarity';
end
if nargin < 5
    windowing = 1;
end

rgb_ind = find_rgb_ind(wl);
bands_sel = false(1,length(wl));
for i = 1:length(rgb_ind)
    bands_sel(rgb_ind(i)-2:2:rgb_ind(i)+2) = 1;
end

I_sel = I(:,:,bands_sel);

theta = 0;
s = 1;

switch transformation
    case 'translation'
        [tform,peak] = findTranslation(I1, I_sel, windowing);
    case 'similarity'
        [tform,peak,theta,s] = findSimilarity(I1, I_sel, windowing);
    case 'rigid'
        [tform,peak,theta,s] = findRigid(I1, I_sel, windowing);
    otherwise
end
degree = theta / (pi/180);

T = tform2myT(tform);
t = T([1,2],3);

end

function T = tform2myT(tform)
T = tform.T';
T = inv(T);
end

%--------------------------------------------------------------
function [tform,peak] = findTranslation(moving,fixed,windowing)

% moving = manageWindowing(moving,windowing);
% fixed  = manageWindowing(fixed,windowing);

[vec,peak] = findTranslationPhaseCorr(moving,fixed);
tform = affine2d([1, 0, 0; 0, 1, 0; vec(1), vec(2), 1]);

end

%------------------------------------------------------------------------
function [tform,peak,theta,s] = findRigid(moving,fixed,windowing)
% A nice block diagram of the pure rigid algorithm appears in:
%   Y Keller, "Pseduo-polar based estimation of large translations rotations and
%   scalings in images", Application of Computer Vision, 2005. WACV/MOTIONS
%   2005 Volume 1. 
%
% This follows directly from the derivation in Reddy, Chatterji.

% Move Moving and Fixed into frequency domain
[M,F] = getFourierMellinSpectra(moving,fixed,windowing);

thetaRange = [0 pi];
Fpolar = images.internal.Polar(F,thetaRange);
Mpolar = images.internal.Polar(M,thetaRange);

Fpolar.resampledImage = manageWindowing(Fpolar.resampledImage,windowing);
Mpolar.resampledImage = manageWindowing(Mpolar.resampledImage,windowing);

% Solve a 1-D phase correlation problem to resolve theta. We already know
% scale. Choose a 1-D profile in our Polar FFT grid parallel to the theta axis.
numSamplesRho = size(Fpolar.resampledImage,1);
rhoCenter = round(0.5+numSamplesRho/2);
vec = findTranslationPhaseCorr(Mpolar.resampledImage(rhoCenter,:,:),Fpolar.resampledImage(rhoCenter,:,:));

% Translation vector is zero based. We want to translate vector
% into one based intrinsic coordinates within the polar grid.
thetaIntrinsic = abs(vec(1))+1;
% We passed a vector to findTranslationPhaseCorr;
rhoIntrinsic   = 1;

% The translation vector implies intrinsic coordinates in the
% Fixed/Moving log-polar grid. We want to convert these intrinsic
% coordinate locations into world coordinates that tell us
% rho/theta.
[theta,~] = intrinsicToWorld(Fpolar,thetaIntrinsic,rhoIntrinsic);

% Use sign of correlation offset to figure out whether rotation
% is positive or negative.
theta = -sign(vec(1))*theta;

% By definition, Scale is 1 for a rigid transformation.
s = 1;

[tform,peak,theta] = solveForTranslationGivenScaleAndRotation(moving,fixed,s,theta,windowing);

end

%----------------------------------------------------------------
function [tform,peak,theta,s] = findSimilarity(moving,fixed,windowing)

% Move Moving and Fixed into frequency domain
[M,F] = getFourierMellinSpectra(moving,fixed,windowing);

% (Reddy,Chatterji) recommends taking advantage of the conjugate
% symmetry of the Fourier-Mellin spectra. All of the unique
% spectral information is in the interval [0,pi].
thetaRange = [0 pi];
Fpolar = images.internal.LogPolar(F,thetaRange);
Mpolar = images.internal.LogPolar(M,thetaRange);

% Use phase-correlation to determine the translation within the
% log-polar resampled Fourier-Mellin spectra that aligns moving
% with fixed.
Fpolar.resampledImage = manageWindowing(Fpolar.resampledImage,windowing);
Mpolar.resampledImage = manageWindowing(Mpolar.resampledImage,windowing);

% Obtain full phase correlation matrix
d = phasecorr2(Mpolar.resampledImage, Fpolar.resampledImage);

% Constrain our search in D to the range 1/4 < S < 4.
d = suppressCorrelationOutsideScale(d,Fpolar,4);

% Find the translation vector in log-polar space.
vec = findTranslationPhaseCorr(d);

% Translation vector is zero based. We want to translate vector
% into one based intrinsic coordinates within the log-polar grid.
thetaIntrinsic = abs(vec(1))+1;
rhoIntrinsic   = abs(vec(2))+1;

% The translation vector implies intrinsic coordinates in the
% Fixed/Moving log-polar grid. We want to convert these intrinsic
% coordinate locations into world coordinates that tell us
% rho/theta.
[theta,rho] = intrinsicToWorld(Fpolar,thetaIntrinsic,rhoIntrinsic);

% Use sign of correlation offset to figure out whether rotation
% is positive or negative.
theta = -sign(vec(1))*theta;

% Use sign of correlation offset to figure out whether or not to invert scale factor
s = rho .^ -sign(vec(2));

[tform,peak,theta] = solveForTranslationGivenScaleAndRotation(moving,fixed,s,theta,windowing);

end


%----------------------------------------------------------------
function [tform,peak,theta] = solveForTranslationGivenScaleAndRotation(moving,fixed,S,theta,windowing)
% There is a 180 degree ambiguity in theta solved in R,Theta space. This
% ambiguity stems from the conjugate symmetry of the Fourier spectrum for real
% valued input images.
%
% This function resolves the ambiguity by forming two resampled versions of moving
% rotated by theta, theta+180, phase correlating each version of the
% resampled image with fixed, and choose the S,Theta that has the highest
% final peak correlation during recovery of translation.
%
% We save 1 FFT2 operation at full scale with the following
% optimizations:
%
% 1) By directly performing the phase correlation here instead of calling
% phasecorr/findTranslationPhaseCorr directly, we save 1 FFT operation by
% not computing the spectra of fixed twice.

theta1 = theta;
theta2 = theta+pi;

tform1 = affine2d([S.*cos(theta1) -S.*sin(theta1) 0; S.*sin(theta1) S.*cos(theta1) 0; 0 0 1]);
tform2 = affine2d([S.*cos(theta2) -S.*sin(theta2) 0; S.*sin(theta2) S.*cos(theta2) 0; 0 0 1]);

[scaledRotatedMoving1,RrotatedScaled1] = imwarp(moving,tform1,'SmoothEdges', true);

% scaledRotatedMoving1 = manageWindowing(scaledRotatedMoving1,windowing);

% This step is equivalent to: 
%   [scaledRotatedMoving2,RrotatedScaled2] = imwarp(moving,tform2)
% We do this to gain efficiency in computing scaledRotatedMoving2,
scaledRotatedMoving2 = rot90(scaledRotatedMoving1,2);
RrotatedScaled2 = imref2d(size(scaledRotatedMoving1),...
                          sort(-RrotatedScaled1.XWorldLimits),...
                          sort(-RrotatedScaled1.YWorldLimits));

% Form 2-D spectra associated with scaledRotatedMoving1, scaledRotatedMoving2, and fixed.
d1 = phasecorr2(scaledRotatedMoving1, fixed);
d2 = phasecorr2(scaledRotatedMoving2, fixed);

% size_moving  = size(scaledRotatedMoving1);
% size_fixed  = size(fixed);
% outSize = size_moving + size_fixed - 1;
% M1 = fft2(scaledRotatedMoving1,outSize(1),outSize(2));
% F  = fft2(fixed,outSize(1),outSize(2));
% M2 = fft2(scaledRotatedMoving2,outSize(1),outSize(2));
% 
% % Form the phase correlation matrix d1 for M1 correlated with F.
% ABConj = F .* conj(M1);
% d1 = ifft2(ABConj ./ abs(eps+ABConj),'symmetric');
% 
% % Form the phase correlation matrix d2 for M2 correlated with F.
% ABConj = F .* conj(M2);
% d2 = ifft2(ABConj ./ abs(eps+ABConj),'symmetric');

% Find the translation vector that aligns scaledRotatedMoving1 with fixed and
% scaledRotatedMoving2 with fixed. Choose S,theta,translation estimate that has
% the highest peak correlation in the final translation recovery step.
[vec1,peak1] = findTranslationPhaseCorr(d1);
[vec2,peak2] = findTranslationPhaseCorr(d2);

if peak1 >= peak2
    vec = vec1;
    tform = tform1;
    RrotatedScaled = RrotatedScaled1;
    peak = peak1;
    theta = theta1;
else
    vec = vec2;
    tform = tform2;
    RrotatedScaled = RrotatedScaled2;
    peak = peak2;
    theta = theta2;
end

% The scale/rotation operation performed prior to the final
% phase-correlation step results in a translation. The translation added
% during scaling/rotation is defined by RrotatedScaled. Form the final
% effective translation by summing the translation added during
% rotation/scale to the translation recovered in the final translation
% step.
finalXOffset  = vec(1) + (RrotatedScaled.XIntrinsicLimits(1)-RrotatedScaled.XWorldLimits(1));
finalYOffset  = vec(2) + (RrotatedScaled.YIntrinsicLimits(1)-RrotatedScaled.YWorldLimits(1));

tform.T(3,1:2) = [finalXOffset, finalYOffset];
% t = [finalXOffset, finalYOffset];

end


%-----------------------------------------------------------------------
function [M,F] = getFourierMellinSpectra(moving,fixed,windowing)

% Move Moving and Fixed into frequency domain
M_size = size(moving);
F_size = size(fixed);
outsize = M_size(1:2) + F_size(1:2) - 1;

% add the constant function
fixed = cat(3, ones(size(fixed,1),size(fixed,2)), fixed);

% Apply windowing function to moving and fixed to reduce aliasing in
% frequency domain.
moving = manageWindowing(moving,windowing);
fixed  = manageWindowing(fixed,windowing);

% Obtain the spectra of moving and fixed: M and F.
M = fft2(moving,outsize(1),outsize(2));
F = fft2(fixed,outsize(1),outsize(2));

% Shift DC of fft to center
F = fft2shift(F);
M = fft2shift(M);

% Form Magnitude Spectra
F1 = [];
cs = combnk(1:size(F,3),2);
for i = 1:size(F,3)
    F1 = cat(3, F1, abs(F(:,:,i)).^2);
end
for i = 1:size(cs,1)
    F1 = cat(3, F1, real(F(:,:,cs(i,1)) .* conj(F(:,:,cs(i,2)))));
end

F = F1;
M = abs(M).^2;

% Apply High-Pass Emphasis filter to each image (Reddy, Chatterji)
H = createHighPassEmphasisFilter(outsize);
H = H.^2;

F = F .* repmat(H, [1,1,size(F,3)]);
M = M .* repmat(H, [1,1,size(M,3)]);

end

function M1 = fft2shift(M)
M1 = [];
for k = 1:size(M,3)
    M1 = cat(3, M1, fftshift(M(:,:,k)));
end
end

%-----------------------------------------------------------
function d = suppressCorrelationOutsideScale(d,Fpolar,scale)
% This function takes a phase correlation matrix that relates the same
% sized log-polar grids Fpolar and Mpolar. We return a phase correlation
% matrix in which we set regions of the phase correlation matrix outside
% the symmetric range (1/scale, scale) to -Inf. This allows us to limit the
% search space of the phase correlation matrix during peak detection so
% that we will never find peaks that correspond to a scale value outside of
% the limits of scale.

[~,logRhoIndex] = worldToIntrinsic(Fpolar,0,scale);
logRhoIndex = floor(logRhoIndex);

% Create mask that is false where S is outside the range (1/scale,scale).
phaseCorrMask = false(size(d));
phaseCorrMask((logRhoIndex+1):(end-logRhoIndex+1),:) = true;

% Constrain our search in D to the range 1/scale < S < scale.
d(phaseCorrMask) = 0;

end

%-----------------------------------------------------------
function [vec,peakVal] = findTranslationPhaseCorr(varargin) %#codegen
%findTranslationPhaseCorr Determine translation using phase correlation.
%
%   [vec,peakVal] = findTranslationPhaseCorr(MOVING, FIXED) estimates the
%   translation of MOVING necessary to align MOVING with the
%   fixed image FIXED. The output VEC is a two element vector of the form
%   [deltaX, deltaY]. The scalar peakVal is the peak value of the phase
%   correlation matrix used to estimate translation.
%
%   [vec,peakVal] = findTranslationPhaseCorr(D) estimates the translation
%   of MOVING necesary to align MOVING with the fixed image FIXED. D is a
%   phase correlation matrix of the form returned by:
%
%       D = phasecorr(fixed,moving).

%   Copyright 2013 The MathWorks, Inc.
%
% Modified by Yuan to use the new phase correlation function
narginchk(1,2)

if nargin == 1
    d = varargin{1};
else
    moving = varargin{1};
    fixed  = varargin{2};
    % Compute phase correlation matrix, D
    d = phasecorr2(moving,fixed);
end

% Use simple global maximum peak finding. Surface fit using 3x3
% neighborhood to refine xpeak,ypeak location to sub-pixel accuracy.
subpixel = true;
[xpeak,ypeak,peakVal] = findpeak(d,subpixel);

% findpeak returns 1 based MATLAB indices. We want 0 based offset for
% translation vector.
xpeak = xpeak-1;
ypeak = ypeak-1;

outSize = size(d);

% Covert peak locations in phase correlation matrix to translation vector
% that defines translation necessary to align moving with fixed.
%
% The translation offsets implied by the phase correlation matrix have the form:
%
% [0, 1, 2,...-2, -1];
%
% The logic below figures out whether we are past the center region of the
% phase correlation matrix in which the offset signs switch from positive to
% negative, i.e. are we closer to the right edge or the left edge?
if xpeak > abs(xpeak-outSize(2))
    xpeak = xpeak-outSize(2);
end

if ypeak > abs(ypeak-outSize(1))
    ypeak = ypeak-outSize(1);
end

% Ensure that we consistently return double for the offset vector and for
% the peak correlation value.
vec = double([xpeak, ypeak]);
peakVal = double(peakVal);
end

%--------------------------------------------
function img = manageWindowing(img,windowing)

if windowing
    img = img .* repmat(createBlackmanWindow(size(img)), [1,1,size(img,3)]);
end

end

%--------------------------------------------
function h = createBlackmanWindow(windowSize)
% Define Blackman window to reduce finite image replication effects in
% frequency domain. Blackman window is recommended in (Stone, Tao,
% McGuire, Analysis of image registration noise due to rotationally
% dependent aliasing).

M = windowSize(1);
N = windowSize(2);

a0 = 7938/18608;
a1 = 9240/18608;
a2 = 1430/18608;

n = 1:N;
m = 1:M;

% Make outer product degenerate if M or N is equal to 1.
h1 = 1;
h2 = 1;
if M > 1
    h1 = a0 - a1*cos(2*pi*m / (M-1)) + a2*cos(4*pi*m / (M-1));
end
if N > 1
    h2 = a0 - a1*cos(2*pi*n / (N-1)) + a2*cos(4*pi*n / (N-1));
end

h = h1' * h2;

end

%---------------------------------------------
function H = createHighPassEmphasisFilter(outsize)
% Defines High-Pass emphasis filter used in Reddy and Chatterji

numRows = outsize(1);
numCols = outsize(2);

x = linspace(-0.5,0.5,numCols);
y = linspace(-0.5,0.5,numRows);

[x,y] = meshgrid(x,y);

X = cos(pi*x).*cos(pi*y);

H = (1-X).*(2-X);

end
