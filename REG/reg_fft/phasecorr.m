function d = phasecorr(A,B)
%PHASECORR Compute phase correlation matrix
%
%   D = PHASECORR(A, B) computes the phase correlation matrix from input
%   2-D images A and B. PHASECORR returns the 'full' correlation matrix
%   similar to normxcorr2.

%   Copyright 2013 The MathWorks, Inc.

size_A  = size(A);
size_B  = size(B);

% Let fft2 zero pad time domain signals such that we form the 'full'
% correlation matrix analogous to the result produced by normxcorr2.
outSize = size_A + size_B - 1;

% Form 2-D spectra of moving and fixed.
% Window each signal prior to taking FFT to help reduce ringing effects in
% frequency domain due to finite image.
A = fft2(A,outSize(1),outSize(2));
B = fft2(B,outSize(1),outSize(2));

% Form phase correlation matrix, d
% Use 'symmetric' option as performance optimization. We expect that the
% input moving images A and B are real valued, so d should always be real
% valued.
ABConj = A .* conj(B);
d = ifft2(ABConj ./ abs(eps+ABConj),'symmetric');



