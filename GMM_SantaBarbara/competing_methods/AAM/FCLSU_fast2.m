function [out,options] = FCLSU_fast2(HIM,M,options)
% Fully Constrained Linear Spectral Unmixing
% Perform a Linear least squares with nonnegativity constraints.
% --------------------------------------------------------------------
% Input:   HIM : input data [nrows x nchannels]
%          M   : set of p endmembers [nchannels x p].
%
% Output:  out : fractions [nrows x p]
%
%
% Copyright (2007) GRNPS group @ University of Extremadura, Spain.
%
% *** Edited version ***


[ns,nb] = size(HIM);
[l,p] = size(M);

Delta = 1/1000; % should be an small value

N = zeros(l+1,p);
N(1:l,1:p) = Delta*M;
N(l+1,:) = ones(1,p);
s = zeros(l+1,1);

OutputImage = zeros(ns,p);

go=(nargin==3);

for i = 1:ns
    s(1:l) = Delta*HIM(i,:)';
    s(l+1) = 1;
    if go==0
        [Abundances,options] = lsqnonneg_fast(N,s);
        go=1;
    else
        Abundances = lsqnonneg_fast(N,s,options);
    end
    OutputImage(i,:) = Abundances;
end


out = OutputImage;