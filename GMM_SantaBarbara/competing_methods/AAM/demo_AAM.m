% demo_AAM

% This demo shows how to use the AAM and MESMA code, with a quick example
% on Gaussian artificial data sets.

% Initialization
d=40;         % number of spectral bands
num=1000;     % number of pixels in input data set
p=3;          % number of spectral libraries
N=[10,20,30]; % sizes of spectral libraries

% Generate data sets
L{1}=randn(d,N(1)); 
L{2}=randn(d,N(2));
L{3}=randn(d,N(3));
x=randn(d,num);

% Call AAM 
[index1, abundances1, reconstruction1, error1]=AAM(x,L,3);

% Call MESMA
[index2, abundances2, reconstruction2, error2]=MESMA_bruteforce (x,L);
