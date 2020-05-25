function [P] = unmix2(data, endmemberMean)
% This function performs NCM-QP unmixing (Normal Compositional Model, Quadratic Programming approach).
% INPUTS 
%  data - double Mat - DxN matrix of image data.
%  endmemberMean - double Mat - DxM matrix of endmember mean value (known).
% OUTPUTS
%   P - double Mat - MxN matrix of abundances/proportions corresponding to N input
%       pixels and M endmembers
%
% This product is Copyright (c) 2014 University of Missouri
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
%
%   1. Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%   2. Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in the
%      documentation and/or other materials provided with the distribution.
%   3. Neither the name of the University nor the names of its contributors
%      may be used to endorse or promote products derived from this software
%      without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY OF MISSOURI AND
% CONTRIBUTORS ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,
% INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED.  IN NO EVENT SHALL THE UNIVERSITY OR CONTRIBUTORS
% BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES,
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
% HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
% OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

options = optimset('Display', 'off', 'LargeScale', 'off');
warning off all

%endmemberMean should be column vectors
X = data;
%number of endmemberMean
M = size(endmemberMean, 2);
%number of pixels
N = size(X, 2);

%set up constraint matrices
A1 = eye(M);
A1 = -1*A1;
b1 = zeros([M, 1]);

A2 = ones([1, M]);
b2 = 1;

A3 = -1*ones([1, M]);
b3 = -1;

A = vertcat(A1, A2);
A = vertcat(A, A3);

b = vertcat(b1, b2);
b = vertcat(b, b3);

H = (2*endmemberMean'*endmemberMean);

parfor i = 1:N
    Xi = X(:,i); 
    F = (-2*Xi'*endmemberMean)';
    P(i,:) = quadprog(H, F, A, b, [], [], [], [], [], options);%QP unmixing
end

P(P<0) = 0;

end

                