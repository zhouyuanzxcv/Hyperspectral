function [P2] = unmix_qpas_correct(data, endmembers,methodFlag)
% This function performs quadratic unmixing given endmember mean values
% SYNTAX : [P2] = unmix_qpas_correct(data, endmembers)
% INPUTS 
%     data - double Mat - DxN matrix of N data points of dimensionality D 
%         (i.e.  N pixels with D spectral bands, each pixel is a column vector).
%     endmembers - double Mat - DxM matrix of M endmembers (endmember mean spectra) of dimensionality D.
%     methodFlag - double Mat - 1x1 integer. 
%               If methodFlag == 1, BCM Spectral QP unmixing;
%               If methodFlag == 3, BCM Spatial QP unmixing.
% OUTPUTS
%     P2 - double Mat - NxM matrix of abundances/proportions corresponding to N input
%         pixels and M endmembers
%
% This product is Copyright (c) 2014 University of Missouri.
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

warning('off', 'all');
options = optimset('Display', 'off', 'LargeScale', 'off');

%endmembers should be column vectors
X = data;

%number of endmembers
M = size(endmembers, 2);
%number of pixels
N = size(X, 2);


%Equation constraint Aeq*x = beq
%All values must sum to 1 (X1+X2+...+XM = 1)
Aeq = ones([1, M]);
beq = 1;

%Boundary Constraints lb >= x >= ub
%All values must be greater than 0 (0 ? X1,0 ? X2,...,0 ? XM)
lb = zeros([M, 1]);
ub = ones([M,1]);

H = 2*(endmembers'*endmembers);

P2 = zeros(N,M);

constraintErrorTrip = 0.01; %off by 0.01% indicates error
%F = ((-2*X'*endmembers)+repmat(gammaVecs,N,1))';

if methodFlag == 1
    h = waitbar(0,'Stage 2/2: QP unmixing...','Name','BCM Spectral Unmixing');
elseif methodFlag == 3
    h = waitbar(0,'Stage 2/2: QP unmixing...','Name','BCM Spatial Unmixing');  
else
    error('not yet implemented');
    P2 = [];
end    
% parfor i = 1:N

for i = 1:N
    F = (-2*X(:,i)'*endmembers)';
%    qpas_ans = qpas(H, F, [], [], Aeq, beq, lb, ub, 0);
    qpas_ans =  quadprog(H, F, [], [], Aeq, beq, lb, [], [], options);
    qpas_sum = sum(qpas_ans);
    constraintPercentError = abs(1-qpas_sum)*100;
    constraintError = (constraintPercentError > constraintErrorTrip); %set if constraints are not met
    if (constraintError)
        %disp('constraint error detected and corrected');
        P2(i,:) = quadprog(H, F, [], [], Aeq, beq, lb, [], [], options);
    else 
        P2(i,:) = qpas_ans;
    end
    
    perc = round(i/N*100);
    if (mod(perc, 20) == 0)
        waitbar(i/N,h,sprintf('Stage 2/2: QP unmixing...%d%%',perc));
    end
end

P2(P2<0) = 0;
close(h);
end

