function [P, F, E, S] = BCM_Spectral(Xim, Parameters, F, E, S)
% This function performs BCM-Spectral unmixing (Beta Compositional Model).
%    Note that this code only deals with input data points defined on [0,1]!
%
% SYNTAX : [P, F, E, S] = BCM_Spectral(Xim, Parameters, F, E, S)
%          or  [P, F, E, S] = BCM_Spectral(Xim, Parameters)
% INPUTS 
%     Xim - double Mat - N1xN2xD matrix of N1xN2 image data points of dimensionality D 
%           (i.e.  N=N1xN2 pixels with D spectral bands). Defined on [0,1].
%     Parameters:    Parameters.
%     F,E,S:         If not given(if nargin<3), they will be computed in this code.
% OUTPUTS
%     P - double Mat - NxM matrix of abundances/proportions corresponding to N input
%         pixels and M endmembers
%     F - double Mat - NxD matrix, Ratio of single Beta parameters
%     E - double Mat - NxD matrix, Expected value of the data (neighbors)              
%     S - double Mat - NxD matrix, Variance of the data (neighbors)
%
% REFERENCE
% X. Du, A. Zare, P. Gader, D. Dranishnikov, 
% “Spatial and Spectral Unmixing Using the Beta Compositional Model,”  
% IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 
% vol. 7, no. 6, pp. 1994-2003, June 2014.
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

Parameters.N = size(Xim,1)*size(Xim,2);
Parameters.M = Parameters.MaxNumEMs;
Parameters.D = size(Xim,3);

Xlist = reshape(Xim, [size(Xim,1)*size(Xim,2), size(Xim,3)]);

% Compute spectral distances and find K spectral neighbors
% Modified by Yuan on 5.17
if 0
    dists = pdist2(Xlist, Xlist);
    for i = 1:Parameters.N
        [~, s] = sort(dists(i,:));
        X(i,:,:) = Xlist(s(1:Parameters.K), :)';
    end
else
    X = zeros(Parameters.N, size(Xlist,2), Parameters.K);
    [~,dist_inds] = pdist2(Xlist, Xlist, 'euclidean', 'Smallest', Parameters.K);
    for i = 1:Parameters.N
        X(i,:,:) = Xlist(dist_inds(:,i), :)';
    end
end

for m = 1:Parameters.M
    M(:,m) = Parameters.BetaParameters(:,1,m)./(Parameters.BetaParameters(:,1,m) + Parameters.BetaParameters(:,2,m));
    if(Parameters.methodFlag == 2)
        V(:,m) = (Parameters.BetaParameters(:,1,m).*Parameters.BetaParameters(:,2,m))./(((Parameters.BetaParameters(:,1,m)+Parameters.BetaParameters(:,2,m)).^2).*(Parameters.BetaParameters(:,1,m)+Parameters.BetaParameters(:,2,m)+1));
    end
end

 h = waitbar(0,'Stage 1/2: Estimate Beta parameters...','Name','BCM Spectral Unmixing');
if(nargin < 3)
    for n = 1:size(X, 1)
        %first, estimate Beta distribution parameters e and f from data.
        for d = 1:Parameters.D
            % Added by Yuan on 5.17.2017
            xmin = min(X(n,d,:));
            xmax = max(X(n,d,:));
            if abs(xmin - xmax) <= 1e-6 % trigger error in betafit
                X(n,d,:) = X(n,d,:) + rand(1,1,Parameters.K)*(10.^-4);
            end
            
            for xx = 1:Parameters.K
                if (X(n,d,xx)<=0)
                    X(n,d,xx)=0+rand(1)*(10.^-3); % Make sure that not all neighbors are 0 to fit beta
                elseif (X(n,d,xx)>=1)
                    X(n,d,xx)=1-rand(1)*(10.^-3); % Make sure that not all neighbors are 1 to fit beta
                end
            end
            ef(d,:) = betafit(squeeze(X(n,d,:)));
        end
        
        F(n,:) = ef(:,1)./ef(:,2);
        E(n,:) = 1./(1./F(n,:) + 1);
        S(n,:) = (ef(:,2)' + 1./(1 + F(n,:))).*(((1-F(n,:)).^3)./F(n,:));
        
        %Display Progress bar
        perc = round(n/size(X, 1)*100);
        waitbar(n/size(X, 1),h,sprintf('Stage 1/2: Estimate Beta parameters...%d%%',perc));
%         n
    end
end
    
close(h);

if(Parameters.methodFlag == 1)  % BCM-Spectral-QP approach
    
    P = unmix_qpas_correct(E', M, Parameters.methodFlag);
    
elseif(Parameters.methodFlag == 2)  % BCM-Spectral-MH approach
 
    %Initialize P
    P = DirichletSample(ones(Parameters.N, Parameters.M));
    

    %Evaluate Samples
    Term1 = (-1/2)*(E - P*M').^2;
    Term2 = (-1/2)*(S - P*V').^2;
    LnLikelihoodNew = sum(Term1 + Term2,2);
    Pold = P;
    lTrackold = LnLikelihoodNew;
    
    h = waitbar(0,'Stage 2/2: MH Sampling...','Name','BCM Spectral Unmixing');
    %Iterate and Sample
    for iter = 2:Parameters.NumberIterations
         
        %Sample P
        [psamples] = DirichletSample(ones(Parameters.N,Parameters.M));
        Term1 = (-1/(2*Parameters.sigM))*(E - psamples*M').^2;
        Term2 = (-1/(2*Parameters.sigV))*(S - psamples*V').^2;
        LnLikelihoodNew = sum(Term1 + Term2,2);
        
        %Evaluate Sample
        acceptRatio = exp(LnLikelihoodNew - lTrackold);
        r = rand(Parameters.N, 1);
        test = repmat(acceptRatio > r, [1, Parameters.M]);
        Pnew = test.*psamples + (1-test).*Pold;
        lTracknew = test(:,1).*LnLikelihoodNew + (1-test(:,1)).*lTrackold;
        
%         %Display and track Log Likelihood values
%         if(mod(iter, 200) == 0);
%             disp(['Iteration ', num2str(iter), ' of ', num2str(Parameters.NumberIterations)]);
%             disp(num2str(sum(lTracknew)));
%         end
        
        %Keep the best P samples
        P( find(lTracknew>lTrackold),:)=Pnew(find(lTracknew>lTrackold),:);
        Pold=P;
        lTrackold = lTracknew;
      
        %Display Progress bar
        perc = round(iter/Parameters.NumberIterations*100);
        waitbar(iter/Parameters.NumberIterations,h,sprintf('Stage 2/2: MH Sampling...%d%%',perc));
    end
  close(h);
else
    error('not yet implemented');
    P = [];
end
end

function [psamples] = DirichletSample(A)

Y        = randg(A) ; %Sample proportions under consideration
v        = sum(Y,2);
psamples = Y./repmat(v, [1, size(A,2)]);

end
