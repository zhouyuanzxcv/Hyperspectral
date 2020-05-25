function [Parameters] = BCMParameters(endmembers)

%   This function sets the parameters to be used during the BCM unmixing algorithm
%
% INPUT
%   endmembers - cell - 1xM cell of endmember samples (known) where M is the number of endmembers.
%               DxNumESamples double matrix within each cell, where D is the
%               number of spectral bands and NumESamples is the number of
%               endmember samples
% OUTPUT
%   Parameters - struct - The struct contains the following fields:
%               % These parameters are user defined (can be modified)
%                   1. NumberIterations: Number of iterations
%                   2. K: Number of nearest neighbors
%                   3. c: Number of clusters for spatial K-means
%                   4. SigV: Weighting on SigmaV
%                   5. SigM: Weighting on SigmaMean
%                   6. s: Scaling parameter for location information
%                   7. ECovariance: diagonal covariance on endmembers, currently same for all endmembers
%               % These parameters are computed from input data and known
%               endmembers (does not need modification)
%                   8. MaxNumEMs/M: Number of endmembers
%                   9. bandnum/D: Number of spectral bands/dimensions
%                  10. N: Number of pixels(input data points)
%                  11. Emean: Mean of endmember samples
%                  12. BetaParameters: parameters for fitting a Beta distribution to each band of each endmember
%
% REFERENCE :
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


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Please modify the parameters in this section as needed  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Parameters.NumberIterations = 20000;%Number of Iterations
Parameters.K = 6; % K nearest neighbor number 
Parameters.c = 3;  % Number of clusters for spatial K-means
Parameters.sigV = 100; % Weighting on SigmaV
Parameters.sigM = 0.001; % Weighting on SigmaMean
Parameters.s = 100; % Scaling paramters for location information

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%   The following Parameters are computed from the input, %%%
%%%%%%%   therefore does NOT need modification              %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Parameters.MaxNumEMs = size(endmembers,2);%Number of endmembers
Parameters.bandnum = size(endmembers{1},1);%Number of bands

for i=1:Parameters.MaxNumEMs
Parameters.Emean(i,:) = mean(endmembers{i}');
end


%%Parameters.BetaParameters - double Mat - Mx2xD matrix
%%Parameters for approximating a Beta distribution to each band of each endmember

% %%***Another option is to set diagonal covariance on endmembers***
% Parameters.ECovariance = 0.0013; % diagonal covariance on endmembers
% for i = 1:Parameters.MaxNumEMs 
%     for j = 1:Parameters.bandnum
%     Parameters.BetaParameters(j,1,i) = (((1-Parameters.Emean(i,j))/Parameters.ECovariance) - (1/Parameters.Emean(i,j))) * (Parameters.Emean(i,j).^2);
%     Parameters.BetaParameters(j,2,i) = Parameters.BetaParameters(j,1,i) * ((1/Parameters.Emean(i,j))-1);
%     end
% end


for i = 1:Parameters.MaxNumEMs 
    for j = 1:Parameters.bandnum    
        phat = betafit(endmembers{1,i}(j,:));
        Parameters.BetaParameters(j,1,i)=phat(1);
        Parameters.BetaParameters(j,2,i)=phat(2);
    end
end
%Parameters for Beta on Endmembers

Parameters.aBeta(:,:) = Parameters.BetaParameters(:,1,:);
Parameters.bBeta(:,:) = Parameters.BetaParameters(:,2,:);
Parameters.aBeta = Parameters.aBeta';
Parameters.bBeta = Parameters.bBeta';

    
Parameters.EmeanBeta = zeros(Parameters.MaxNumEMs,Parameters.bandnum);
Parameters.ESigmaBeta = zeros(Parameters.MaxNumEMs,Parameters.bandnum);
for i = 1:Parameters.MaxNumEMs 
    for j = 1:Parameters.bandnum
        Parameters.EmeanBeta(i,j) = Parameters.aBeta(i,j) / (Parameters.aBeta(i,j) + Parameters.bBeta(i,j));
        Parameters.ESigmaBeta(i,j) = (Parameters.aBeta(i,j) *  Parameters.bBeta(i,j)) / (((Parameters.aBeta(i,j) + Parameters.bBeta(i,j)).^2) *(Parameters.aBeta(i,j)+Parameters.bBeta(i,j)+1));
    end
end


 Parameters.ECovariance = mean(mean(Parameters.ESigmaBeta));


end