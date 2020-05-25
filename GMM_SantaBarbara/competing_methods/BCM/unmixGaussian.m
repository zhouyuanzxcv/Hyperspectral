function [Pbest] = unmixGaussian(X,Parameters)

% This function performs NCM-MH unmixing (Normal Compositional Model, MH Sampling approach).
%     Update Proportions Part ONLY (given endmembers), Gaussian Sampling method unmix 
% INPUTS 
%   X - double Mat - NxD image data. 
%   Parameters - struct - Parameters.
% OUTPUTS
%   Pbest - double Mat - MxN matrix of abundances/proportions corresponding to N input
%       pixels and M endmembers
%
% Author: Xiaoxiao Du, Alina Zare
% University of Missouri, Department of Electrical and Computer Engineering
% Email Address: xdy74@mail.missouri.edu; zarea@missouri.edu
% Created: August 2013
% Latest Revision: November 3, 2014
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

 NumberIterationsforUpdateP = Parameters.NumberIterations;
 NPts = size(X,1); % number of data points
 M = size(Parameters.Emean,1); %number of endmembers
 D = size(Parameters.Emean,2); %number of bands
 Parameters.AlphaPropVector = ones(1, Parameters.MaxNumEMs);%Parameters for Sampling Proportions, currently all ones (symmetric dirichlet)
 
 P(:,:,1) = DirichletSample(Parameters.AlphaPropVector,size(X,1)); %Initialize P

 onesLengthAlpha = ones(1, Parameters.MaxNumEMs);
 
 E = Parameters.Emean;
 LogLikelihoodOld          = ComputeLogLikelihoodAll(X, E, P, Parameters, NPts);%Initialize LogLikelihoodOld
 LogLikelihoodTraceBest    = LogLikelihoodOld;
 Pbest = P;

%%
for iteration = 2:NumberIterationsforUpdateP+1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Update Proportions                       %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Y                     = randg(Parameters.AlphaPropVector(1), NPts, Parameters.MaxNumEMs) ; %Sample proportions under consideration
        v                     = sum(Y,2);
        samples               = Y./v(:, onesLengthAlpha);
        LogLikelihoodNew      = ComputeLogLikelihoodAll(X, E(:,:,1), samples, Parameters, NPts);
        Ratios                = exp(LogLikelihoodNew - LogLikelihoodOld(:,1));
        rands                 = rand(NPts,1);
        Vals                  = rands < Ratios;
        Vrep                  = Vals(:, onesLengthAlpha);
        P(:,:,1)              = samples.*Vrep  + P( :,:,1).*(1-Vrep);
        LogLikelihoodOld(:,1) = (1-Vals).*LogLikelihoodOld(:,1) + (Vals).*LogLikelihoodNew;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
         LogLikelihoodNewTrace     = LogLikelihoodOld ;

        for i =1:NPts
            if  LogLikelihoodNewTrace(i,:) > LogLikelihoodTraceBest(i,:)
            Pbest(i,:) = P(i,:);
            LogLikelihoodTraceBest(i,:) = LogLikelihoodNewTrace(i,:);
            end
        end
        
%         %Display and track Log Likelihood values
%         if(mod(iteration, 200) == 0);
%             disp(['Iteration ', num2str(iteration), ' of ', num2str(NumberIterationsforUpdateP)]);
%             disp(['sum of LogLikelihoodOld = ' num2str(sum(LogLikelihoodOld))]);
%         end

end

end

function [LogLikelihoodAll] = ComputeLogLikelihoodAll(X, E, P, Parameters, NumPoints)
%compute log likelihood
NumSets = size(E,3);
LogLikelihoodAll = zeros(NumPoints, NumSets);
C = (-1/2)*1./sum(P.*P*Parameters.ECovariance,2);
N = -1*log(squeeze(sum(P.*P*Parameters.ECovariance,2)).*size(X,2));
for i = 1:NumSets
    Recon = P(:,:,i)*E(:,:,i);
    D = (X - Recon)';
    D = sum(D.*D);
    LogLikelihoodAll(:,i) = N(:,i) + C(:,1,i).*D';
end

end

function [samples] = DirichletSample(alpha, numSamples)
%Generate samples from Dirichlet distribution
    if(size(alpha,1) == 1)
        Y = randg(repmat(alpha, [numSamples,1]));
        samples = Y./repmat(sum(Y')', [1, length(alpha)]);
    else
        Y = randg(alpha);
        samples = Y./repmat(sum(Y')', [1, size(alpha,2)]);
    end
end

