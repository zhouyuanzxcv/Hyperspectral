function [P] = BCM(Xim, Parameters, MethodFlag)

% BCM (Beta Compositional Model) Spatial and Spectral Unmixing Algorithm
%       Unmixes Input Data given known Endmember samples 
% 
% REFERENCE :
% X. Du, A. Zare, P. Gader, D. Dranishnikov, 
% “Spatial and Spectral Unmixing Using the Beta Compositional Model,”  
% IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 
% vol. 7, no. 6, pp. 1994-2003, June 2014.
%
% SYNTAX : [P] = BCM(Xim, Parameters, MethodFlag)
%
% INPUTS:
%     Xim - double Mat - N1xN2xD matrix of N1xN2 image data points of dimensionality D 
%           (i.e.  N=N1xN2 pixels with D spectral bands).
%           Note that BCM only takes in input data between 0 and 1!
%   endmembers - double Mat - known endmember samples (provided)
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
%               *Parameters can be modified in [Parameters] = BCMParameters(endmembers) function.
%   MethodFlag - double Mat - 1x1 integer, choose from 1, 2, 3, 4, where each number indicates one BCM approach:
%                       MethodFlag==1 => BCM-Spectral-QP (BCM-Spectral unmixing: Quadratic Programming approach)
%                       MethodFlag==2 => BCM-Spectral-MH (BCM-Spectral unmixing: Metropolis-Hastings Sampling approach)
%                       MethodFlag==3 => BCM-Spatial-QP  (BCM-Spatial unmixing: Quadratic Programming approach)
%                       MethodFlag==4 => BCM-Spatial-MH  (BCM-Spatial unmixing: Metropolis-Hastings Sampling approach)
% OUTPUTS:
%   P - double Mat - NxM matrix of abundances/proportions corresponding to N input
%       pixels and M endmembers
% OTHER m-FILES REQUIRED: 
%   unmix_qpas_correct.m, BCM_Spectral.m, BCM_Spatial.m, BCMParameters.m,     
%   Matlab Statistics Toolbox, Matlab Optimization Toolbox
%
% Author: Xiaoxiao Du, Alina Zare
% University of Missouri, Department of Electrical and Computer Engineering
% Email Address: xdy74@mail.missouri.edu; zarea@missouri.edu
% Created: August 2013
% Latest Revision: October 13, 2014
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


Parameters.methodFlag = MethodFlag;

if Parameters.methodFlag == 1 || Parameters.methodFlag == 2 %BCM-Spectral
    [P, F, E, S] = BCM_Spectral(Xim, Parameters);
    
elseif Parameters.methodFlag == 3 || Parameters.methodFlag == 4 %BCM-Spatial
    [P, F, E, S] = BCM_Spatial(Xim, Parameters);
    
else
    error('not yet implemented');
    P = []; 
end


end