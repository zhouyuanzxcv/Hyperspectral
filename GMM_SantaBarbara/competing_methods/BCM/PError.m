function [PErrorOutput, PErrorOutputmean, PErrorOutputvar] = PError(N,Ptrue,Pestimate)

% This function computes the proportion erros per pixel per endmember 
%       given true and estimated proportion values.
%
% SYNTAX : [PErrorOutput, PErrorOutputmean, PErrorOutputvar] = PError(N,Ptrue,Pestimate)
%
% INPUTS 
%     N:             Total number of data points. N=N1*N2 in hyperspectral
%                    image data.
%     Ptrue:         True proportion values. N*M, where M is the number of
%                    endmembers.
%     Pestimate:     Estimated proportion values. N*M.
% OUTPUTS
%     PErrorOutput:       PError per pixel
%     PErrorOutputmean:   Mean of PErrorOutput, PError per pixel per
%                         endmember. (This is the result in the Tables in the paper)
%     PErrorOutputvar:    Variance of PErrorOutput across endmembers.
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

  M = size(Ptrue,2);
  t =   (Ptrue- Pestimate).^2;
  PErrorOutput = sum(sum(t)) / N; % PError per pixel
  

PErrorOutputmean = PErrorOutput/M; %Mean of PErrorOutput, PError per pixel per endmember
PErrorOutputvar= var(PErrorOutput); %Variance of PErrorOutput

%%Plot histogram across iterations for one N. FOR DEMO PURPOSES.
% figure(100);
% plot(PErrorOutput);
% xlabel('N');ylabel('PError');
% title('Plot of PErrorOutput');
% figure(200);
% hist(PErrorOutput);
% xlabel('PError');ylabel('# of Data Points');
% title('Histogram of PErrorOutput');

end