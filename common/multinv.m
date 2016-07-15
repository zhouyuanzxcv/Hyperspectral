function X = multinv(M)
% function X = multinv(M)
%
% PURPOSE: Inverse each 2D slice of an array (M) with arbitrary dimensions 
%          support.
%
% INPUT
%   M  : n_D array (m x m x [p x q x ...]), for all possible m=1,2,3,... 
%        and optional higher dimensions.
%   
% OUTPUT
%   X  : n_D array (m x m x [p x q x  ...]), with same size as M.
%
% Inverse every 2D slice (the first two dimensions of M) for multi-dimension
% array M.
%   M(:,:,p,q,...) * X(:,:,p,q,...) = repmat(eye(m),[1,1,p,q,...])
%
% NOTE 1 -- This function may use a large amount of memory for huge array. 
%           Test before usage.
%
% NOTE 2 -- Underdetermined system (more unknowns than equation)
%   The solution is basic solution obtained with sparse mldivide
%   which is not the same as basic solution when calling for full matrix.
%
% See also: multiprod
%
% Author: Xiaodong Qi <i2000s@hotmail.com>
% History: original 26-Apr-2011
% Inspired by Bruno Luong's MultiSolver code to solve a linear equ system.
%
% Copyright (c) 2011, Xiaodong Qi
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

sn = size(M);
m=sn(1);
n=sn(2);
if m~=n
   error('multinv: The first two dimensions of M must be m x m slices.');
end 
p=prod(sn(3:end));
M=reshape(M,[m,n,p]);

% Build sparse matrix and solve
I = reshape(1:m*p,m,1,p);
I = repmat(I,[1 n 1]); % m x n x p
J = reshape(1:n*p,1,n,p);
J = repmat(J,[m 1 1]); % m x n x p
M = sparse(I(:),J(:),M(:));
clear I J
RHS = repmat(eye(m),[p,1]);
X = M \ RHS;
clear RHS M
X = reshape(X, [n p m]);
X = permute(X,[1,3,2]);
X = reshape(X,[n,m,sn(3:end)]);

end % multinv

