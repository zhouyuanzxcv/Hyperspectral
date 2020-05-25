function [x,options] = lsqnonneg_fast(C,d,options)
%LSQNONNEG Linear least squares with nonnegativity constraints.
%
%function [x,resnorm,resid,exitflag,output,lambda] =
%lsqnonneg(C,d,x0,options)
%
%   X = LSQNONNEG(C,d) returns the vector X that minimizes NORM(C*X - d)
%   subject to X >= 0. C and d must be real.
%
%   X = LSQNONNEG(C,d,X0) uses X0 as the starting point if all(X0 > 0);
%   otherwise the default is used. The default start point is the 
%   origin (the default is used when X0==[] or when only two input 
%   arguments are provided). 
%
%   X = LSQNONNEG(C,d,X0,OPTIONS) minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, an argument
%   created with the OPTIMSET function.  See OPTIMSET for details.  Used
%   options are Display and TolX. (A default tolerance TolX of 
%   10*MAX(SIZE(C))*NORM(C,1)*EPS is used). 
%   
%   [X,RESNORM] = LSQNONNEG(...) also returns the value of the squared 2-norm of 
%   the residual: norm(C*X-d)^2.
%
%   [X,RESNORM,RESIDUAL] = LSQNONNEG(...) also returns the value of the  
%   residual: C*X-d.
%   
%   [X,RESNORM,RESIDUAL,EXITFLAG] = LSQNONNEG(...) returns an EXITFLAG that 
%   describes the exit condition of LSQNONNEG.  
%   If EXITFLAG is:
%     1 then LSQNONNEG converged with a solution X.
%     0 then the iteration count was exceeded. Increasing the tolerance
%       (OPTIONS.TolX) may lead to a solution.
%  
%   [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT] = LSQNONNEG(...) returns a structure
%   OUTPUT with the number of steps taken in OUTPUT.iterations and the type 
%   of algorithm used in OUTPUT.algorithm.
%
%   [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA] = LSQNONNEG(...) returns 
%   the dual vector LAMBDA  where LAMBDA(i) <= 0 when X(i) is (approximately) 0 
%   and LAMBDA(i) is (approximately) 0 when X(i) > 0.
% 
%   See also LSCOV, SLASH.

%   L. Shure 5-8-87
%   Revised, 12-15-88,8-31-89 LS, 5-26-98 MAB.
%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 1.6 $  $Date: 1998/08/25 15:18:02 $

% Reference:
%  Lawson and Hanson, "Solving Least Squares Problems", Prentice-Hall, 1974.

if nargin<3
    %disp('Calculating optimset');
    
    defaultopt = optimset('display','final','TolX','10*eps*norm(C,1)*length(C)');
    if ~isreal(C) || ~isreal(d), error('C and d must be real.'); end
    options = [];
    options = optimset(defaultopt,options);
    printtype = optimget(options,'display');
    tol = optimget(options,'tolx');
    
    % In case the defaults were gathered from calling: optimset('fminsearch'):
    c=C;
    if ischar(tol)
        tol = eval(tol);
    end
    
    switch printtype
        case {'none','off'}
            verbosity = 0;
        case 'iter'
            warning('''iter'' value not valid for ''Display'' parameter for lsqnonneg');
            verbosity = 2;
        case 'final'
            verbosity = 1;
        otherwise
            error('Bad value for options parameter: ''Display''');
    end
end
tol = optimget(options,'tolx');
c=C;
if ischar(tol)
    tol = eval(tol);
end

[m,n] = size(C);
P = zeros(1,n);
Z = 1:n;
x = P';

ZZ=Z;
resid = d-C*x;
w = C'*(resid);

% set up iteration criterion
outeriter = 0;
iter = 0;
itmax = 3*n;
exitflag = 1;

% outer loop to put variables into set to hold positive coefficients
while any(Z) && any(w(ZZ) > tol)
   outeriter = outeriter + 1;
   [wt,t] = max(w(ZZ));
   t = ZZ(t);
   P(1,t) = t;
   Z(t) = 0;
   PP = find(P);
   ZZ = find(Z);
   nzz = size(ZZ);
   CP(1:m,PP) = C(:,PP);
   CP(:,ZZ) = zeros(m,nzz(2));
   z = pinv(CP)*d;
   z(ZZ) = zeros(nzz(2),nzz(1));
   % inner loop to remove elements from the positive set which no longer belong
   while any((z(PP) <= tol))
      iter = iter + 1;
      if iter > itmax
         if verbosity 
            warnstr = sprintf('Exiting: Iteration count is exceeded, exiting LSQNONNEG.', ...
               '\n','Try raising the tolerance (OPTIONS.TolX).');
            disp(warnstr);
         end
         exitflag = 0;
         output.iterations = outeriter;
         resnorm = sum(resid.*resid);
         x = z;
         lambda = w;
         return
      end
      QQ = find((z <= tol) & P');
      alpha = min(x(QQ)./(x(QQ) - z(QQ)));
      x = x + alpha*(z - x);
      ij = find(abs(x) < tol & P' ~= 0);
      Z(ij)=ij';
      P(ij)=zeros(1,length(ij));
      PP = find(P);
      ZZ = find(Z);
      nzz = size(ZZ);
      CP(1:m,PP) = C(:,PP);
      CP(:,ZZ) = zeros(m,nzz(2));
      z = pinv(CP)*d;
      z(ZZ) = zeros(nzz(2),nzz(1));
   end
   x = z;
   resid = d-C*x;
   w = C'*(resid);
end

lambda = w;
resnorm = sum(resid.*resid);
output.iterations = outeriter;
output.algorithm = 'active-set using svd';

verbosity =0;

if verbosity > 0
   disp('Optimization terminated successfully.');   
end

