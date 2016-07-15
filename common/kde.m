function [logp, s] = kde(train, test, s)

%KDE Kernel Density Estimation.
%    This function computes a kernel density estimator from a set of examples,
%    by placing Gaussian kernels (with identical masses) on each training data
%    point and adjusting the "widths" to maximize the sum of leave-one-out log
%    densities. In multivariate data, either isotropic or diagonal covariance
%    matrix Gaussians are used.
%
%    usage: [logp, s] = kde(train, test, s)
%
%    inputs:  train (n by D) is a matrix of n training points in D dimensions
%             test  (N by D) is a matrix of test points
%             s     scalar or (1 by D) is a start guess for std devs (optional)
%
%    outputs: logp  (1 by N) is the vector of test log densities
%             s     final std devs
% 
%    The size of the initial guess for s indicates whether isotropic or
%    diagonal covariance is desired. If no initial guess is supplied, diagonal
%    is assumed. The method uses a Newton scheme in the log of s, and usually
%    converges in very few iterations. The Newton steps are checked to ensure
%    that a reasonable fraction (half) of the expected improvement is achieved;
%    otherwise smaller steps are tried. The computational complexity is order
%    (nD)^2 + nND, memory requirement order Dn(N+n). The algorithm may
%    encounter numerical problems if initial guess is too small.
%
%    (C) Copyright Carl Edward Rasmussen, July 4th 2000.

[N, D] = size(test); [n, D] = size(train);  % get number of cases and dimension
if nargin == 2                                % if no start guess is given then
  s = std(train)/(n^(0.5/D));      % use scaled empirical axis aligned std devs
end
P = length(s);                                    % number of parameters to fit

t = repmat(train,[1,1,n]);
if P == 1                                    % if we are fitting a single width
  c = sum((permute(t,[1,3,2])-permute(t,[3,1,2])).^2,3);
else                                                 % else multiple parameters
  c = (permute(t,[1,3,2])-permute(t,[3,1,2])).^2;
end

G = 1; TINY = 1e-10;
ec = exp(-sum(c./repmat(permute(2*s.^2,[1,3,2]),[n,n,1]),3));
secd = repmat(sum(ec-eye(n),2),[1,P]);
f_old = sum(log(secd(:,1)))-n*sum(log(s));
while max(abs(G)) > TINY
  x = shiftdim(sum(repmat(ec,[1,1,P]).*c,1))./secd;
  DE = sum(x,1)./s.^2;
  xx = repmat(x,[1 1 P]);
  for i=1:P        % rewriting this loop as matrix expr would cost too much mem
    DDE(i,1:P) = sum(shiftdim(sum(c.*repmat(ec.*c(:,:,i),[1,1,P])))./secd);
  end
  DDE = (DDE - shiftdim(sum(xx.*permute(xx,[1,3,2]))))./(s'*s).^2;
  [v, l] = eig(DDE-2*diag(DE));                               % eigs of Hessian
  l = max(abs(diag(l)),min(l(:)/TINY));    % control sign and magnitude of eigs
  G = (n*D/P-DE)*v*diag(-1./l)*v';                        % compute Newton step
  G = G/max(1, sqrt(G*G'));                       % don't take too large a step
  eta = 1; s_old = s;
  while eta == 1 | f_new < f_old - eta*G*(n*D/P-DE')/2 - TINY    % improvement?
    s = s_old.*exp(G*eta);
    eta = eta/2;           % if we fail, then try smaller step next time around
    ec = exp(-sum(c./repmat(permute(2*s.^2,[1,3,2]),[n,n,1]),3));
    secd = repmat(sum(ec-eye(n),2),[1,P]);
    f_new = sum(log(secd(:,1)))-n*sum(log(s));
  end
  f_old = f_new;                   % remember function value for next iteration
end

x = repmat(2*s.^2,[n,D/P]);
logp = zeros(N, 1);
for i=1:N          % rewriting this loop as matrix expr would cost too much mem
  cc = sum((repmat(test(i,:),[n,1])-train).^2./x,2);
  hh = min(cc);
  logp(i) = log(mean(exp(hh-cc)))-hh;
end

logp = logp - D*log(2*pi)/2 - D*sum(log(s))/P;            % normalize densities
