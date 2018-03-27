function H1 = solve_for_H1(X, Y1, options)
if nargin < 3
    options = [];
end
Y1Y = Y1'*Y1;
lambda = options.lambda * max(Y1Y(:)); % 1e-3 * N
L = calc_Laplacian_for_H1(size(Y1,2));
H1 = inv(Y1Y + lambda*L)*Y1'*X;


%% The NMF with positivity constraint on H1 does not give correct SRF
% function H1 = solve_for_H1(X, options)
% % requires options.Y1, options.Y1Y, options, Y1Yinv
% H1 = max(options.Y1Yinv * (options.Y1' * X), 0);
% Y1X = options.Y1' * X;
% old_H1 = H1;
% mag = max(H1(:));
% thresh = 1e-6 * mag;
% 
% vals = sum(sum((X - options.Y1*H1).^2));
% for i = 1:500
%     H1 = H1 .* Y1X ./ (options.Y1Y * H1);
%     if mod(i, 10) == 0 && mean(abs(old_H1(:) - H1(:))) < thresh
%         break;
%     end
%     old_H1 = H1;
%     
%     vals = [vals, sum(sum((X - options.Y1*H1).^2))];
% end
% 
% figure,plot(vals);
% 

