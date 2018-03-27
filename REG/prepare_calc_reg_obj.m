function options = prepare_calc_reg_obj(I1, wl, options)
%PREPARE_LSQ_REG Summary of this function goes here
%   Detailed explanation goes here
Y1 = calc_Y1(I1, wl, options);
Y1Y = Y1'*Y1;
options.Y1 = Y1;

% options.U = calc_U(I1, wl);
% options.Y1Y = Y1Y;
% options.Y1Yinv = inv(options.Y1Y);

lambda = options.lambda * max(Y1Y(:));
L = calc_Laplacian_for_H1(size(Y1,2));

options.Y1_invY1LambdaL = Y1*inv(Y1Y + lambda*L);



function Y1 = calc_Y1(I, wl, options)
I_sel = select_relevant_bands(I,wl,options.srf_range);

Y1 = reshape_hsi(I_sel);
Y1 = [ones(size(Y1,1),1),Y1];
