function val = calc_reg_metric(I1, I2, options)
%CALC_REG_METRIC Summary of this function goes here
%   Detailed explanation goes here
% I1 is the fixed image, I2 is the moving image

% val = eval_obj_fun_nn(I1, I2, options);
if isfield(options,'reg_metric') && ~isempty(options.reg_metric)
    if strcmp(options.reg_metric, 'MI histogram') || ...
        strcmp(options.reg_metric, 'MI kdensity')
        val = eval_obj_fun_mi(I1, I2, options);
    end
else
    val = eval_obj_fun_lsq(I1, I2, options);
end

function val = eval_obj_fun_lsq(I1, I2, options)
X = reshape_hsi(I2);
%% old lsq formulation
% val = sum(sum((X - options.U*(options.U'*X)).^2));

%% NMF solution for H1
% H1 = solve_for_H1(X, options);
% val = sum(sum((X - options.Y1 * H1).^2));

%% new lsq formulation
val = sum(sum((X - options.Y1_invY1LambdaL * (options.Y1' * X)).^2));

function val = eval_obj_fun_mi(I1, I2, options)
rgb_ind = find_rgb_ind(options.wl);
img1 = I2(:,:,1);
img2 = I1(:,:,rgb_ind(1));
if max(img1(:)) <= 1
    img1 = img1*255;
end
if max(img2(:)) <= 1
    img2 = img2*255;
end

index = (1:size(I1,1)*size(I1,2))';
if strcmp(options.reg_metric, 'MI histogram')
    [M2,NM2]=compute_Hist_Metric2(img1,img2,index);
elseif strcmp(options.reg_metric, 'MI kdensity')
    M2 = calc_mi_kdensity(img1,img2,index);
end
val = M2;

function val = eval_obj_fun_nn(I1, I2, options)
L1 = options.L1;
L2 = options.L2;
Y = reshape_hsi(I2);
val1 = zeros(size(Y,1),1);
val2 = zeros(size(Y,1),1);
for k = 1:length(L1)
    val1 = val1 + sum((L1{k} * Y).^2, 2);
    val2 = val2 + sum((L2{k} * Y).^2, 2);
end
val = sum(val1) / sum(val2);

function val = eval_obj_fun_nn1(I1, I2, options)
J1 = options.J1;
J2 = options.J2;
K = length(J1{1});
Y = reshape_hsi(I2);
val1 = 0;
val2 = 0;
for i = 1:length(J1)
    val1 = val1 + sum(sum((repmat(Y(i,:),[K 1]) - Y(J1{i},:)) .^ 2));
    val2 = val2 + sum(sum((repmat(Y(i,:),[K 1]) - Y(J2{i},:)) .^ 2));
end
val = val1 / val2;

function val = eval_obj_fun_lle(I1, I2, options)
Y = reshape_hsi(I2);
val = sum(sum((options.L * Y).^2));

function val = eval_obj_fun_mani_area(I1, I2, options)
I3 = cat(3,I1,I2);
val = calc_area(I3);
val1 = 1;
for k = 1:size(I2,3)
    val1 = val1 * (calc_area(I2(:,:,k)))^(1/size(I3,3));
end
val = val / val1;

function val = eval_obj_fun_mani_area1(I1, I2, options)
I3 = cat(3,I1,I2);
val = calc_area(I3);
val1 = calc_area(I2);
val = val / val1;

function val = eval_obj_fun_mani_area2(I1, I2, options)
I3 = cat(3,I1,I2);
I4 = I3 / 255;
val = calc_area(I4);

function val = eval_obj_fun_intensity_difference(I1, I2, options)
I1 = I1(:,:,[54,28,11]);
val = mean(abs(I2(:) - I1(:)));



