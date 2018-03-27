function [degree,t,s,sigma] = reg_init(rgb1, I1, wl, options)
% figure,imshow(uint8(rgb1));
% figure,imshow(uint8(retrieve_rgb(I1,wl)*255));
rgb1 = double(rgb1);
I1 = I1 * 255;

s = options.s;
sigma = 2*s;
s = [s,s]'; % s is the initial scaling

I2 = transform(rgb1, create_T(0,[0,0]), s, sigma, ...
    floor(size(rgb1,2)/s(1)), floor(size(rgb1,1)/s(2)), options);

init_method = parse_param(options,'init_method','pc');
switch init_method
    case 'pc' % use phase correlation and least squares refinement
        disp('Use phase correlation to do initial registration.');
        % s1 is the scale based on the reduced color image, where T includes s1
        [T,degree_pc,t_pc,s1_pc] = reg_phasecorr(I2, I1, wl, options);
        
        if options.show_figure
            I3 = transform(I2, T, [1,1], 0.01, size(I1,2), size(I1,1));
            figure,imshow(uint8(I3));
        end
        
        disp('Refine the initial result by the registration metric.');
        s1_pc = [s1_pc,s1_pc];
        %     options.reg_metric = 'MI histogram';
        [T,degree,t,s1] = reg_refine(degree_pc,t_pc,s1_pc,I2,I1,options);
        
        %     I3 = transform(I2,T,[1,1],0.01,size(I1,2),size(I1,1));
        %     figure,imshow(uint8(I3));
        
        %     options.reg_metric = [];
        
    case 'lsq' % only use refinement with BCD
        [T,degree,t,s1] = reg_refine(0,[0,0],[1,1],I2,I1,options);
    case 'mi' % use matlab routine (multiresolution mutual information)
        [optimizer, metric] = imregconfig('multimodal');
        metric.NumberOfHistogramBins = 20;
        optimizer.InitialRadius = 1e-3;
        transformType = 'similarity';
        rgb_ind = find_rgb_ind(wl);
        input = I2(:,:,1);
        template = I1(:,:,rgb_ind(1));
        tform = imregtform(input,template,transformType,optimizer,metric);
        [T,degree,t,s1] = tform2myT(tform);
    otherwise
end

% save('result_init.mat','T','degree','t','s1','s','sigma');

if options.show_figure % show transformed initial condition in the low scale
    I3 = transform(I2,T,[1,1],0.01,size(I1,2),size(I1,1));
    figure('name','initial condition in low scale');
    imshow(uint8(I3));
end

[t,s] = localTrans2GlobalTrans(T,s);

% I4 = transform(rgb1,create_T(degree,t),[s,s],lambda,size(I1,2),size(I1,1));
% figure,imshow(uint8(I4));

function [T,degree,t,s] = reg_refine(degree,t,s,moving,fixed,options)
t_start = tic;

pa = ParameterAnalysis();
pa.StepSizeMode = 2;
pa.CheckValueInListSearchOptions = 0;
list_params = {'degree','tx','ty','sx','sy'};
range = [-1e6,1e6;-1e6,1e6;-1e6,1e6;-1e6,1e6;-1e6,1e6]';
fcn_run = @(options1) calc_val(options1,moving,fixed,options);
step_size = [1, 1, 1, 0.02, 0.02];

% pa.NeighborhoodMode = 7*ones(1,5);
% pa.Algorithm = 'BrutalForce';
% pa.BrutalForceDepth = 3;

pa.NeighborhoodMode = 5*ones(1,5);
pa.Algorithm = 'BlockCoordinateDescent';
pa.BlockCoordinateDescentDepth = 5;
pa.BlockCoordinateDescentGroup = {1,[2,3],[4,5]};

options1 = [];
options1.degree = degree;
options1.tx = t(1);
options1.ty = t(2);
options1.sx = s(1);
options1.sy = s(2);

options1 = pa.autoParamSelection(fcn_run, options1, list_params, ...
    step_size, range);
degree = options1.degree;
t = [options1.tx, options1.ty]';
s = [options1.sx, options1.sy]';
T = create_T(degree,t,'',s);

disp(['Elapsed time for initial refinement is ', num2str(toc(t_start))]);

function val = calc_val(options1,I,I1,options)
degree = options1.degree;
t = [options1.tx, options1.ty]';
s = [options1.sx, options1.sy]';
T = create_T(degree,t,'',s);
sigma = 0.01; % make sure the PSF occupies only 1 pixel
val = eval_obj_fun(I, I1, T, [1,1], sigma, options);


function [t,s] = localTrans2GlobalTrans(T, s)
T1 = create_T(0,[0,0],'',s) * T;
t = T1([1,2],3)';
s1 = sqrt(sum(T1([1,2],1).^2));
s2 = sqrt(sum(T1([1,2],2).^2));
s = [s1,s2];


function [T,degree,t,s] = reg_phasecorr(I1, I, wl, options)
% figure,imshow(uint8(I1));
% figure,imshow(uint8(retrieve_rgb(I,wl)*255));
% lambda = 0.01;
t_start = tic;

if 0 % use the original phase correlation
    if size(I,3) > 3
        rgb_ind = find_rgb_ind(wl);
    else
        rgb_ind = 1;
    end

    tform = imregcorr(I1(:,:,1), I(:,:,rgb_ind(1)), 'similarity', 'window', 1);
    [T,degree,t,s] = tform2myT(tform);

%     I2 = imwarp(I1, tform, 'OutputView', imref2d(size(I)));
%     figure,imshow(uint8(I2));
elseif 1 % use the new phase correlation
    windowing = 1;
    [tform, degree, s] = reg_hyper_fft(I1, I, wl, 'similarity', windowing);
    [T,degree1,t,s1] = tform2myT(tform);
%     I2 = imwarp(I1, tform, 'OutputView', imref2d(size(I)));
%     figure,imshow(uint8(I2));
    
%     I3 = transform(I1, tform2myT(tform), [1,1], 0.01, size(I,2), size(I,1));
%     figure,imshow(uint8(I3));
else % use the translation version with brutal force of angle 
    S = 1;
    degrees = (-20:20);
    ts = zeros(length(degrees),2);
    vals = zeros(1,length(degrees));
    textprogressbar('Use phase correlation to calculate initial parameters: ');
    for i = 1:length(degrees)
        degree = degrees(i);
        theta1 = degree*pi/180;
        tform1 = affine2d([S.*cos(theta1) -S.*sin(theta1) 0; S.*sin(theta1) S.*cos(theta1) 0; 0 0 1]);
    
        [I1_rotated,~] = imwarp(I1,tform1,'SmoothEdges', true);

        [tform, degree, s] = reg_hyper_fft(I1_rotated, I, wl, 'translation', 0);
        [T,degree1,t,s1] = tform2myT(tform);
        ts(i,:) = t';
        vals(i) = eval_obj_fun(I1_rotated, I, T, [1,1], 0.01, options);
        textprogressbar(i/length(degrees)*100);
    end
    vals(vals==0) = Inf;
    textprogressbar('done');
    
    [~,ind] = min(vals);
    t = ts(ind,:);
    degree = degrees(ind);
end

disp(['Elapsed time for phase correlation is ', num2str(toc(t_start))]);


