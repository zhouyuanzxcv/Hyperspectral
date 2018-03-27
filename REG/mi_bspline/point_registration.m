function [O_trans,Spacing,Xreg]=point_registration(sizeI,Xmoving,Xstatic,Options)
% This function creates a 2D or 3D b-spline grid, which transform space to
% fit a set of points Xmoving to set of corresponding points in  Xstatic. 
% Usefull for:
% - For image-registration based on corresponding landmarks like in
% Sift or OpenSurf (see Mathworks). 
% - 2D and 3D Spline based Data gridding and surface fitting
% - Smooth filtering of 2d / 3D Point data.
%
%   [O_trans,Spacing,Xreg]=point_registration(sizeI,Xstatic,Xmoving,Options);
%
% Inputs,
%   sizeI : The size of the (virtual) image/space which will be warped
%            With the b-spline grid
%   Xmoving : List with 2D or 3D points N x 2, or N x 3, these points will be
%            warped be the fitted b-sline grid to transform to the 
%            static points in Xstatic
%   Xstatic : List with 2D or 3D points N x 2, or N x 3, corresponding with Xmoving
%   Options : Struct with options, see below.
%
% Outputs,
%   Grid: The b-spline controlpoints, can be used to transform another
%       image in the same way: I=bspline_transform(Grid,I,Spacing);
%   Spacing: The uniform b-spline knot spacing
%   Xreg: The points in Xmoving transformed by the fitted b-spline grid.
%
% Options,
%   Options.Verbose: Display Debug information 0,1 or 2
%   Options.MaxRef : Maximum number of grid refinements steps (default 5)
%
%
% % Example. Image Warp 2D based on landmarks,
%  % Load corresponding landmarks
%    load('images/starpoints.mat');
%  % Load the images
%    I1=im2double(imread('images/star1.png')); 
%    I2=im2double(imread('images/star2.png')); 
%  % Fit the bspline grid to the corresponding landmarks
%    options.Verbose=true;
%    [O_trans,Spacing]=point_registration(size(I1),[x2(:) y2(:)],[x1(:) y1(:)],options);
%  % Transform the 2D image  
%    Ireg=bspline_transform(O_trans,I1,Spacing,3);
%  % Show the result
%    figure,
%     subplot(1,3,1),imshow(I1); title('Moving Image'); 
%     hold on;for i=1:10, plot([y1(i) y2(i)],[x1(i) x2(i)],'b'); end
%     subplot(1,3,2),imshow(I2); title('Static Image');
%     subplot(1,3,3),imshow(Ireg); title('Registered Image')
%     hold on;for i=1:10, plot([y1(i) y2(i)],[x1(i) x2(i)],'b'); end 
%  % Show b-spline grid
%    Igrid=make_grid_image(Spacing,size(I1));
%    figure, 
%     subplot(1,2,1), imshow(Igrid)
%     Ireg=bspline_transform(O_trans,Igrid,Spacing,3);
%     subplot(1,2,2), imshow(Ireg)
%     hold on;for i=1:10, plot([y1(i) y2(i)],[x1(i) x2(i)],'b'); end
%
% % Example, B-spline Fitting to Grid Sparse Data
%   % Load image
%     I=im2double(rgb2gray(imread('images/lena.jpg')));
%   % Select 3000 points from image
%     x=round(rand(3000,1)*255)+1; y=round(rand(3000,1)*255)+1;
%     z=I(sub2ind(size(I),x,y))+16; % Z is equal to intensity
%   % Create an image with the 3000 pixels
%     T=zeros(256,256); T(sub2ind([256 256],x,y))=1; T=I.*T;
%   % Fit a B-spline grid
%     z2=16*ones(size(x));
%     options.MaxRef=5;
%     options.Verbose=true;
%     [O_trans,Spacing,X3]=point_registration([256 256 32],[x y z2],[x y z],options);
%   % Create a uniform pixel grid
%     [xc,yc]=ndgrid(1:0.5:256,1:0.5:256); xc=xc(:); yc=yc(:); zc=16*ones(size(xc));
%   % Transform the points with the b-spline grid
%     J = bspline_trans_points_double(O_trans,Spacing,[xc yc zc]);
%   % Transform z-coordinate back to intensity image
%     J = reshape(J,[511 511 3]); J=J(:,:,3)-16;
%   % Compare with matlab griddata
%     K = reshape(griddata(x,y,z,xc,yc,'cubic')-16,[511 511]);
%   % Show the results
%    figure,
%     subplot(1,3,1),imshow(T); title('original');
%     subplot(1,3,2),imshow(J); title('b-spline fit');
%     subplot(1,3,3),imshow(K); title('Matlab Griddata');
%
% % Example, 3D
%    Xmoving = round(rand(10,3)*50+1);
%    Xstatic = Xmoving + round(rand(10,3)*20+1);
%    options.Verbose=true;
%    Options.MaxRef=7;
%    % Fit a bspline grid to transform points Xstatic into Xmoving
%    [O_trans,Spacing]=point_registration([50 50 50],Xmoving,Xstatic,options);
%    % Transforming some other point with the b-spline grid
%    X3 = round(rand(10,3)*50+1);
%    X3t = bspline_trans_points_double(O_trans,Spacing,X3);
% 
% Example, 2D Diffeomorphic Warp
%    Xstatic=[1 1;
%      1 128;
%      64+32 64
%      64-32 64
%      128 1;
%      128 128];
% 
%   Xmoving=[1 1;
%      1 128;
%      64-32 64
%      64+32 64
%      128 1;
%      128 128];
%   option=struct; options.MaxRef=4;
%   sizeI=[128 128]; 
%   [O_trans,Spacing]=point_registration(sizeI,Xstatic,Xmoving,options);
%   [O_trans,Spacing]=MakeDiffeomorphic(O_trans,Spacing,sizeI);
%
%   Igrid=make_grid_image(Spacing*2,sizeI);
%   [Ireg,B]=bspline_transform(O_trans,Igrid,Spacing,3);
%   figure, imshow(Ireg)
%
%  Function is written by D.Kroon University of Twente (August 2010)
  
% add all needed function paths
add_function_paths;

% Disable warning
warning('off', 'MATLAB:maxNumCompThreads:Deprecated')

% Process inputs
defaultoptions=struct('Verbose',false,'MaxRef',5);
if(~exist('Options','var')), Options=defaultoptions;
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags), if(~isfield(Options,tags{i})), Options.(tags{i})=defaultoptions.(tags{i}); end, end
    if(length(tags)~=length(fieldnames(Options))),
        warning('register_images:unknownoption','unknown options found');
    end
end

switch(size(Xstatic,2))
    case 2
        % Calculate max refinements steps
        MaxItt=min(floor(log2(sizeI(1:2)/2)));

        % set b-spline grid spacing in x and y direction
        Spacing=[2^MaxItt 2^MaxItt];
    case 3
        % Calculate max refinements steps
        MaxItt=min(floor(log2(sizeI(1:3)/2)));

        % set b-spline grid spacing in x,y and z direction
        Spacing=[2^MaxItt 2^MaxItt 2^MaxItt];
end

% Make an initial uniform b-spline grid
O_ref = make_init_grid(Spacing,sizeI);

% Calculate difference between the points
R=Xstatic-Xmoving;

% Initialize the grid-update needed to b-spline register the points
O_add=zeros(size(O_ref));

% Loop through all refinement itterations
for i=1:Options.MaxRef
    if(Options.Verbose)
        disp('.');
        disp(['Iteration : ' num2str(i) '/' num2str(Options.MaxRef)]);
        disp(['Grid size : ',num2str(size(O_ref))]);
    end
    
    % Make a b-spline grid which minimizes the difference between the
    % corresponding points
    O_add=bspline_grid_fitting(O_add,Spacing,R,Xmoving);
    
    % Warp the points
    Xreg=bspline_trans_points_double((O_ref+O_add),Spacing,Xmoving);
    
    % Calculate the remaining difference between the points
    R=Xstatic-Xreg;
    if(Options.Verbose)
        err=sqrt(sum(R.^2,2));
        disp(['Mean Distance : ',num2str(mean(err))]);
    end
    
    if(i<Options.MaxRef)
        % Refine the update-grid and reference grid
        switch(size(Xstatic,2))
            case 2
                O_add = refine_grid(O_add,Spacing,sizeI(1:2));
                [O_ref ,Spacing]=refine_grid(O_ref,Spacing,sizeI(1:2));
            case 3
                O_add = refine_grid(O_add,Spacing,sizeI);
                [O_ref ,Spacing]=refine_grid(O_ref,Spacing,sizeI);
        end
    end
end
% The final transformation grid, is the reference grid with update-grid
% added
O_trans=O_ref+O_add;

function add_function_paths()
try
    functionname='point_registration.m';
    functiondir=which(functionname);
    functiondir=functiondir(1:end-length(functionname));
    addpath([functiondir '/functions'])
    addpath([functiondir '/functions_affine'])
    addpath([functiondir '/functions_nonrigid'])
catch me
    disp(me.message);
end

 