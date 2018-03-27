function [ output_args ] = test_image_registration( input_args )
%TEST_IMAGE_REGISTRATION Summary of this function goes here
%   Detailed explanation goes here
Imoving=imread('images/lenag1.png');
Istatic=imread('images/lenag2.png');

% Register the images
[Ireg,O_trans,Spacing,M,B,F] = image_registration(Imoving,Istatic);

% Show the registration result
figure,
subplot(2,2,1), imshow(Imoving); title('moving image');
subplot(2,2,2), imshow(Istatic); title('static image');
subplot(2,2,3), imshow(Ireg); title('registerd moving image');
% Show also the static image transformed to the moving image
Ireg2=movepixels(Istatic,F);
subplot(2,2,4), imshow(Ireg2); title('registerd static image');

% Show the transformation fields
figure,
subplot(2,2,1), imshow(B(:,:,1),[]); title('Backward Transf. in x direction');
subplot(2,2,2), imshow(F(:,:,2),[]); title('Forward Transf. in x direction');
subplot(2,2,3), imshow(B(:,:,1),[]); title('Backward Transf. in y direction');
subplot(2,2,4), imshow(F(:,:,2),[]); title('Forward Transf. in y direction');

% Calculate strain tensors
E = strain(F(:,:,1),F(:,:,2));
% Show the strain tensors
figure,
subplot(2,2,1), imshow(E(:,:,1,1),[-1 1]); title('Strain Tensors Exx');
subplot(2,2,2), imshow(E(:,:,1,2),[-1 1]); title('Strain Tensors Exy');
subplot(2,2,3), imshow(E(:,:,2,1),[-1 1]); title('Strain Tensors Eyx');
subplot(2,2,4), imshow(E(:,:,2,2),[-1 1]); title('Strain Tensors Eyy');


end

