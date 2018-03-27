function T1 = create_T(degree, b, flag, s)
% Create transform matrix. If flag = 'inv', we are to rotate and translate
% the image by degree being in [0, 180], b = [bx, by]
if nargin < 3
    flag = '';
end
if nargin < 4
    s = 1;
end

if length(s) == 1
    sx = s;
    sy = s;
else
    sx = s(1);
    sy = s(2);
end

theta = degree * pi / 180;
T1 = [sx*cos(theta) -sy*sin(theta) b(1);sx*sin(theta) sy*cos(theta) b(2);0 0 1];
if strcmp(flag, 'inv')
    T1 = inv(T1);
end

