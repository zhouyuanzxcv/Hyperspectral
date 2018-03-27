function [T,degree,t,s] = tform2myT(tform)
T = tform.T';
T = inv(T);

t = T([1,2],3);

% assume s1 and s2 appear on the first and the second row of T
s = sqrt(T(1:2,1).^2 + T(1:2,2).^2);
R = T(1:2,1:2) ./ repmat(s,[1,2]);
theta = atan2(R(2,1), R(1,1));
degree = theta * 180/pi;

