function [M2,NM2]=compute_Hist_Metric2(img1,img2,index)

img1 = double(img1)/255;
img2 = double(img2)/255;
img1 = img1 - 0.01*(img1==1);
img2 = img2 - 0.01*(img2==1);

N = length(index);
sN = 20000;
gray_level = 64;%round(N^(1/(2+alpha)));
s1 = img1(index);
s2 = img2(index);
rand('state',1);
% if sampling    
%     sindex = round(rand(sN,1)*(N-1))+1;    
%     s1 = s1(sindex);
%     s2 = s2(sindex);
%     gray_level = round(sN^(1/(2+alpha)));
% end

bitLen = floor(log2(gray_level)-0.01) + 1;

sd1 = floor(s1*gray_level);
sd2 = floor(s2*gray_level);

h1 = getHist(sd1,bitLen);
h2 = getHist(sd2,bitLen);

h12 = getHist([sd1 sd2],bitLen);

e1 = -sum(h1.*log(h1));
e2 = -sum(h2.*log(h2));

e12 = -sum(h12.*log(h12));

M2 = e12 - e1- e2;
NM2 = M2/e12;
