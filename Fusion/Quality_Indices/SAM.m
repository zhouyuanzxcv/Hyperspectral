function [angle_SAM,map] = SAM(ref,tar)
%--------------------------------------------------------------------------
% Spectral angle mapper (SAM)
%
% USAGE
%   [angle_SAM,map] = SAM(msoriginal,msfused)
%
% INPUT
%   ref : reference HS data (rows,cols,bands)
%   tar : target HS data (rows,cols,bands)
%
% OUTPUT
%   angle_SAM : average SAM in degree (scalar)
%   map       : 2-D map of SAM
%
%--------------------------------------------------------------------------
[rows,cols,bands] = size(tar);
prod_scal = dot(ref,tar,3); 
norm_orig = dot(ref,ref,3);
norm_fusa = dot(tar,tar,3);
prod_norm = sqrt(norm_orig.*norm_fusa);
prod_map = prod_norm;
prod_map(prod_map==0)=eps;
map = acos(prod_scal./prod_map);
prod_scal = reshape(prod_scal,rows*cols,1);
prod_norm = reshape(prod_norm, rows*cols,1);
z=find(prod_norm==0);
prod_scal(z)=[];prod_norm(z)=[];
angolo = sum(sum(acos(prod_scal./prod_norm)))/(size(prod_norm,1));
angle_SAM = real(angolo)*180/pi;

end