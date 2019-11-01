function ERGAS = ERGAS(I,I_Fus,Resize_fact)
%--------------------------------------------------------------------------
% Erreur relative globale adimensionnelle de synthese (ERGAS)
%
% USAGE
%   ERGAS = ERGAS(I,I_Fus,Resize_fact)
%
% INPUT
%   I          : reference HS data (rows,cols,bands)
%   I_Fus      : target HS data (rows,cols,bands)
%   Resize_fact: resize factor
%
% OUTPUT
%   ERGAS : ERGAS (scalar)
%
%--------------------------------------------------------------------------
I = double(I);
I_Fus = double(I_Fus);

Err=I-I_Fus;
ERGAS=0;
for iLR=1:size(Err,3),
    ERGAS=ERGAS+mean2(Err(:,:,iLR).^2)/(mean2((I(:,:,iLR))))^2;   
end

ERGAS = (100/Resize_fact) * sqrt((1/size(Err,3)) * ERGAS);

end