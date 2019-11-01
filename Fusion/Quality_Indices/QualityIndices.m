function Out = QualityIndices(I_HS,I_REF,ratio)
%--------------------------------------------------------------------------
% Quality Indices
%
% USAGE
%   Out = QualityIndices(I_HS,I_REF,ratio)
%
% INPUT
%   I_HS  : target HS data (rows,cols,bands)
%   I_REF : reference HS data (rows,cols,bands)
%   ratio : GSD ratio between HS and MS imagers
%
% OUTPUT
%   Out.cc   : CC
%   Out.sam  : SAM
%   Out.rmse : RMSE
%   Out.ergas: ERGAS
%
%--------------------------------------------------------------------------
[rows,cols,bands] = size(I_REF);

%Remove border from the analysis
I_HS  = I_HS(ratio+1:rows-ratio,ratio+1:cols-ratio,:);
I_REF = I_REF(ratio+1:rows-ratio,ratio+1:cols-ratio,:);

cc = CC(I_HS,I_REF);
Out.cc = mean(cc);
Out.cc_std = std(cc);
[angle_SAM,map] = SAM(I_HS,I_REF);
Out.sam = angle_SAM;
Out.rmse = RMSE(I_HS,I_REF);
Out.ergas = ERGAS(I_HS,I_REF,ratio);

disp(['CC   : ', num2str(Out.cc), ' (std: ',num2str(Out.cc_std),')']);
disp(['SAM  : ' num2str(Out.sam)]);
disp(['RMSE : ' num2str(Out.rmse)]);
disp(['ERGAS: ' num2str(Out.ergas)]);
