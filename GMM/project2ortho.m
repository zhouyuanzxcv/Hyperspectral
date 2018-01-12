function [mus_p,sigmas_p] = project2ortho(mus,sigmas,c,E)
[rows,cols] = size(mus);
mus_p = cell(rows,cols);
sigmas_p = cell(rows,cols);
B = size(E,2);

for i = 1:rows
    for j = 1:cols
        K_j = size(mus{i,j},1);
        mus_p{i,j} = zeros(K_j,B);
        sigmas_p{i,j} = zeros(B,B,K_j);

        for k = 1:size(mus{i,j},1)
            mus_p{i,j}(k,:) = (mus{i,j}(k,:) - c)*E;
            sigmas_p{i,j}(:,:,k) = E'*sigmas{i,j}(:,:,k)*E;
        end
    end
end
