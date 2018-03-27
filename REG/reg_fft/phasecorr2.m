function d = phasecorr2(I1, I)
% new phase correlation for the hyperspectral image
size_I1 = size(I1);
size_I = size(I);
outSize = size_I1(1:2) + size_I(1:2) - 1;
I = cat(3, ones(size(I,1), size(I,2)), I);
ds = [];
N = prod(outSize(1:2));

lambda = 0;
epsilon = 0;

d1_num = size(I,3);
L = calc_Laplacian_for_H1(d1_num);
vecL = L(:);
B = vecL * vecL';

A1 = zeros(outSize(1), outSize(2), d1_num^2);
a = zeros(outSize(1), outSize(2), d1_num);
b = zeros(outSize(1), outSize(2), d1_num);
Fis = zeros(outSize(1), outSize(2), d1_num);

for k = 1:size(I1,3)
    F1 = fft2(I1(:,:,k), outSize(1), outSize(2));
        
    for i = 1:size(I,3)
        Fi = fft2(I(:,:,i), outSize(1), outSize(2));
        Fi1 = Fi .* conj(F1) ./ (abs(F1+eps).^2);
        
        a(:,:,i) = real(Fi1);
        b(:,:,i) = imag(Fi1);
        Fis(:,:,i) = Fi1;
    end
    
    % calculate the kronecker product manually
    for i = 1:d1_num
        for j = 1:d1_num
            A1(:,:,(i-1)*d1_num + j) = a(:,:,i).*a(:,:,j) + b(:,:,i).*b(:,:,j);
        end
    end
    
    A = reshape(A1, [N, d1_num^2]);
    A1A = A'*A;
    A1A_mag = max(abs(A1A(:)));
    lambda1 = lambda * A1A_mag;
    epsilon1 = epsilon * A1A_mag;
    q = (A1A + lambda1*B + epsilon1*eye(d1_num^2, d1_num^2)) \ (A'*ones(N,1));
    h = sum(reshape(q, [d1_num,d1_num]),1);
    
    Fi1 = zeros(outSize(1), outSize(2));
    for i = 1:length(h)
        Fi1 = Fi1 + h(i) * Fis(:,:,i);
    end
    
    Fi1 = Fi1./abs(Fi1);

    d = ifft2(Fi1, 'symmetric');
    
    if max(-d(:)) > max(d(:))
        d = -d;
    end

    ds = cat(3, ds, d);
end

% peaks = zeros(1,size(ds,3));
% for k = 1:length(peaks)
%     peaks(k) = max(max(ds(:,:,k)));
% end
% [~,ind] = max(peaks);
% d = ds(:,:,ind);
d = mean(ds, 3);



