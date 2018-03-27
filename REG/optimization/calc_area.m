function A = calc_area(I)
% The derivative uses central difference and only applies to the interior
% points
der_x = calc_der_x(I(2:end-1,:,:));
der_y = calc_der_y(I(:,2:end-1,:));
G11 = sum(der_x.^2, 3) + 1;
G12 = sum(der_x.*der_y, 3);
G22 = sum(der_y.^2, 3) + 1;
A = sum(sum(sqrt(G11.*G22 - G12.^2)));
% A = sum(sum(G11.*G22 - G12.^2));

function der_x = calc_der_x(I)
der_x = (1/2) * (I(:,3:end,:) - I(:,1:end-2,:));

function der_y = calc_der_y(I)
der_y = (1/2) * (I(3:end,:,:) - I(1:end-2,:,:));
