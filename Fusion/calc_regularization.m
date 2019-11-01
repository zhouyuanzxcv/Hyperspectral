function L = calc_regularization(I1, options)
% eta = parse_param(options,'eta',0.05);
k = parse_param(options,'K_neighbor',3);
% step = parse_param(options,'window_step',10);
% window = parse_param(options,'window_size',30);
% tau = parse_param(options,'tau',10);
r2 = parse_param(options,'r2',15);

% dummy_beta = 0.01; % not used since 0 is input as beta2
% [L,~] = calc_Laplacians(I1, eta, dummy_beta, 0);
if strcmp(options.regularization, 'LLMM')
%     L = L + tau*calc_manifold_LLE_constraint(I1,k,step,window);
    L = calc_manifold_LLE_constraint_radius(I1,k,1) + ...        
        calc_manifold_LLE_constraint_radius(I1,k,r2);
end

% L = L + calc_local_LMM_constraint(I1);
% L = L + calc_linear_operator(I1);
% L = L + calc_manifold_constraint(I1,5,eta);

function L = calc_manifold_LLE_constraint_radius(I,k,radius)
% add random noise for removing itself in distance computation
I = I + randn(size(I))*(1e-9*max(I(:)));

[rows,cols,b] = size(I);
N = rows*cols;
Y = reshape(I, [N,b]);
epsilon = 1e-4;
% epsilon = epsilon / k * b; % to adjust to the scales of two terms

J2 = cell(1,N);
A2 = cell(1,N);

[js_1,is_1] = meshgrid((1-ceil(radius):1+ceil(radius)),(1-ceil(radius):1+ceil(radius)));
js_1 = js_1(:);
is_1 = is_1(:);
% remove radius that is larger than preset
r_outside_inds = sum((js_1-1).^2 + (is_1-1).^2, 2) > radius^2;
js_1(r_outside_inds) = [];
is_1(r_outside_inds) = [];

radius1 = ceil(radius);
for i = 1:rows
    for j = 1:cols
        x_y = [js_1 + (j-1), is_1 + (i-1)];
        % remove points outside the image domain
        if ~(i > radius1 && j > radius1 && i < rows - radius1 && j < cols - radius1)
            x_y(x_y(:,1) < 1 | x_y(:,2) < 1 | x_y(:,1) > cols | x_y(:,2) > rows, :) = [];
        end
        inds = sub2ind([rows,cols], x_y(:,2), x_y(:,1));
        curr = sub2ind([rows,cols], i, j);
        dist2 = sum((Y(inds,:) - repmat(Y(curr,:),[length(inds),1])).^2, 2);
        [d,idx] = mink(dist2, k+1);
        % remove the point itself in the distance computed closest points
        nhbr_inds = inds(idx(2:end));
        J2{curr} = [nhbr_inds',curr];
    end
end

for i = 1:length(J2)
    inds = J2{i}(1:end-1); % (1:k)
    k1 = length(inds);
    curr = J2{i}(end); % k+1
    Z = repmat(Y(curr,:),[k1,1]) - Y(inds,:);
    C = Z*Z' + epsilon*eye(k1);
    wn = C\ones(k1,1);
    wn = wn/sum(wn);
    A2{i} = [wn;-1]';
end

L = cellarr2sparse(J2,A2,N,N);
L = L'*L;

% 
% function L = calc_manifold_LLE_constraint(I,k,step,window)
% % add random noise for removing itself in distance computation
% I = I + randn(size(I))*(1e-9*max(I(:)));
% 
% [rows,cols,b] = size(I);
% N = rows*cols;
% Y = reshape(I, [N,b]);
% epsilon = 1e-4;
% 
% J2 = cell(1,N);
% A2 = cell(1,N);
% 
% extra = (window - step) / 2;
% 
% 
% for i = 1:step:rows
%     for j = 1:step:cols
%         i_end = min(i + step - 1, rows);
%         j_end = min(j + step - 1, cols);
%         [js,is] = meshgrid((j:j_end),(i:i_end));
%         ind_J2 = sub2ind([rows,cols], is(:), js(:));
%         Y_center = Y(ind_J2,:);
%         % select window
%         win_row_start = max(i-extra, 1);
%         win_col_start = max(j-extra, 1);
%         rows_sel = win_row_start : min(i+step-1+extra,rows);
%         cols_sel = win_col_start : min(j+step-1+extra,cols);
%         I_sel = I(rows_sel, cols_sel, :);
%         % find nearest neighbors
%         Y_sel = reshape(I_sel, [size(I_sel,1)*size(I_sel,2), size(I_sel,3)]);
%         [Dist,ind_nearest] = pdist2(Y_sel,Y_center,'euclidean','Smallest',k+1);
%         ind_nearest(1,:) = []; % remove itself in distance computation
%         [rows_nearest,cols_nearest] = ind2sub([size(I_sel,1), ...
%             size(I_sel,2)], ind_nearest(:)); 
%         endm_y_x = [rows_nearest, cols_nearest];        
%         % convert endmember position to index based on the whole image
%         endm_y_x1 = endm_y_x - 1 + ...
%             repmat([win_row_start,win_col_start], [size(endm_y_x,1),1]);
%         endm_ind = sub2ind([rows,cols], endm_y_x1(:,1), endm_y_x1(:,2));
%         J2_val = [reshape(endm_ind, [k,length(ind_J2)])', ind_J2];
%         J2_val = mat2cell(J2_val, ones(1,size(J2_val,1)), size(J2_val,2));
%         J2(ind_J2) = J2_val;
%     end
% end
% 
% for i = 1:length(J2)
%     inds = J2{i}(1:k);
%     curr = J2{i}(k+1);
%     Z = repmat(Y(curr,:),[k,1]) - Y(inds,:);
%     C = Z*Z' + epsilon*eye(k);
%     wn = C\ones(k,1);
%     wn = wn/sum(wn);
%     A2{i} = [wn;-1]';
% end
% 
% L = cellarr2sparse(J2,A2,N,N);
% L = L'*L;


%% The following codes are obsolete
% function L = calc_manifold_constraint(I,k,eta)
% [rows,cols,B] = size(I);
% N = rows*cols;
% Y = reshape(I, [rows*cols,B]);
% eta = eta*sqrt(B);
% 
% [A,~,~] = calc_adjacency_graph_cached(Y,k);
% [I,J,Z] = find(A);
% Z = exp(-Z.^2/(2*eta^2));
% W = sparse(I,J,Z,N,N);
% 
% D = sparse((1:N),(1:N),sum(W),N,N);
% L = D - W;
% 
% 
% function L = calc_local_LMM_constraint(I)
% % Do linear unmixing locally and apply abundances to force linear relation
% [rows,cols,b] = size(I);
% N = rows*cols;
% 
% J2 = cell(1,N);
% A2 = cell(1,N);
% 
% step = 30;
% window = step*3;
% num_endm = 6;
% 
% extra = (window - step) / 2;
% 
% for i = 1:step:rows
%     for j = 1:step:cols
%         i_end = min(i + step - 1, rows);
%         j_end = min(j + step - 1, cols);
%         [js,is] = meshgrid((j:j_end),(i:i_end));
%         ind_J2 = sub2ind([rows,cols], is(:), js(:));
%         % select window
%         win_row_start = max(i-extra, 1);
%         win_col_start = max(j-extra, 1);
%         rows_sel = win_row_start : min(i+step-1+extra,rows);
%         cols_sel = win_col_start : min(j+step-1+extra,cols);
%         I_sel = I(rows_sel, cols_sel, :);
%         % do linear unmixing
%         [endm_y_x,abund] = local_linear_unmixing(I_sel, num_endm);
%         % convert endmember position to index based on the whole image
%         endm_y_x1 = endm_y_x - 1 + ...
%             repmat([win_row_start,win_col_start], [size(endm_y_x,1),1]);
%         endm_ind = sub2ind([rows,cols], endm_y_x1(:,1), endm_y_x1(:,2));
%         J2_val = [repmat(endm_ind', [length(ind_J2),1]), ind_J2];
%         J2_val = mat2cell(J2_val, ones(1,size(J2_val,1)), size(J2_val,2));
%         J2(ind_J2) = J2_val;
%         % find center area to be unmixed
%         rows_center_local = is - win_row_start + 1;
%         cols_center_local = js - win_col_start + 1;
%         center_inds_local = sub2ind([size(I_sel,1),size(I_sel,2)], ...
%             rows_center_local(:), cols_center_local(:));
%         % assign values in sparse matrix
%         A2_val = abund(center_inds_local,:);
%         A2_val = [A2_val,ones(size(A2_val,1),1)*-1];
%         A2_val = mat2cell(A2_val, ones(1,size(A2_val,1)), size(A2_val,2));
%         A2(ind_J2) = A2_val;
%     end
% end
% 
% L = cellarr2sparse(J2,A2,N,N);
% L = L'*L;
% 
% 
% 
% function [endm_y_x, A] = local_linear_unmixing(I_sel, num_endm)
% [rows,cols,b] = size(I_sel);
% N = rows*cols;
% lambda = sqrt(b)*0.01;
% epsilon = 1e-4;
% 
% Y_sel = reshape(I_sel, [rows*cols, b]);
% 
% % try endmember extraction algorithm
% [idx, C] = kmeans(Y_sel, num_endm);
% endm_indices = zeros(num_endm, 1);
% for i = 1:num_endm
%     ind_class_i = find(idx == i);
%     Y1 = Y_sel(ind_class_i,:);
%     [~,ind_min] = min(sum((Y1 - repmat(C(i,:), [size(Y1,1),1])).^2, 2));
%     endm_indices(i) = ind_class_i(ind_min);
% end
% % [endm_indices] = PPI(Y_sel, num_endm, 1000);
% % M = Y_sel(endm_indices,:);
% % 
% % % solve for abundances
% % Y1 = [Y_sel, ones(N,1)*lambda];
% % M1 = [M, ones(size(M,1),1)*lambda];
% % 
% % A = Y1*M1'*inv(M1*M1'+epsilon*eye(size(M1,1)));
% A = zeros(size(Y_sel,1), num_endm);
% for i = 1:num_endm, A(idx==i,i) = 1; end
% 
% % convert indices to row, col
% [I,J] = ind2sub([rows,cols], endm_indices);
% endm_y_x = [I,J];
% 
% 
% 
% function L = calc_linear_operator(I)
% % Use LLE in a 8-neighborhood graph
% % calculate neighborhood based on 8-neighbor system
% epsilon = 0.0001;
% [rows,cols,B] = size(I);
% Y = reshape(I, [rows*cols,B]);
% 
% N = rows * cols;
% N1 = (rows-2) * (cols-2);
% J2 = cell(1,N1);
% A2 = cell(1,N1);
% ind_J2 = 1;
% for j = 2:cols-1
%     for i = 2:rows-1
%         y = [i-1,i,i+1,i-1,i+1,i-1,i,i+1];
%         x = [j-1,j-1,j-1,j,j,j+1,j+1,j+1];
% %         y = [i-1,i,i+1,i];
% %         x = [j,j+1,j,j-1];
%         K = length(y);
%         inds = sub2ind([rows,cols],y,x);
%         curr = sub2ind([rows,cols],i,j);
%         Z = repmat(Y(curr,:),[K,1]) - Y(inds,:);
%         C = Z*Z' + epsilon*eye(K);
%         wn = C\ones(K,1);
%         wn = wn/sum(wn);
%         J2{ind_J2} = [inds,curr]';
%         A2{ind_J2} = [wn;-1];
%         ind_J2 = ind_J2 + 1;
%     end
% end
% 
% L = cellarr2sparse(J2,A2,N1,N);
% L = L'*L;

