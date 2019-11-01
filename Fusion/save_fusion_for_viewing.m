function save_fusion_for_viewing(I,wl,bbl,save_file)
M = 4;
names = cell(1,M);
for j = 1:M
    names{j} = ['endmember ',num2str(j)];
end
[rows,cols,B] = size(I);

R_gt = zeros(M,B);
A_gt = zeros(rows,cols,M);
I_gt = A_gt2I_gt(A_gt);

rgb = uint8(sqrt(retrieve_rgb(I,wl))*255);
gt_colors = distinguishable_colors(M);
tick_labels = names;

save(save_file,'I','R_gt','A_gt','names','wl','rgb','I_gt',...
    'gt_colors','tick_labels','bbl');

end