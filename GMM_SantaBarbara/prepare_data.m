function [Is,wl,bbl,A_gt_total,gt_names,endmembers,endmembers_r] = prepare_data(dataset)
%PREPARE_DATA Summary of this function goes here
%   Detailed explanation goes here
% dataset can be a string "sba_im_16_lib_16" or "sba_syn_im_16_lib_16"
if isa(dataset, 'function_handle')
    [Is,wl,bbl,A_gt_total,gt_names,endmembers,endmembers_r] = dataset();
    return;
end

if strcmp(dataset(1:7), 'sba_syn')
    load('analysis/santabarbara_data_synthetic.mat');
    type = dataset(9:end);
elseif strcmp(dataset(1:7), 'sba_im_')
    if 0 
        load('analysis/santabarbara_data.mat');
    else
        % since github does not allow files larger than 50MB, the original
        % data is split into 6 parts
        filename = '../data/SantaBarbara/santabarbara_data_part_';
        for part_idx = (0:5)
            load([filename,num2str(part_idx)]);
        end
        I4s = [I4s_part1, I4s_part2, I4s_part3, I4s_part4, I4s_part5];
    end
    
    A_gt_total16 = A_gt_total;
    A_gt_total4 = A_gt_total;
    type = dataset(5:end);
end

switch type
    case 'im_16_lib_16'
        Is = I16s;
        endmembers = endmembers16;
        endmembers_r = endmembers16_r;
        A_gt_total = A_gt_total16;
    case 'im_16_lib_4'
        Is = I16s;
        endmembers = endmembers4;
        endmembers_r = endmembers4_r;
        A_gt_total = A_gt_total16;
    case 'im_4_lib_16'
        Is = I4s;
        endmembers = endmembers16;
        endmembers_r = endmembers16_r;
        A_gt_total = A_gt_total4;
    case 'im_4_lib_4'
        Is = I4s;
        endmembers = endmembers4;
        endmembers_r = endmembers4_r;
        A_gt_total = A_gt_total4;
    otherwise
end

end

