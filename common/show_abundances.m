function [imhs, fig_h] = show_abundances(A2, m, n, title, plot_row, plot_col, options)
if nargin < 2
    if ndims(A2) ~= 3
        return;
    end
end

if nargin < 4
    title = 'Abundances';
end

if ndims(A2) == 3
    [m,n,M] = size(A2);
    A2 = reshape(A2, [m*n,M]);
end

endmember_num = size(A2,2);
if nargin < 5 || plot_row == 0 || plot_col == 0
    [row,col] = auto_plot_size(m,n,endmember_num);
else
    row = plot_row;
    col = plot_col;
end

if nargin < 7
    options = [];
end

update = 0;
imhs = [];
fig_h = [];
show_abundance_histogram = 0;
color_map = 'jet';

if nargin > 5 && isstruct(options)
    arg_set = fieldnames(options);
    for i = 1:length(arg_set)
        eval([arg_set{i},'=options.',arg_set{i},';']);
    end
end


if 0
    figure;
    for i = 1:endmember_num
        subplot(row,col,i);
        I = reshape(A2(:,i),[m n]);
        imshow(I);
        xlabel(['(',char(96+i),') endmember ',num2str(i)]);
    end
    colormap(gca,color_map);
else
    if ~update
        fig_h = figure('name',title);
        for i = 1:endmember_num
            subtightplot(row,col,i);
            I = reshape(A2(:,i),[m n]);
            imhs(i) = imshow(I);
            colormap(gca,color_map);
        end
        if show_abundance_histogram
            figure('name','Abundance histogram');
            for i = 1:endmember_num
                subplot(row,col,i);
                I = reshape(A2(:,i),[m n]);
                hist(I(:));
                xlim([0 1]);
            end
        end
    else
        for i = 1:endmember_num
            I = reshape(A2(:,i),[m n]);
            set(imhs(i),'CData',I);
        end
    end


%     margin = round(0.05*min(m,n));
%     I = ones(row*m + (row-1)*margin, col*n + (col-1)*margin,3,'double');
%     for i = 1:endmember_num
%         r = ceil(i/col);
%         c = mod(i-1,col) + 1;
%         I1 = reshape(A2(:,i),[m n]);
%         I2 = ind2rgb(uint8(round(I1*255)),colormap('jet'));
%         start_r = (r-1)*m+1 + (r-1)*margin;
%         start_c = (c-1)*n+1 + (c-1)*margin;
%         I(start_r:start_r+m-1,start_c:start_c+n-1,:) = I2;
%     end
%     figure;
%     imshow(I);
end
end

function [plot_rows,plot_cols] = auto_plot_size(im_rows,im_cols,endmember_num)
rows_candidate = (1:endmember_num);
vals = zeros(size(rows_candidate));
for i = 1:length(rows_candidate)
    rows = rows_candidate(i);
    cols = ceil(endmember_num / rows);
    height = rows * im_rows;
    width = cols * im_cols;
    vals(i) = abs(height - width) / mean([height,width]);
end
[~,ind] = min(vals);
plot_rows = rows_candidate(ind);
plot_cols = ceil(endmember_num / plot_rows);
end