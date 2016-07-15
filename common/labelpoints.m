function h = labelpoints (xpos, ypos, labels, varargin)
%  h = labelpoints (xpos, ypos, labels, position, buffer, adjust_axes, {parameters})
%
%   Given x and y position vectors (xpos, ypos) and given 
%     a vector of labels this script will label all data points 
%     and output the label handles in vector, h. 
%      
%   'xpos' and 'ypos' are required to be the same length. However, if all labels
%     fall along the same horizontal or vertical line, the function will accept 
%     a single number for that xpos or ypos (see examples).  
%
%   'labels' should be a cell or numerical array of the 
%     same length as xpos or ypos.  Alternatively it can be a singleton
%     that will be replicated for all labels (see examples).  
%
%   'position' (optional) describes the position of the labels
%     relative to their locations by entering one of the following
%     abbreviated compass directions in single quotes (default is 'NW').
%     N, S, E, W, NE, NW, SE, SW, center
%
%   'buffer' (optional) is a number between 0:1 that adds distance between
%     the label and the plotted point where 0 (default) is none and 1 is
%     1/10th of the axis.  Ignored for 'center' position.
%
%   'adjust_axes' (optional, default=0): depending on the positioning of labels, 
%     some may fall beyond the axis limits. adjust_axes = 1 will readjust xlim 
%     & ylim slightly so labels are not beyond axis limits.
%
%   This function includes three optional parameters to label only the outliers.
%     This can come in handy when there are many points but only the outliers
%     should be labeled. (may require stats toolbox)  (see examples)
%       'outliers_SD', N   -  will only label points that are greater than N 
%                             standard deviations from the median of xpos or ypos.
%       'outliers_Q', N    -  will only label points that are greater than N
%                             times the interquartile range of xpos or ypos.
%       'outliers_N', N    -  will calculate the distance of each point from the median 
%                             point(xpos,ypos) and will label the N-furthest points.
%                             Alternatively, N can be a decimal to label N% of the points. 
%
%   The following parameters may also be entered in any order
%     'FontSize', N -   font size of all labels
%     'Color', S    -   font color of all labels ('y' 'm' 'c' 'r' 'g' 'b' 'w' 'k')
%     'rotation', N -   will rotate the labels N degrees about label-center
%                       (positive is counterclockwise).
%
%
%   Examples:  
%      Fake Data:
%       x = 1:10;  y=rand(1,10); 
%       scatter(x,y)
%       labs = {'a' 'b' 'c' 'd' 'middle' 'f' 'g' 'h' 'i' 'last'};
%   
%      Label Examples
%       txt_h = labelpoints(x, y, labs);              
%       txt_h = labelpoints(x, y, labs, 'E', 0.15); 
%       txt_h = labelpoints(x, y, labs, 'E', 0.15, 1);
%       txt_h = labelpoints(x, y, labs, 'W', 0.15, 1, 'FontSize', 14, 'Color', 'r');
%       txt_h = labelpoints(x, y, labs, 'W', 0.15, 1, 'FontSize', 12, 'Color', 'm', 'rotation', 45);
% 
%      Also works for     
%       labs = [1:1:10];            
%       labs = {'Sofia' '' '' '' 'Bucharest' '' '' 'Belgrade' '' 'Ankara'}
%       labs = '*';              
%       labs = 'string';
%
%      When all labels share same xpos or ypos
%         boxplot(1:10)
%         labelpoints(0.8, [3, 5.5, 8], {'25%' '50%' '75%'}, 'center');
%
%      Outlier Examples
%       Fake Data:
%         x = [rand(1,30), rand(1,8)*2];
%         y = [rand(1,30), rand(1,8)*2];
%         scatter(x, y)
%         labs = 1:38;
%        labelpoints(x, y, labs, 'N', 0.1, 1, 'outliers_N', 5);                  %will label 5 furthest points from median
%        labelpoints(x, y, labs, 'N', 0.1, 1, 'outliers_N', 0.1, 'Color', 'r');  %will label 10% of furthest points from median    
%        labelpoints(x, y, labs, 'N', 0.1, 1, 'outliers_SD', 2);                 %will label all points > 2 SD from median
%        labelpoints(x, y, labs, 'N', 0.1, 1, 'outliers_Q', 1.5);                %will label points greater that 1.5 x IQR
% 
%   Alternative use: 
%     Density Distributions:
%       x = randn(1,100); y = 1:100;  
%       scatter(x,y)
%       labs = '|';
%       txt_h = labelpoints(x, 8, labs, 'center');
%
%     Single Labeling
%       x = 2004:2013;  y=rand(1,10); 
%       plot(x,y,'-o')
%       labs = 'acquisition';
%       labelpoints(x(3), y(3), labs, 'N', 0.2, 1);
%       labelpoints(2008.5, min(ylim), {['labelpoints.m   ', datestr(now, 'mm/dd/yy')]}, 'N', 0.3, 0, 'fontsize', 12, 'color', 'm');
%
%     Use labels instead of markers
%       x = randn(1,15); y = randn(1,15);
%       labs = char('a'+(1:15)-1)';
%       labelpoints(x, y, labs, 'center', 0, 1, 'color', 'b');
% 
%  140331 v.1 
%  141115 v.2
% Copyright (c) 2014, Adam Danz
%All rights reserved

% source: http://www.mathworks.com/matlabcentral/fileexchange/46891-labelpoints

% Changes history
%   11/02/14    if labels are entered as char, they are convered to cell
%   11/03/14    varargin added to accomodate auto_outlier feature
%   11/03/14    convert inputs to row vector if they are column vectors
%   11/04/14    now position in put is not case sensitive
%   11/04/14    now 1 label can be entered to label multiple points
%   11/04/14    now fontsize and color can be specified by params
%   11/05/14    changed 'outlier' to 'outlier_SD' and added 'outlier_Q'
%   11/07/14    added option to rotate test
%   11/15/14    added 'outliers_N' input option and cleaned up text rotation section.
%   11/19/14    curr_extent is not always a cell.  fixed.  
%   11/20/14    when outliers_N is selected N is min(N, lenght(xpos));
%   11/21/14    removes entire point and label when there is an 'inf' value


%%
% Check Class of 'labels'
    %If 'labels' are numberical, convert to cell
    if isnumeric(labels) == 1
        labels = num2cell(labels); 
    end

    % if 'labels' are char, convert to cell
    if ischar(labels)
        labels = cellstr(labels);
    end
    
% if all labels share the same xpos or ypos (only 1 value entered in 1 of the position vectors)
    if length(xpos)==1 && length(ypos)>1
        xpos = repmat(xpos, size(ypos));
    elseif length(ypos)==1 && length(xpos)>1
        ypos = repmat(ypos, size(xpos));
    end
    
% if only one lable is entered for all points, replicate it
    if length(labels)==1 && length(xpos) > 1
        labels = repmat(labels, [1, length(xpos)]);
    end
    
% ensures xpos, ypos, and labels are all row vectors 
    if iscolumn(xpos);      xpos = xpos';       end
    if iscolumn(ypos);      ypos = ypos';       end
    if iscolumn(labels);    labels = labels';   end

%check that x, y, and labels are same length
    if isequal(length(xpos), length(ypos), length(labels)) == 0
        error('xpos, ypos, and labels must all be the same length unless using one input for labels.')
    end
    
%if an 'inf' value is entered, this will remove that entire point and label
    xinf = find(xpos==inf);
    yinf = find(ypos==inf);
    findinf = [xinf yinf];
    if ~isempty(findinf)
        xpos(findinf)=[];
        ypos(findinf)=[];
        labels(findinf) = [];
    end       

%Validate inputs and optional parameters
    validPositions = {'N' 'NE' 'E' 'SE' 'S' 'SW' 'W' 'NW' 'center'};
    checkPosition = @(x) any(validatestring(x, validPositions));

    validColors = {'y' 'm' 'c' 'r' 'g' 'b' 'w' 'k'};
    checkColors = @(x) any(validatestring(x, validColors));

    p = inputParser;
    p.FunctionName = mfilename;
    addRequired(p, 'xpos', @isnumeric);
    addRequired(p, 'ypos', @isnumeric);
    addRequired(p, 'labels');
    addOptional(p, 'position', 'NW', checkPosition);
    addOptional(p, 'buffer', 0, @isnumeric);
    addOptional(p, 'adjust_axes', 0, @isnumeric);

    addParamValue(p, 'outliers_SD', 3, @isnumeric);
    addParamValue(p, 'outliers_Q', 1.5, @isnumeric);
    addParamValue(p, 'outliers_N', 1, @isnumeric);

    addParamValue(p, 'FontSize', 10, @isnumeric);
    addParamValue(p, 'Color', 'k', checkColors);
    addParamValue(p, 'rotation', 0, @isnumeric);
    parse(p, xpos, ypos, labels, varargin{:})

mfile = [mfilename,'.m'];

%calculate buffer
    a = axis/10;% I've somewhat arbitrarily divided by 10 to make 'buffer' more sensitive
    u1 = 0;     %x offset
    u2 = 0;     %y offset

%assign position
    switch upper(p.Results.position) 
        case 'E',       va = 'middle'; ha = 'left';         u1 = a(2)-a(1);         
        case 'W',       va = 'middle'; ha = 'right';        u1 = (a(2)-a(1))*-1;
        case 'N',       va = 'bottom'; ha = 'center';       u2 = a(4)-a(3);
        case 'S',       va = 'top'; ha = 'center';          u2 = (a(4)-a(3))*-1;
        case 'NE',      va = 'bottom'; ha = 'left';         u1 = (a(2)-a(1))/2;     u2 = (a(4)-a(3))/2;
        case 'NW',      va = 'bottom'; ha = 'right';        u1 = (a(2)-a(1))*-0.5;  u2 = (a(4)-a(3))/2;
        case 'SE',      va = 'top'; ha = 'left';            u1 = (a(2)-a(1))/2;     u2 = (a(4)-a(3))*-0.5;
        case 'SW',      va = 'top'; ha = 'right';           u1 = (a(2)-a(1))*-0.5;  u2 = (a(4)-a(3))*-0.5;
        case 'CENTER',  va = 'middle'; ha = 'center';    
    end

%Factor in buffer
    u1 = u1*p.Results.buffer;
    u2 = u2*p.Results.buffer;


%If outliers parameters are selected
    if sum(strcmp(varargin, 'outliers_SD')) == 1 || sum(strcmp(varargin, 'outliers_Q')) == 1 || sum(strcmp(varargin, 'outliers_N')) == 1
        if sum(strcmp(varargin, 'outliers_SD')) == 1
            outlier_idx = logical(abs(xpos - median(xpos)) > p.Results.outliers_SD*std(xpos)  |  abs(ypos - median(ypos)) > p.Results.outliers_SD*std(ypos)); %index of outliers
            
        elseif sum(strcmp(varargin, 'outliers_Q')) == 1
            xbounds = [prctile(xpos,25) - p.Results.outliers_Q * iqr(xpos) , prctile(xpos, 75) + p.Results.outliers_Q * iqr(xpos)];   %[lower upper] bounds of outliers
            ybounds = [prctile(ypos,25) - p.Results.outliers_Q * iqr(ypos) , prctile(ypos, 75) + p.Results.outliers_Q * iqr(ypos)];   %[lower upper] bounds of outliers
            outlier_idx = logical(ypos<ybounds(1) | ypos>ybounds(2) |  xpos<xbounds(1) | xpos>xbounds(2));
            
        elseif sum(strcmp(varargin, 'outliers_N')) == 1
            if p.Results.outliers_N<1;  
                N = round(length(xpos) * p.Results.outliers_N); 
            else
                N = min(p.Results.outliers_N, length(xpos));        %ensures that user cannot label more outliers than coordinates.
            end
             medianpoint = repmat([median(xpos) median(ypos)], [length(xpos),1]);
             paired = horzcat(xpos', ypos');
             distances = (((medianpoint(:,1)-paired(:,1)).^2)  +  ((medianpoint(:,2)-paired(:,2)).^2)).^(1/2);       %all distances from median
             [~, idx] = sort(distances, 'descend');
             outlier_idx = false(1,length(xpos));
             outlier_idx(idx(1:N))=1; 
        end
        
        xpos = xpos(outlier_idx);               
        ypos = ypos(outlier_idx);
        labels = labels(outlier_idx);

        if any(outlier_idx) == 0;           %dispay msg if there are no outliers to label
            disp(['There are no outliers to label in ', mfile,'.'])
            disp('Change outlier value for less sensitivity; See help file.'); 
        end
    end

%Label points
    h = text(xpos+u1 , ypos+u2, labels, 'VerticalAlignment',va, 'HorizontalAlignment',ha, 'FontSize', p.Results.FontSize, 'color', p.Results.Color);
    
%Rotate text if specified
    if sum(strcmp(varargin, 'rotation')) == 1               %if rotation parameter is called in inputs
        xl = xlim;      yl = ylim;                          %In case text rotation auto adjusts axes.
        curr_extent = get(h, 'extent');                     %Need to store current center point of all labels since text rotation relocates position
        if iscell(curr_extent); cell2mat(curr_extent); end
        hold on
        curr_position = [curr_extent(:,1)+(curr_extent(:,3)/2),curr_extent(:,2)+(curr_extent(:,4)/2)];          %uses extent to locate center of label
        set(h, 'rotation', p.Results.rotation, 'VerticalAlignment','middle', 'HorizontalAlignment','center');  	%note: text rotation changes alignment which is why they need to be centered back to specifications.
        for i = 1:length(h)                                 %after rotation, reposition lables back to desired location 
            set(h(i), 'position', curr_position(i,:))
        end
        set(gca, 'xlim', xl); set(gca, 'ylim', yl);         %In case text rotation auto adjusts axes.
    end     
    
%Determine if any labels go beyond axis limits and adjust if desired  (adjust_axes = 0 or 1)
    if p.Results.adjust_axes == 1   &&   ~isempty(h)    
        x_adj = sign(u1+0.0000001);                 %the addition is to avoid '0'
        y_adj = sign(u2+0.0000001);                 %the addition is to avoid '0'

        labelextent = get(h, 'extent');
        if isequal(class(labelextent),'cell')
           labelextent = cat(1, labelextent{:});
        end
        xl = xlim;      yl = ylim;
        lablimX = [min(labelextent(:,1)), max(labelextent(:,1)+(labelextent(:,3).*x_adj))] +u1;
        lablimY = [min(labelextent(:,2)), max(labelextent(:,2)+(labelextent(:,4).*y_adj))] +u2;

        xlim([min(min(xl), min(lablimX)), max(max(xl), max(lablimX))])
        ylim([min(min(yl), min(lablimY)), max(max(yl), max(lablimY))])
    end
     
end


%% Notes
% a video (not mine) explaining this method:  http://blogs.mathworks.com/videos/2012/05/30/how-to-label-a-series-of-points-on-a-plot-in-matlab/
% Text properties :  http://www.mathworks.com/help/matlab/ref/text_props.html
% info on input parsing:  http://www.mathworks.com/help/matlab/matlab_prog/parse-function-inputs.html
%                   and   http://www.mathworks.com/help/matlab/ref/inputparser-class.html
% Outlier info : https://docs.oracle.com/cd/E17236_01/epm.1112/cb_statistical/frameset.htm?ch07s02s10s01.html

