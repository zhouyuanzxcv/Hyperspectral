function [data_mean,data_diff,md,sd] = bland_altman_x_gt(data1,data2)
% Function to generate Bland Altman plots. Barry Greene, September 2008
% Bland, J.M., Altman, D.G. 'Statistical methods for assessing agreement ...
% between two methods of clinical measurement'(1986) Lancet, 1 (8476), pp. 307-310.
%
% Inputs: data1: ground truth
%         data2: estimated
%
% Produces Bland Altman plot with mean difference and mean difference +/-
% 2*SD difference lines.

[m,n] = size(data1);
if(n>m)
    data1 = data1';
end

if(size(data1)~=size(data2))
    error('Data matrices must be the same size')
end

data_mean = data1;  % Mean of values from each instrument 
data_diff = data2 - data1;              % Difference between data from each instrument
md = mean(data_diff);               % Mean of difference between instruments 
sd = std(data_diff);                % Std dev of difference between instruments 

% figure;
plot(data_mean,data_diff,'ok','MarkerSize',3,'LineWidth',1);   % Bland Altman plot
hold on; x = [min(data_mean),max(data_mean)];
y = md*ones(1,2); plot(x,y,'-k'); % Mean difference line  
y = md+2*sd*ones(1,2); 
plot(x,y,'--k'); % Mean plus 2*SD line  
% text(x(2),y(2),'+2 SD'); 
y = md-2*sd*ones(1,2); 
plot(x,y,'--k'); % Mean minus 2*SD line   

% text(x(2),y(2),'-2 SD'); 
% grid on
% title('Bland Altman plot','FontSize',9)
% xlabel('Mean of two measures','FontSize',8)
% ylabel('Difference between two measures','FontSize',8)