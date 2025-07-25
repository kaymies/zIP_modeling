function [p_mean, p_SD] = plotSDArea(x,y_mean,y_SD,ifig,option,lineSpecs)
%% PLOTSDAREA Plots a shaded area showing 1 SD (or SE) from mean.
% 
% Inputs:
% x      - independent variable
% y_mean - dependent variable mean value
% y_SD   - dependent variable standard deviation (SD) or standard error (SE)
% ifig   - figure number to plot on; if empty, create a new figure
% option - plotting options:
%          0: plot mean+/-SD as a shaded area only (default)
%          1: plot mean as a line only
%          2: plot mean+/-SD area and mean line
% lineSpecs - line specifications
%           color: 1x3 vector for color
%           width: line width of mean
%           alpha: transparency of SD area
%
% Outputs:
% p_mean - plot handle for mean line
% p_SD   - plot handle for SD or SE area
%
%
% Author: Rika Sugimoto-Dimitrova (rikasd@mit.edu)
% 2021-03-19
% 
% Modified:
% 2021-04-12 - added plotting options
% 2022-03-02 - added lineSpecs
% 2025-07-02 - cleaned up code for post-stroke zIP paper

p_mean = [];
p_SD = [];

if nargin > 3 && ~isempty(ifig)
    figure(ifig);hold on;
else
    figure;hold on;
end

if nargin < 5 || isempty(option)
    if nargin > 2 && ~isempty(y_SD)
        option = 0;
    else
        option = 1;
    end
end

if nargin < 6 || isempty(lineSpecs)
    lineSpecs.color = [0 0 0];
    lineSpecs.width = 1.2;
    lineSpecs.alpha = 0.5;
    lineSpecs.style = '-';
end

if option == 0 || option == 2
    
    if size(x,1) > size(x,2)
        x = x';
        y_mean = y_mean';
        y_SD = y_SD';
    end
    
    if size(y_SD,1) == 2
        if y_SD(1,1) < y_SD(2,1)
            lb = y_SD(1,:);
            ub = y_SD(2,:);
        else
            lb = y_SD(2,:);
            ub = y_SD(1,:);   
        end
    else   
        lb = y_mean-y_SD;
        ub = y_mean+y_SD;
    end
    
    for i = 1:size(lb,1)
        p_SD = fill([x,flip(x)],[ub(i,:),flip(lb(i,:))],lineSpecs.color,...
            'FaceAlpha',lineSpecs.alpha,'EdgeAlpha',0);
    end
    
    if option == 2
        ax = gca;
        p_mean = plot(x,y_mean,lineSpecs.style,'Linewidth',lineSpecs.width,...
            'Color',lineSpecs.color);
    end
    
else
    p_mean = plot(x,y_mean,'Linewidth',lineSpecs.width);
end

end