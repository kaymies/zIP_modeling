function [fignum, p_ellipse, p_circ] = drawStiffnessEllipse(K, fignum, lineSpecs)
%% DRAWSTIFFNESSELLIPSE Draw ellipse-circle representation of stiffness matrix
%
% Inputs:
% K - 2x2 matrix
% fignum - (optional) figure number on which to plot ellipse
% lineSpecs - (optional) struct specifying the following line parameters:
%              - color: color of plot lines
%              - width: width of plot lines
%              - style: style of plot lines ('-','--',etc.)
%
%
% Author: Rika Sugimoto-Dimitrova (rikasd@mit.edu)
% 2021-03-14
% 
% Modified:
% 2025-07-02 - cleaned up code for post-stroke zIP paper

if nargin > 1 && ~isempty(fignum)
    figH = figure(fignum);hold on;
else
    figH = figure;hold on;
    fignum = figH.Number;
end

if nargin < 3 || isempty(lineSpecs)
    linecolor = zeros(1,3);
    linewidth = 1.2;
    linestyle = '-';
else
    linecolor = lineSpecs.color;
    linewidth = lineSpecs.width;
    linestyle = lineSpecs.style;
end
  
K_S = (K + K')/2;
K_A = (K - K')/2;
[V_sym,D_sym] = eig(K_S);
D_anti = K_A(1,2);

if ~isreal(D_sym)
    warning('Complex eigenvalues');
    V_sym = real(V_sym);
    D_sym = real(D_sym);
end

[~,max_index] = max([D_sym(1,1),D_sym(2,2)]);
[~,min_index] = min([D_sym(1,1),D_sym(2,2)]);
axis1_x = [(-1)*V_sym(1,max_index)*D_sym(max_index,max_index), V_sym(1,max_index)*D_sym(max_index,max_index)];
axis1_y = [(-1)*V_sym(2,max_index)*D_sym(max_index,max_index), V_sym(2,max_index)*D_sym(max_index,max_index)];            
axis2_x = [(-1)*V_sym(1,min_index)*D_sym(min_index,min_index), V_sym(1,min_index)*D_sym(min_index,min_index)];
axis2_y = [(-1)*V_sym(2,min_index)*D_sym(min_index,min_index), V_sym(2,min_index)*D_sym(min_index,min_index)];

% Plot ellipse
a = norm([axis1_x(2),axis1_y(2)]);
b = norm([axis2_x(2),axis2_y(2)]);

th = 0:0.01:(2*pi);

x_sym = a.*cos(th);
y_sym = b.*sin(th);

r_anti = abs(D_anti(1,1));
x_anti = r_anti.*cos(th);
y_anti = r_anti.*sin(th);

th_rot = atan2(axis1_y(2), axis1_x(2));
R = [cos(th_rot) -sin(th_rot); sin(th_rot) cos(th_rot)];

ellipse_points = R*[x_sym;y_sym];

p_ellipse = plot(ellipse_points(1,:),ellipse_points(2,:),...
    linestyle,'Color',linecolor,'Linewidth',linewidth);
        
p_circ = plot(x_anti,y_anti,...
    linestyle,'Color',linecolor,'Linewidth',linewidth);

axis equal

xticks(yticks);

end