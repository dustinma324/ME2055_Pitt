function [x, y, h, startx, starty] = meshGeneration(N,Dimension)
% defining dx and dy
h = Dimension(1)/(N-1);

% creating mesh
xg = linspace(0,Dimension(1),N);
yg = linspace(0,Dimension(2),N);
[x y] = meshgrid(xg,yg);

startx = 0:2*h:1;
starty = 0:2*h:1;
end