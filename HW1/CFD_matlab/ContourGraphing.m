clear all; close all; clc;
X = importdata('../xGrid.csv');
Y = importdata('../yGrid.csv');
T = importdata('../Temperature.csv');

for i=1:5
    for j = 9:11
        X(i,j) = nan;
        Y(i,j) = nan;
        T(i,j)= nan;
    end
end


%contourf(X,Y,T,'Showtext','on')
s=surf(X,Y,T);colorbar;s.FaceColor = 'interp';view(0,90);s.EdgeColor = 'none'
xlabel('X')
ylabel('Y')