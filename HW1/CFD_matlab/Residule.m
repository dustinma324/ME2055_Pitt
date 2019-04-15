clear all; close all; clc;
%% Residual Vs. Iteration
JX = [250 500 750 1000 1100 1200 1300 1400 1500]; %1600 1700 1800];
JY = [2.08923817 0.35758389 0.05613233 0.00818617 0.00449353 0.00192563...
    0.00096278 0.00032092 0.00016046];
GX = [250 500 750 1000];
GY = [0.782520 0.019531 0.000325 0];

figure();
hold on
plot(JX,JY,'-ro','LineWidth',2)
plot(GX,GY,'-bx','LineWidth',2)
hold off

ylabel('R = L2 Norm_{new} - L2 Norm_{old}')
xlabel('Number of Iteration')
axis auto
legend('Jacobi','Gauss-Seidel','Location','NorthWest')
grid on

%% Residual Vs. Mesh Size
JX = [2 3 4 5];
JY = [0.000078 0.000052 0.000078 0.000093];

GX = [2 3 4 5];
GY = [0.000080 0.000053 0.000080 0.000095];

figure();
hold on
loglog(JX,JY,'-ro','LineWidth',2)
loglog(GX,GY,'-bx','LineWidth',2)
hold off

ylabel('R = L2 Norm_{new} - L2 Norm_{old}')
xlabel('Mesh Refinment 2^{n}')
axis auto
legend('Jacobi','Gauss-Seidel','Location','NorthWest')
grid on

%% Heat Generation Vs. Mesh Quality
JX = [1 2 3 4 5];
JY = [2.4730 ];

GX = [1 2 3 4 5];
GY = [1956.77 2167.32 2215.14 2234.85 2245.86];

figure();
hold on
plot(JX,JY,'-ro','LineWidth',2)
plot(GX,GY,'-bx','LineWidth',2)
hold off

ylabel('Heat Generation Rate, q [W]')
xlabel('Mesh Refinment 2^{n}')
axis auto
legend('Jacobi','Gauss-Seidel','Location','NorthWest')
grid on
