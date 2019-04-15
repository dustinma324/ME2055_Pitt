clear all; close all; clc;
% Graphing CFL vs Residual
dt = [0.006 0.005 0.004 0.002 0.001 0.0005 0.00005 0.000001];
CFL = [0.486 0.405 0.243 0.162 0.081 0.041 0.004 0.000081];

RE = [1.00 1.28 2.05 2.57 3.14 3.44 3.81 3.85];
RI = [1.32 1.54 2.18 2.63 3.16 3.45 3.81 3.85];

subplot(2,1,1)
hold on
%plot(CFL,RE,'-bs','LineWidth',2)
%plot(CFL,RI,'-r^','LineWidth',2)
loglog(CFL,RE,'-bs','LineWidth',2)
loglog(CFL,RI,'-r^','LineWidth',2)
hold off
grid on
grid minor
legend('Explicit','Implicit','Location','northeast')
title({'100 Iteration Optimization','dt = [6, 5, 4, 2, 1, 0.5, 0.05 0.001] E-3 [s]'})
xlabel('Courrant Number')
ylabel('Residual')

dx = [3 4 5 6 7 8 9 10]
RE = [5.6E-1 6.34E-1 7.06E-1 7.73E-1 8.35E-1 8.93E-1 9.48E-1 1.0];
RI = [1.01 9.96e-1 1.03 1.08 1.14 1.20 1.26 1.32];
subplot(2,1,2)
hold on
loglog(dx,RE,'-bs','LineWidth',2)
loglog(dx,RI,'-r^','LineWidth',2)
hold off
grid on
grid minor

title('100 Iteration Optimization')
xlabel('Num Points')
ylabel('Residual')

