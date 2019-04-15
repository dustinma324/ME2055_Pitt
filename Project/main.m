% ME 2055 - Computational Fluid Dyanics
% Dustin (Ting-Hsuan) Ma
clear all; close all; clc;

%% Mesh Properties
N = 60;
l = 1;
Dim = [l l]; % X, Y

[x y h startx starty] = meshGeneration(N,Dim);
%% Time Properties
%dt = 0.0000125; %re800
dt = 0.00015;

%% Flow Properties
Utop = 100;
nu = 0.1;

Re = Utop*l/nu;
CFL = Utop*dt/h;

disp("2D Lid Driven Cavity")
disp("=====================")
fprintf("Dimensions [%d %d]\nGrid [%d %d]\nRe = %4.2f\n",Dim(1),Dim(2),N,N,Re);
fprintf("CFL = %2.4f\n",CFL);
%% Simulation Properties
eps = 1e-4;

%% Memory Allocation
sF   = zeros(N,N);
w    = zeros(N,N);
u    = zeros(N,N);
v    = zeros(N,N);

%% Boundary Condition
[sF w] = boundaryCondition(sF,w,N,h,Utop);

%% Iteration/Plotting
err = 1;
iter = 0;
timeTotal = 0;

figure(1)    
view([0 0 1]);
colormap jet;


while err > eps
    %iter
    [sF w err] = FTCS_GS(sF,w,h,dt,N,nu,Utop);
    contourf(x,y,sF,10);
    pause(eps)
    iter = iter + 1;
    timeTotal = timeTotal + dt;
end

contourf(x,y,sF,15);
xlabel("x")
ylabel("y")
%title(['Stream Function, Re = ',num2str(Re)])

figure(2)
contourf(x,y,w,100);
xlabel("x")
ylabel("y")
%title(['Vorticity, Re = ',num2str(Re)])


mag = zeros(N,N);
[u v] = veloctiyBC (u,v,sF,Utop,h,N);
for i = 2:N-1
    for j = 2:N-1
        mag(j,i) = sqrt(u(j,i).^2+v(j,i).^2);
        u(j,i) = u(j,i)/mag(j,i);
        v(j,i) = v(j,i)/mag(j,i);
    end
end
figure(3)
streamline(x,y,u,v,startx,starty)
xlabel("x")
ylabel("y")
%title(['U,V, Re = ',num2str(Re)])

disp("SIMULATION COMPLETE")
fprintf("Steady State at iteration# %d, time = %4.2f(s)\n",iter,timeTotal)