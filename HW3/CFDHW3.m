clear all; close all; clc;

% ME 2055 - CFD
% Dustin (Ting-Hsuan) Ma
% Homework #3

%% User Define Method and Variables
tags = ["Explicit FTCS","Implicit FTCS","1st Order Upwind","MacCormack Scheme"];
method = 3;
MAXITER = 100;

nu = 1; % nu = 0 = invisid, nu > 0 = diffusive
U = 1;
L = 1;
eps = 1e-6;

if(nu == 0)
    fprintf("Invisid solver,");
elseif (nu > 0)
    fprintf("Diffusive solver,");
else
    fprintf("Not a valid viscosity input.\n");
end
fprintf(" using %s\n",tags(method))
fprintf("--------------------------------------------\n")

%% Memory Allocation Vlues
num = 11;
GhostLayer = 2;

%% Dependent Variables
dx = L/(num-1);
dt = 0.001;
period = L/U;
x = linspace(0,L,num);

%% Stability Parameters
CFL = U*dt/dx;
F = nu*dt/(dx*dx);
Re = U*dx/nu;

if (F > 0.5 || CFL > 1)
    fprintf("Please pick different dx and dt Values\n")
end
fprintf("CFL = %3.3f, F = %3.3f, Re = %3.3f, Period = %3.3f\n",CFL,F,Re,period)


%% Solvers
switch method
    case 1 % FTCS, Explicit
        % Initialization
        T = zeros(1,num+GhostLayer);
        Tnew = zeros(1,num+GhostLayer);
        T(2:end-1) = 1;
        
        % Boundary Conditions
        T = DirchletBC(T,method);
        T = NeumannBC(T,method);
        
        for iter = 1:1:MAXITER
            for i = 2:(numel(T)-1)
                Tnew(i) = -0.5*CFL*(T(i+1)-T(i-1))+F*(T(i-1)-2*T(i)+T(i+1))+T(i);
            end
            Tnew = DirchletBC(Tnew,method); %Applying BC Dirichlet
            Tnew = NeumannBC(Tnew,method); % Applying BC Neumann
            T = Tnew;   % Array swapping
        end
        plot(x,T(2:end-1),'-bs','LineWidth',2)
        
    case 2 % FTCS, Implicit
        % Initialization
        T = zeros(num,1);
        Tnew = zeros(num,1);
        T(1:end) = 1;
        
        % Boundary Conditions
        T = DirchletBC(T,method);
        T = NeumannBC(T,method);
        
        % Tridiagonal Matrix
        A = zeros(num,num);
        A(1,1) = 1;
        A(num,num) = 1;
        A(1,2) = 1; %(0.5*CFL+F);           % alpha
        A(num,num-1) = 1; %(0.5*CFL+F);    % gamma
        for i=2:num-1
            A(i,i-1) = (0.5*CFL-F);    % gamma
            A(i,i) = (1+2*F);           % beta
            A(i,i+1) = -(0.5*CFL+F);     % alpha
        end
        
        for count = 1:1:MAXITER
            Tnew = A\T;
            Tnew = DirchletBC(Tnew,method); %Applying BC Dirichlet
            Tnew = NeumannBC(Tnew,method); % Applying BC Neumann
            T = Tnew;   % Array swapping
        end
        plot(x,T,'-*b','LineWidth',2)
        
    case 3 % FT, 1st Order Upwind
        % Initialization
        T = zeros(1,num+GhostLayer);
        Tnew = zeros(1,num+GhostLayer);
        T(2:end-1) = 1;
        
        % Boundary Conditions
        T = DirchletBC(T,method);
        T = NeumannBC(T,method);
        
        for iter = 1:1:MAXITER
            for i = 2:(numel(T)-1)
                Tnew(i) = (CFL+F)*T(i-1)+(1-CFL-2*F)*T(i)+F*T(i+1);
            end
            Tnew = DirchletBC(Tnew,method); %Applying BC Dirichlet
            Tnew = NeumannBC(Tnew,method); % Applying BC Neumann
            T = Tnew;   % Array swapping
        end
        plot(x,T(2:end-1),'-bs','LineWidth',2)
    case 4 % MacCormack Scheme
        
end
%% Plotting Variables
xlabel("x")
ylabel("Temperature")
legend(tags(method),'location','northwest')
%axis([0 1 0 1.1])
grid on
grid minor

%% Functions

function T = DirchletBC(T,method)
switch method
    case 1
        T(2) = T(end-1);
    case 2
        T(1) = T(end);
    case 3
        T(2) = T(end-1);
    case 4
        
end
end

function T = NeumannBC(T,method)
switch method
    case 1
        T(end-2) = T(3);
        %T(3) = T(end-2); % This way gets rid of the skewness
    case 2
        T(end-1) = T(2);   
    case 3
        T(end-2) = T(3);
    case 4
        
end
end
