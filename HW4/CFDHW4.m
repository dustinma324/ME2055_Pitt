% ME2055 - CFD Homework #4
% Boundary Layer Flow over Flat Plate
% Dustin (Ting-Hsuan) Ma

clear all; clc; close all;
%% Pick Solving Method
fprintf("1: Fourth Order Runge-Kutta\n2: Forward Euler\n")
n = input('Enter a number : ');
fprintf("===================================================\n")
switch n        
    case 1  % Fourth Order Runge-Kutta Method
        fprintf("Using Fourth Order Runge-Kutta Method\n");
        % Change These Values
        deltaEta = 0.2;
        y0guess1 = [0 0 0.1];
        y0guess2 = [0 0 1];

        % Breaking 3rd Order down to 3 coupled 1st order ODE
        f1 = @(y1,y2,y3) y2;
        f2 = @(y1,y2,y3) y3;
        f3 = @(y1,y2,y3) -0.5*y1*y3;
                
        % Initial Values
        h = deltaEta;
        eta = 0:deltaEta:10;
        
        % Convergence Criteria
        iter = 1;
        eps = 1e-14;
        err = 1;
        while err > eps
            if iter == 1
                y1(1) = y0guess1(1); y2(1) = y0guess1(2); y3(1) = y0guess1(3);
                for i = 1:(numel(eta)-1)
                    [output] = RungeKutta(f1,f2,f3,y1,y2,y3,h);
                    y1(i+1) = output(end,1);
                    y2(i+1) = output(end,2);
                    y3(i+1) = output(end,3);
                end
                m1 = y2(end);
                % =========================================================
                y1 = []; y2 = []; y3 = [];
                y1(1) = y0guess2(1); y2(1) = y0guess2(2); y3(1) = y0guess2(3);
                for i = 1:(numel(eta)-1)
                    [output] = RungeKutta(f1,f2,f3,y1,y2,y3,h);
                    y1(i+1) = output(end,1);
                    y2(i+1) = output(end,2);
                    y3(i+1) = output(end,3);
                end
                m2 = y2(end);
            else
                y1 = []; y2 = []; y3 = [];
                y0guess2(3) = (y0guess2(3)-y0guess1(3))/(m2-m1)*(1.0-m1)+y0guess1(3);
                y1(1) = y0guess2(1); y2(1) = y0guess2(2); y3(1) = y0guess2(3);
                for i = 1:(numel(eta)-1)
                    [output] = RungeKutta(f1,f2,f3,y1,y2,y3,h);
                    y1(i+1) = output(end,1);
                    y2(i+1) = output(end,2);
                    y3(i+1) = output(end,3);
                end
                m2 = y2(end);
                err = abs(1.0 - m2);
            end
            iter = iter + 1;
        end
        
        plot(eta,y1,'--k',eta,y2,'-r',eta,y3,'-.b','LineWidth',2)
        legend("f","fp","fpp")
        limit = zeros(1,numel(eta))+1;
        %plot(eta,limit,'--k',eta,y2,'-r','LineWidth',2)
        %ylim([0 1.2]);
        %ylim([0 5]);
        xlabel('\eta', 'fontsize', 16);
        %ylabel("f'=u/U_{\infty}",'fontsize',16);
        ylabel("f,fp,fpp",'fontsize',16);
        
        fprintf("f'' = %f\n",y3(1));
        
    case 2  % Forward Euler Mathod
        fprintf("Using Forward Euler Method\n");
        % Change These Values
        deltaEta = 0.2;
        y0guess1 = [0 0 0.1];
        y0guess2 = [0 0 1];

        % Breaking 3rd Order down to 3 coupled 1st order ODE
        f1 = @(y1,y2,y3) y2;
        f2 = @(y1,y2,y3) y3;
        f3 = @(y1,y2,y3) -0.5*y1*y3;
                
        % Initial Values
        h = deltaEta;
        eta = 0:deltaEta:10;
        
        % Convergence Criteria
        iter = 1;
        eps = 1e-14;
        err = 1;
        while err > eps
            if iter == 1
                y1(1) = y0guess1(1); y2(1) = y0guess1(2); y3(1) = y0guess1(3);
                for i = 1:(numel(eta)-1)
                    [output] = FEuler(f1,f2,f3,y1,y2,y3,h);
                    y1(i+1) = output(end,1);
                    y2(i+1) = output(end,2);
                    y3(i+1) = output(end,3);
                end
                m1 = y2(end);
                % =========================================================
                y1 = []; y2 = []; y3 = [];
                y1(1) = y0guess2(1); y2(1) = y0guess2(2); y3(1) = y0guess2(3);
                for i = 1:(numel(eta)-1)
                    [output] = FEuler(f1,f2,f3,y1,y2,y3,h);
                    y1(i+1) = output(end,1);
                    y2(i+1) = output(end,2);
                    y3(i+1) = output(end,3);
                end
                m2 = y2(end);
            else
                y1 = []; y2 = []; y3 = [];
                y0guess2(3) = (y0guess2(3)-y0guess1(3))/(m2-m1)*(1.0-m1)+y0guess1(3);
                y1(1) = y0guess2(1); y2(1) = y0guess2(2); y3(1) = y0guess2(3);
                for i = 1:(numel(eta)-1)
                    [output] = FEuler(f1,f2,f3,y1,y2,y3,h);
                    y1(i+1) = output(end,1);
                    y2(i+1) = output(end,2);
                    y3(i+1) = output(end,3);
                end
                m2 = y2(end);
                err = abs(1.0 - m2);
            end
            iter = iter + 1;
        end
        
        plot(eta,y1,'--k',eta,y2,'-r',eta,y3,'-.b','LineWidth',2)
        legend("f","fp","fpp")
        limit = zeros(1,numel(eta))+1;
        %plot(eta,limit,'--k',eta,y2,'-r','LineWidth',2)
        %ylim([0 1.2]);
        xlabel('\eta', 'fontsize', 16);
        ylabel("f'=u/U_{\infty}",'fontsize',16);
        
        fprintf("f'' = %f\n",y3(1));
                
    otherwise
        fprintf("Invalid input, please select a solver method\n");
        return;
end

%% Post processing Graphs
for i = 1:numel(y2)
    err(i) = (1-y2(i))^2;
end
figure(2);
loglog(eta,err,'LineWidth',2);
title("Residual Vs. \eta")
xlabel('log(\eta)')
ylabel('log(error)')

%% Fric Coeff and Drag
nu = 1;
rho = 1;
Uinf = 1;
mu = rho*nu;
L = 10;
x = 0:deltaEta:L;
for i = 1:numel(x)
    REx = Uinf*x(i)/nu;
    tau(i) = y3(1)*mu*Uinf*sqrt(Uinf/(nu*x(i)));
    cf(i) = tau(i)/((rho/2)*Uinf^2);
    cfREx(i) = 0.664/sqrt(REx);
end

syms w;
func = @(w) 0.332*mu*Uinf^(3/2)*nu^(-1/2)*w.^(-1/2);
D = integral(func,0,L);
cfm = D/(0.5*rho*Uinf^2*L);
cfmREl = 1.32/sqrt(Uinf*L/nu);

relerr = abs(cfmREl-cfm)/cfmREl * 100;
fprintf("cfmREl = %f, cfm = %f, error = %f\n",cfmREl,cfm,relerr)

figure(3);
plot(x,cfREx,'-k',x,cf,'--r','LineWidth',2)
ylabel("Skin-friction coefficient");
xlabel("x")
ylim([0 2])
legend('RK4: \tau_{w}(x)/\rho/2*{U_{\infty}^2}',"0.664/sqrt(Re_{x})",'location', 'northeast')

% Boundary Layer Thickness
for i = 1:numel(y2)
    if y2(i) >= 0.99*Uinf
        fprintf("delta99 = %f\n",eta(i));
        break;
    end
end