%% Foward Euler Method
clear all; close all; clc;
values = [0.5 0.2 0.09 0.07 0.05];
for change = 1:numel(values)
values = [0.5 0.2 0.09 0.07 0.05];
% Change These Values
deltaEta = values(change);
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
fprintf("FE iter = %d\n",iter)

%plot(eta,y1,'--k',eta,y2,'-r',eta,y3,'-.b','LineWidth',2)
figure(1);
limit = zeros(1,numel(eta))+1;
if change == 1
    plot(eta,limit,'--k',eta,y2,'LineWidth',1.5)
else
    plot(eta,y2,'LineWidth',1.5)
end
hold on
ylim([0 1.2]);
xlabel('\eta', 'fontsize', 16);
ylabel("f'=u/U_{\infty}",'fontsize',16);

for i = 1:numel(y2)
   err(i) = (1-y2(i))^2;
end
figure(2);
loglog(eta,err,'LineWidth',1.75);
hold on
title("FE optimum \Delta \eta")
xlabel('log(\eta)')
ylabel('log(error)')
clear all;
end
figure(1);legend("f^{'}(\infty)",'\Delta\eta=0.5','\Delta\eta=0.2','\Delta\eta=0.09','\Delta\eta=0.07','\Delta\eta=0.05','location','southeast')
figure(2);legend('\Delta\eta=0.5','\Delta\eta=0.2','\Delta\eta=0.09','\Delta\eta=0.07','\Delta\eta=0.05','location','southwest')