function dydt = blausius(t,x)
dydt = zeros(3,1);

dydt(1) = x(2);
dydt(2) = x(3);
dydt(3) = -0.5*x(1)*x(3);
end
