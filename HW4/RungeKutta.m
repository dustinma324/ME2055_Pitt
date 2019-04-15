function output = RungeKutta(f1,f2,f3,y1,y2,y3,h)
% Runge-Kutta Constants
k1 = [f1(y1(end),y2(end),y3(end)),  f2(y1(end),y2(end),y3(end)), f3(y1(end),y2(end),y3(end))];
k2 = [f1(y1(end)+k1(1)*h/2,y2(end)+k1(2)*h/2,y3(end)+k1(3)*h/2),    f2(y1(end)+k1(1)*h/2,y2(end)+k1(2)*h/2,y3(end)+k1(3)*h/2), f3(y1(end)+k1(1)*h/2, y2(end)+k1(2)*h/2, y3(end)+k1(3)*h/2)];
k3 = [f1(y1(end)+k2(2)*h/2,y2(end)+k2(2)*h/2,y3(end)+k2(3)*h/2),    f2(y1(end)+k2(1)*h/2,y2(end)+k2(2)*h/2,y3(end)+k2(3)*h/2), f3(y1(end)+k2(1)*h/2, y2(end)+k2(2)*h/2, y3(end)+k2(3)*h/2)];
k4 = [f1(y1(end)+k3(1)*h,y2(end)+k3(2)*h,y3(end)+k3(3)*h),    f2(y1(end)+k3(1)*h,y2(end)+k3(2)*h,y3(end)+k3(3)*h), f3(y1(end)+k3(1)*h,y2(end)+k3(2)*h,y3(end)+k3(3)*h)];

% Solving 1st order ODE functions
y1(end+1) = y1(end) + h/6*(k1(1)+2*k2(1)+2*k3(1)+k4(1));
y2(end+1) = y2(end) + h/6*(k1(2)+2*k2(2)+2*k3(2)+k4(2));
y3(end+1) = y3(end) + h/6*(k1(3)+2*k2(3)+2*k3(3)+k4(3));

output = [y1' y2' y3'];
end
