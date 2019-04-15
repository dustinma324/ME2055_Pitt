function output = FEuler(f1,f2,f3,y1,y2,y3,h)
g = [f1(y1(end),y2(end),y3(end)),f2(y1(end),y2(end),y3(end)),f3(y1(end),y2(end),y3(end))];

y1(end+1) = y1(end) + h*g(1);
y2(end+1) = y2(end) + h*g(2);
y3(end+1) = y3(end) + h*g(3);

output = [y1' y2' y3'];
end