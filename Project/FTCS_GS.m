function [sF w err] = FTCS_GS(sF,w,h,dt,N,nu,Utop)
% Advancing vorticity in time
f = zeros(N,N);
for i = 2:N-1
    for j = 2:N-1
        f(j,i) = -((sF(j+1,i)-sF(j-1,i))*(w(j,i+1)-w(j,i-1)))/(4*h*h)...
            +((sF(j,i+1)-sF(j,i-1))*(w(j+1,i)-w(j-1,i)))/(4*h*h)...
            +nu*(w(j,i+1)+w(j,i-1)+w(j+1,i)+w(j-1,i)-4*w(j,i))/(h*h);
        
        w(j,i) = dt*f(j,i) + w(j,i);
    end
end

[sF w] = boundaryCondition(sF,w,N,h,Utop);
before = norm(sF,2);

for i = 2:N-1
    for j = 2:N-1
        sF(j,i) = 0.25*(h*h*w(j,i)+sF(j,i+1)+sF(j,i-1)+sF(j+1,i)+sF(j-1,i));
    end
end

after = norm(sF,2);
err = abs(before - after);

end