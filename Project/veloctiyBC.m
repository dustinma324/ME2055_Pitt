function [u v] = veloctiyBC (u,v,sF,Utop,h,N)

for i = 1:N
u(N,i) = Utop;
u(1,i) = Utop*5;

v(N,i) = 0;
v(1,i) = 0;
end

for j = 1:N
u(j,N) = 0;
u(j,1) = 0;

v(j,N) = 0;
v(j,1) = 0;
end

for i = 2:N-1
    for j = 2:N-1
        u(j,i) = (sF(j+1,i)-sF(j-1,i))/(2*h);
        v(j,i) = (-sF(j,i+1)+sF(j,i-1))/(2*h);
    end
end

end