function [sF w] = boundaryCondition(sF,w,N,h,Utop)
sF(:,1)    = 0;  % West
sF(:,N)    = 0;  % East
sF(1,:)    = 0;  % South
sF(N,:)    = 0;  % North

for j = 1:N
    w(j,N) = -2*sF(j,N-1)/h^2;            %east
    w(j,1) = -2*sF(j,2)/h^2;              %west
end

for i = 1:N
    w(N,i) = -2*(sF(N-1,i)/h^2)-2*Utop/h;   %north
    w(1,i) = -2*sF(2,i)/h^2   +2*Utop*5/h;                %south
end

end