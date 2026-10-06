function [P,nloc] = recmesh2412(a,b,shift,nElems)
% makes a rectangular mesh with four node plate elements Q4 (2422)  
n = nElems.x;
m = nElems.y;
%

x = shift.x;
y = shift.y;
% Geospace for spacing
xk=geospace(0+x,a+x,n+1,1)';
yk=geospace(0+y,b+y,m+1,1)';
% generate nodal coordinates in xy-plane
P = [];
for k=1:m+1
    P = [P; xk yk(k).*ones(n+1,1)];  
end
% generate elememnt connectivity
nloc = [];
for l = 1:1:m
  for k = 1:1:n
    loc = [k k+1 k+2+n n+k+1] + (l-1)*(n+1);
    nloc = [nloc; loc];
  end
end

