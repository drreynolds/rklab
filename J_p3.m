function J = J_p3(t, y)
% usage: J = J_p3(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2012
% All Rights Reserved

% extract problem data
global Pdata;
a  = Pdata.a; 
b  = Pdata.b; 
d1 = Pdata.d1;
d2 = Pdata.d2;
m  = Pdata.m;
dx = Pdata.dx;

% extract solution components
u = y(1:m);
v = y(m+1:2*m);

% initialize Jacobian blocks terms
Juu = sparse([],[],[],m,m,3*m);
Juv = sparse([],[],[],m,m,m);
Jvu = sparse([],[],[],m,m,m);
Jvv = sparse([],[],[],m,m,3*m);

% diffusion components
for j=2:m-1
   Juu(j,j-1) = Juu(j,j-1) + d1/dx/dx;
   Juu(j,j)   = Juu(j,j) - 2*d1/dx/dx;
   Juu(j,j+1) = Juu(j,j+1) + d1/dx/dx;
end
for j=2:m-1
   Jvv(j,j-1) = Jvv(j,j-1) + d2/dx/dx;
   Jvv(j,j)   = Jvv(j,j) - 2*d2/dx/dx;
   Jvv(j,j+1) = Jvv(j,j+1) + d2/dx/dx;
end

% reaction components
for j=2:m-1
   Juu(j,j) = Juu(j,j) - (b+1) + 2*u(j)*v(j);
end
for j=2:m-1
   Juv(j,j) = Juv(j,j) + u(j)^2;
end
for j=2:m-1
   Jvv(j,j) = Jvv(j,j) - u(j)^2;
end
for j=2:m-1
   Jvu(j,j) = Jvu(j,j) + b - 2*u(j)*v(j);
end

% combine together into output
J = [Juu, Juv; Jvu, Jvv];

% end function
