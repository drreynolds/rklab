function stab_region(A,b,box,fig,fmt)
% Usage: stab_region(A,b,box,fig,fmt)
%
% Inputs:
%    A is a Butcher table matrix
%    b is a Butcher table gluing coefficients
%    box = [xl, xr, yl, yr] is the bounding box for the sub-region
%          of the complex plane in which to perform the test
%    fig = figure handle to use
%    fmt = plot format string to use
%
% We consider the RK stability function
%    R(eta) = 1 + eta*b'*((I-eta*A)\e)
%
% We sample the values in 'box' within the complex plane, plugging
% each value into |R(eta)|, and plot the contour of this function
% having value 1.
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2016, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

% extract the components of the Butcher table
s = length(b);
I = eye(s);
e = ones(s,1);

% ensure that b is a row vector
b = reshape(b,1,s);

% set mesh of sample points
xl = box(1);
xr = box(2);
yl = box(3);
yr = box(4);
N = 1000;
x = linspace(xl,xr,N);
y = linspace(yl,yr,N);

% evaluate |R(eta)| for each eta in the mesh
for j=1:N
  for i=1:N
    eta = x(i) + y(j)*sqrt(-1);
    R(j,i) = abs(1 + eta*b*((I - eta*A)\e));
  end
end

% create contour
%c = contourc(x,y,R,[1 1]);
c = contourc(x,y,R,[1+eps 1+eps]);

% assemble plot
figure(fig);
hold on
idx = 1;
while idx < size(c,2)
  cols = idx+1:idx+c(2,idx);
  X = c(1,cols);
  Y = c(2,cols);
  plot(X, Y, fmt)
  idx = idx + c(2,idx) + 1;
end
hold off

% end of function
