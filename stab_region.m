function [X,Y] = stab_region(A,b,box)
% Usage: [X,Y] = stab_region(A,b,box)
% 
% Inputs:
%    A is a Butcher table matrix
%    b is a Butcher table gluing coefficients
%    box = [xl, xr, yl, yr] is the bounding box for the sub-region
%          of the complex plane in which to perform the test
%
% Outputs:
%    X is an array of real components of the stability boundary
%    Y is an array of imaginary components of the stability boundary
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
X = c(1,2:c(2,1)+1);
Y = c(2,2:c(2,1)+1);
  
% end of function
