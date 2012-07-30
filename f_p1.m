function dy = f_p1(t, y)
% usage: dy = f_p1(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2012
% All Rights Reserved

% model parameters
global Pdata;
ep = Pdata.ep;

% extract variables
u = y(1);
v = y(2);

% form the ODE RHS
du = v;
dv = (v - v*u^2)/ep - u;
dy = [du; dv];

% end function
