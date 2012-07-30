function J = J_p2(t, y)
% usage: J = J_p2(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2012
% All Rights Reserved

% model parameters
global Pdata;
a  = Pdata.a; 
b  = Pdata.b; 
ep = Pdata.ep;

% form the ODE Jacobian
J = [-(y(3)+1) + 2*y(1)*y(2),   y(1)*y(1),     -y(1);
        y(3) - 2*y(1)*y(2),    -y(1)*y(1),      y(1);
             -y(3),               0,       -1/ep - y(1)];

% end function
