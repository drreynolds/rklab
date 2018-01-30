function dy = f_p2(t, y)
% usage: dy = f_p2(t, y)
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

% form the ODE RHS
dy = [a - (y(3)+1)*y(1) + y(1)*y(1)*y(2);
      y(3)*y(1) - y(1)*y(1)*y(2);
      (b-y(3))/ep - y(3)*y(1)];

% end function
