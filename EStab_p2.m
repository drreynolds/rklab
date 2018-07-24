function dt = EStab_p2(t, y)
% usage: dt = EStab_p2(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2012
% All Rights Reserved

% compute the Jacobian ODE RHS
global Pdata;
J = [-(y(3)+1) + 2*y(1)*y(2),   y(1)*y(1),     -y(1);
     y(3) - 2*y(1)*y(2),    -y(1)*y(1),      y(1);
     -y(3),               0,       -1/Pdata.ep - y(1)];

% determine the largest eigenvalue magnitude
lam = max(abs(eig(J)));

% assume explicit stability region includes Euler stability region
% (this assumes that the eigenvalue is in fact negative).
dt = 1/lam;

% end function
