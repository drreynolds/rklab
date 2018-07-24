function dt = EStab_p1(t, y)
% usage: dt = EStab_p1(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2012
% All Rights Reserved

% compute the Jacobian ODE RHS
global Pdata;
J = [0, 1; -1, (1-y(1)^2)/Pdata.ep];

% determine the largest eigenvalue magnitude
lam = max(abs(eig(J)));

% assume explicit stability region includes Euler stability region
% (this assumes that the eigenvalues are in fact negative).
dt = 1/lam;

% end function
