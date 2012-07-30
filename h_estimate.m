function h = h_estimate(Y1, Y2, h_old, r_tol, a_tol, p, reset)
% Usage: h = h_estimate(Y1, Y2, h_old, r_tol, a_tol, p, reset)
%
% Adaptive time step estimation routine, that attempts to guarantee a 
% local truncation error satisfying the bound
%    norm(local error,inf) < norm(rtol*y_i + atol_i,inf)
%
% Inputs:  Y1 -- first time-evolved solution (more accurate)
%          Y2 -- second time-evolved solution (same size as Y1)
%       h_old -- previous time step size (scalar)
%       r_tol -- desired relative tolerance (scalar)
%       a_tol -- desired absolute tolerance (vector, size of Y1)
%           p -- order of accuracy for method
%       reset -- flag to denote reset of history
%
% Output:   h -- new time step
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% November 2012
% All Rights Reserved

if (reset == 1)
   h = 1;
   return;
end

% set variables
safety = 0.9;
dt_growth = 10;
alpha = 1.0/p;

% compute new error estimate ratio, bound from below 
Eratio = max(norm((Y1 - Y2)./(r_tol*Y1 + a_tol),inf), eps);

% if Eratio == 0, set a maximum h value
if (Eratio == 0.0)
   h = 1000;
   return;
end

% compute updated time step
h = safety * h_old * Eratio^(-alpha);

% enforce growth factor
h = min(dt_growth*h_old, h);


% end of function