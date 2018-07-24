function [tvals,Y,ns,nf] = solve_RMIS_MIS_38(fs,ff,tvals,Y0,rtol,atol,hmin,hmax,hinit,hfinit)
% usage: [tvals,Y,ns,nf] = solve_RMIS_MIS_38(fs,ff,tvals,Y0,rtol,atol,hmin,hmax,hinit,hfinit)
%
% Adaptive time step RMIS-3/8 and MIS-3/8, explicit+explicit
% multirate Runge-Kutta method for the vector-valued ODE problem
%     y' = fs(t,Y) + ff(t,Y), t >= t0, y in R^n,
%     Y0 = [y1(t0), y2(t0), ..., yn(t0)]'.
% The individual time steps are performed using the step_RMIS_MIS.m
% function; this routine sets the inner and outer Butcher tables to
% the "3/8-Rule" table, and calls step_RMIS_MIS.m in a loop to fill
% the output arrays.
%
% We perform step size adaptivity for the outer method using the
% RMIS+MIS embedding; we perfom step size adaptivity for the inner
% method using it's built-in embedding (if possible; otherwise we
% use Richardson extrapolation).
%
% Inputs:
%     fs     = function handle for (slow) ODE RHS
%     ff     = function handle for (fast) ODE RHS
%     tvals  = array of desired output times, [t0, t1, t2, ..., tN]
%     Y0     = solution vector at start of step (column vector of length n)
%     rtol   = desired relative error of solution  (scalar)
%     atol   = desired absolute error of solution  (vector or scalar)
%     hmin   = minimum internal time step size (hmin <= t(i)-t(i-1), for all i)
%     hmax   = maximum internal time step size (hmax >= hmin)
%     hinit  = initial slow time step size (hmin <= hinit <= hmax)
%     hfinit = initial fast time step size (hfinit < hinit)
%
% Outputs:
%     tvals  = the same as the input array tvals
%     Y      = [y(t0), y(t1), y(t2), ..., y(tN)], where each
%               y(t*) is a column vector of length n.
%     ns     = number of 'slow' time steps taken by method
%     nf     = number of 'fast' time steps taken by method
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% July 2018
% All Rights Reserved

% set 3/8-Rule Butcher table   Bo = [co Ao; qo bo ]
B = [  0    0    0    0    0;
      1/3  1/3   0    0    0;
      2/3 -1/3   1    0    0;
       1    1   -1    1    0;
       4   1/8  3/8  3/8  1/8];
p = 3;  % embedding order of accuracy

% initialize output arrays
N = length(tvals)-1;
n = length(Y0);
Y = zeros(n,N+1);
Y(:,1) = Y0;

% initialize diagnostics
ns = 0;
nf = 0;

% set the solver parameters
h_reduce = 0.1;          % failed step reduction factor
h_safety = 0.9;          % adaptivity safety factor
h_growth = 10;           % adaptivity growth bound
ONEMSM   = 1-sqrt(eps);  % coefficients to account for
ONEPSM   = 1+sqrt(eps);  %   floating-point roundoff
ERRTOL   = 1.1;          % upper bound on allowed step error
                           %   (in WRMS norm)
% initialize temporary variables
t = tvals(1);
Jf = @(t,y) 0;  % no Jacobian required for explicit+explicit methods

% set initial time step sizes
h = hinit;
hf = hfinit;

% iterate over output time steps
for tstep = 2:length(tvals)

   % loop over internal time steps to get to desired output time
   while (t < tvals(tstep+1)*ONEMSM)

      % bound internal time step
      h = max([h, hmin]);            % enforce minimum time step size
      h = min([h, hmax]);            % maximum time step size
      h = min([h, tvals(tstep)-t]);  % stop at output time
      Fdata.h = h;
      Fdata.yold = Y0;

      % reset step failure flag
      st_fail = 0;

      % call RMIS stepper to do the work, increment counters
      [Ynew,Yerr,m,hf] = step_RMIS_MIS(fs,ff,Jf,t,Y0,B,B,h,hf,rtol,atol);
      ns = ns + 1;
      nf = nf + m;

      % estimate error in current step
      err_step = max(norm(Yerr./(rtol*abs(Ynew) + atol),inf), eps);

      % if error too high, flag step as a failure (will be recomputed)
      if (err_step > ERRTOL*ONEPSM)
         st_fail = 1;
      end

      % if step was successful (i.e. error acceptable)
      if (st_fail == 0)

         % update solution and time for last successful step
         Y0 = Ynew;
         t  = t + h;

         % use error estimate to adapt the time step
         h_old = h;
         if (err_step == 0.0)     % no error, set max possible
            h = tvals(end)-t;
         else                     % set next h (I-controller)
            h = h_safety * h_old * err_step^(-1.0/p);
         end

         % enforce maximum growth rate on step sizes
         h = min(h_growth*h_old, h);

     % if error test failed
     else

        % if already at minimum step, just return with failure
        if (h <= hmin)
           error('Cannot achieve desired accuracy.\n  Consider reducing hmin or increasing rtol.\n');
        end

        % otherwise, reset guess, reduce time step, retry solve
        h = h * h_reduce;

     end  % end logic tests for step success/failure

   end  % end while loop attempting to solve steps to next output time

   % store updated solution in output array
   Y(:,tstep+1) = Ynew;

end  % time output loop

% end function
