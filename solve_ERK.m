function [tvals,Y,nsteps,h] = solve_ERK(fcn,StabFn,tvals,Y0,B,rtol,atol,hmin,hmax,hinit)
% usage: [tvals,Y,nsteps,h] = solve_ERK(fcn,StabFn,tvals,Y0,B,rtol,atol,hmin,hmax,hinit)
%
% Adaptive time step explicit Runge-Kutta solver for the
% vector-valued ODE problem
%     y' = F(t,Y), t in tvals, y in R^m,
%     Y(t0) = [y1(t0), y2(t0), ..., ym(t0)]'.
%
% Inputs:
%     fcn    = string holding function name for F(t,Y)
%     StabFn = string holding function name for stability constraint on F
%     tvals  = [t0, t1, t2, ..., tN]
%     Y0     = initial value array (column vector of length m)
%     B      = Butcher matrix for IRK coefficients, of the form
%                 B = [c A;
%                      q b;
%                      p b2 ]
%              Here, c is a vector of stage time fractions (s-by-1),
%                    A is a matrix of Butcher coefficients (s-by-s),
%                    q is an integer denoting the method order of accuracy,
%                    b is a vector of solution weights (1-by-s),
%                    p is an integer denoting the embedding order of accuracy,
%                    b2 is a vector of embedding weights (1-by-s),
%              The [p, b2] row is optional.  If that row is not
%              provided the method will default to taking fixed
%              step sizes of size hmin.
%     rtol   = desired relative error of solution  (scalar)
%     atol   = desired absolute error of solution  (vector or scalar)
%     hmin   = minimum internal time step size (hmin <= t(i)-t(i-1), for all i)
%     hmax   = maximum internal time step size (hmax >= hmin)
%     hinit  = initial internal time step size (hmin <= hinit <= hmax)
%
% Outputs:
%     tvals  = the same as the input array tvals
%     y      = [y(t0), y(t1), y(t2), ..., y(tN)], where each
%               y(t*) is a column vector of length m.
%     nsteps = number of internal time steps taken by method
%     h      = last internal step size
%
% Note: to run in fixed-step mode, call with hmin=hmax as the desired
% time step size.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% July 2018
% All Rights Reserved

% determine whether adaptivity is desired
adaptive = 0;
if (abs(hmax-hmin)/abs(hmax) > sqrt(eps))
   adaptive = 1;
end

% if adaptivity enabled, determine approach for error estimation,
% and set the lower-order of accuracy accordingly
[Brows, Bcols] = size(B);
embedded = 0;
p = 0;
if (hmax > hmin)      % check whether adaptivity is desired
   if (Brows > Bcols)
      if (max(abs(B(Bcols+1,2:Bcols))) > eps)   % check for embedding coeffs
         embedded = 1;
         p = B(Bcols+1,1);
      end
   end
end
if (embedded == 0)
   p = B(Bcols,1);
end

% initialize output arrays
N = length(tvals);
m = length(Y0);
Y = zeros(m,N);
Y(:,1) = Y0;

% initialize diagnostics
h_a = 0;       % number of accuracy-limited time steps
h_s = 0;       % number of stability-limited time steps
a_fails = 0;   % total accuracy failures

% set the solver parameters
h_reduce = 0.1;          % failed step reduction factor
h_safety = 0.9;          % adaptivity safety factor
h_growth = 10;           % adaptivity growth bound
h_stable = 0.5;          % fraction of stability step to take
ONEMSM   = 1-sqrt(eps);  % coefficients to account for
ONEPSM   = 1+sqrt(eps);  %   floating-point roundoff
ERRTOL   = 1.1;          % upper bound on allowed step error
                           %   (in WRMS norm)

% initialize temporary variables
t = tvals(1);
Ynew = Y0;

% set initial time step size
h = hinit;

% initialize work counter
nsteps = 0;

% iterate over output time steps
for tstep = 2:length(tvals)

   % loop over internal time steps to get to desired output time
   while (t < tvals(tstep)*ONEMSM)

      % bound internal time step
      h = max([h, hmin]);            % enforce minimum time step size
      h = min([h, hmax]);            % maximum time step size
      h = min([h, tvals(tstep)-t]);  % stop at output time

      % reset step failure flag
      st_fail = 0;

      % compute updated solution and error estimate (if possible)
      if (adaptive)
         if (embedded)
            [Ynew,Yerr] = ERKstep_embedded(fcn, Y0, t, h, B);
         else
            [Ynew,Yerr] = ERKstep_Richardson(fcn, Y0, t, h, B);
         end
      else
         [Ynew] = ERKstep_basic(fcn, Y0, t, h, B);
      end

      % increment number of internal time steps taken
      nsteps = nsteps + 1;

      % if time step adaptivity enabled, check step accuracy
      if (adaptive)

         % estimate error in current step
         err_step = max(norm(Yerr./(rtol*Ynew + atol),inf), eps);

         % if error too high, flag step as a failure (will be recomputed)
         if (err_step > ERRTOL*ONEPSM)
            a_fails = a_fails + 1;
            st_fail = 1;
         end

      end

      % if step was successful (i.e. error acceptable)
      if (st_fail == 0)

         % update solution and time for last successful step
         Y0 = Ynew;
         t  = t + h;

         % for adaptive methods, use error estimate to adapt the time step
         if (adaptive)

            h_old = h;
            if (err_step == 0.0)     % no error, set max possible
               h = tvals(end)-t;
            else                     % set next h (I-controller)
               h = h_safety * h_old * err_step^(-1.0/p);
            end

            % enforce maximum growth rate on step sizes
            h = min(h_growth*h_old, h);

         % otherwise, just use the fixed minimum input step size
         else
            h = hmin;
         end

         % limit time step by explicit stability condition
         hstab = h_stable * StabFn(t, Ynew);

         % keep statistics on how many steps are accuracy vs stability limited
         if (h < hstab)
            h_a = h_a + 1;
         else
            h_s = h_s + 1;
         end
         h = min([h, hstab]);

     % if error test failed
     else

        % if already at minimum step, just return with failure
        if (h <= hmin)
           error('Cannot achieve desired accuracy.\n  Consider reducing hmin or increasing rtol.\n');
           return
        end

        % otherwise, reset guess, reduce time step, retry solve
        Ynew = Y0;
        h    = h * h_reduce;
        h_a  = h_a + 1;

     end  % end logic tests for step success/failure

   end  % end while loop attempting to solve steps to next output time

   % store updated solution in output array
   Y(:,tstep) = Ynew;

end  % time step loop

% end solve_ERK function
end



function [y,yerr] = ERKstep_embedded(fcn, y0, t0, h, B)
% usage: [y,yerr] = ERKstep_embedded(fcn, y0, t0, h, B)
%
% Inputs:
%    fcn = ODE RHS function, f(t,y)
%    y0  = solution at beginning of time step
%    t0  = 'time' at beginning of time step
%    h   = step size to take
%    B   = Butcher table to use
%
% Outputs:
%    y     = new solution at t0+h
%    yerr  = error vector

   % extract ERK method information from B
   [Brows, Bcols] = size(B);
   s = Bcols - 1;        % number of stages
   c = B(1:s,1);         % stage time fraction array
   b = (B(s+1,2:s+1))';  % solution weights (convert to column)
   A = B(1:s,2:s+1);     % RK coefficients
   d = (B(s+2,2:s+1))';  % embedding coefficients

   % initialize storage for RHS vectors
   k = zeros(length(y0),s);

   % loop over stages
   for stage=1:s

      % construct stage solution and evaluate RHS
      %    zi = y_n + h*sum_{j=1}^{i-1} (A(i,j)*f(zj))
      z = y0;
      for j=1:stage-1
         z = z + h*A(stage,j)*k(:,j);
      end

      % construct new stage RHS
      k(:,stage) = fcn(t0+h*c(stage),z);

   end

   % compute new solution and error estimate
   %    ynew = yold + h*sum(b(j)*fj)
   y  = y0 + h*k*b;
   yerr = h*k*(b-d);

% end of function
end



function [y] = ERKstep_basic(fcn, y0, t0, h, B)
% usage: [y] = ERKstep_basic(fcn, y0, t0, h, B)
%
% Inputs:
%    fcn = ODE RHS function, f(t,y)
%    y0  = solution at beginning of time step
%    t0  = 'time' at beginning of time step
%    h   = step size to take
%    B   = Butcher table to use
%
% Outputs:
%    y     = new solution at t0+h

   % extract ERK method information from B
   [Brows, Bcols] = size(B);
   s = Bcols - 1;        % number of stages
   c = B(1:s,1);         % stage time fraction array
   b = (B(s+1,2:s+1))';  % solution weights (convert to column)
   A = B(1:s,2:s+1);     % RK coefficients

   % initialize storage for RHS vectors
   k = zeros(length(y0),s);

   % loop over stages
   for stage=1:s

      % construct stage solution and evaluate RHS
      %    zi = y_n + h*sum_{j=1}^{i-1} (A(i,j)*f(zj))
      z = y0;
      for j=1:stage-1
         z = z + h*A(stage,j)*k(:,j);
      end

      % construct new stage RHS
      k(:,stage) = fcn(t0+h*c(stage),z);

   end

   % compute new solution and error estimate
   %    ynew = yold + h*sum(b(j)*fj)
   y  = y0 + h*k*b;

% end of function
end



function [y,yerr] = ERKstep_Richardson(fcn, y0, t0, h, B)
% usage: [y,yerr] = ERKstep_Richardson(fcn, y0, t0, h, B)
%
% Inputs:
%    fcn = ODE RHS function, f(t,y)
%    y0  = solution at beginning of time step
%    t0  = 'time' at beginning of time step
%    h   = step size to take
%    B   = Butcher table to use
%
% Outputs:
%    y     = new solution at t0+h
%    yerr  = error vector

   % extract ERK method information from B
   [Brows, Bcols] = size(B);
   s = Bcols - 1;        % number of stages
   c = B(1:s,1);         % stage time fraction array
   b = (B(s+1,2:s+1))';  % solution weights (convert to column)
   A = B(1:s,2:s+1);     % RK coefficients

   % initialize storage for RHS vectors
   k = zeros(length(y0),s);

   % First compute solution with a single step
   for stage=1:s

      % construct stage solution and evaluate RHS
      %    zi = y_n + h*sum_{j=1}^{i-1} (A(i,j)*f(zj))
      z = y0;
      for j=1:stage-1
         z = z + h*A(stage,j)*k(:,j);
      end

      % construct new stage RHS
      k(:,stage) = fcn(t0+h*c(stage),z);

   end

   % compute full-step solution
   %    ynew = yold + h*sum(b(j)*fj)
   y1 = y0 + h*k*b;


   % Second compute solution with two half steps
   for stage=1:s

      % construct stage solution and evaluate RHS
      %    zi = y_n + h*sum_{j=1}^{i-1} (A(i,j)*f(zj))
      z = y0;
      for j=1:stage-1
         z = z + h/2*A(stage,j)*k(:,j);
      end

      % construct new stage RHS
      k(:,stage) = fcn(t0+h/2*c(stage),z);

   end
   y2 = y0 + h/2*k*b;
   for stage=1:s

      % construct stage solution and evaluate RHS
      %    zi = y_n + h*sum_{j=1}^{i-1} (A(i,j)*f(zj))
      z = y2;
      for j=1:stage-1
         z = z + h/2*A(stage,j)*k(:,j);
      end

      % construct new stage RHS
      k(:,stage) = fcn(t0+h/2*(1+c(stage)),z);

   end
   y2 = y2 + h/2*k*b;


   % Compute Richardson extrapolant and error estimate
   y = 2*y2-y1;
   yerr = y-y2;

% end of function
end
