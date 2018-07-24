function [tvals,Y,nsteps,lits,h] = solve_DIRK(fcn,Jfcn,tvals,Y0,B,rtol,atol,hmin,hmax,hinit)
% usage: [tvals,Y,nsteps,lits,h] = solve_DIRK(fcn,Jfcn,tvals,Y0,B,rtol,atol,hmin,hmax,hinit)
%
% Adaptive time step diagonally-implicit Runge-Kutta solver for the
% vector-valued ODE problem
%     y' = F(t,Y), t in tvals, y in R^m,
%     Y(t0) = [y1(t0), y2(t0), ..., ym(t0)]'.
%
% Inputs:
%     fcn    = string holding function name for F(t,Y)
%     Jfcn   = string holding function name for Jacobian of F, J(t,Y)
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
%     lits   = number of linear solves required by method
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
c_fails = 0;   % total convergence failures
a_fails = 0;   % total accuracy failures

% set the solver parameters
h_reduce   = 0.1;          % failed step reduction factor
h_safety   = 0.9;          % adaptivity safety factor
h_growth   = 10;           % adaptivity growth bound
ONEMSM     = 1-sqrt(eps);  % coefficients to account for
ONEPSM     = 1+sqrt(eps);  %   floating-point roundoff
ERRTOL     = 1.1;          % upper bound on allowed step error
                           %   (in WRMS norm)

% initialize temporary variables
t = tvals(1);
Ynew = Y0;

% set initial time step size
h = hinit;

% initialize work counters
nsteps = 0;
lits   = 0;

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

      % call stepper routine to take the step and compute error
      % estimate (if applicable)
      if (adaptive)
         if (embedded)
            [Ynew,Yerr,cfail,lin] = DIRKstep_embedded(fcn, Jfcn, Y0, t, h, B);
         else
            [Ynew,Yerr,cfail,lin] = DIRKstep_Richardson(fcn, Jfcn, Y0, t, h, B);
         end
      else
         [Ynew,cfail,lin] = DIRKstep_basic(fcn, Jfcn, Y0, t, h, B);
      end

      % increment linear iteration and number of internal time steps counters
      nsteps = nsteps + 1;
      lits = lits + lin;

      % check for nonlinear convergence/divergence
      if (cfail ~= 0)
         st_fail = 1;
         c_fails = c_fails + 1;
      end

      % if stages succeeded and time step adaptivity enabled, check step accuracy
      if ((st_fail == 0) && adaptive)

         % estimate error in current step
         err_step = max(norm(Yerr./(rtol*Ynew + atol),inf), eps);

         % if error too high, flag step as a failure (will be be recomputed)
         if (err_step > ERRTOL*ONEPSM)
            a_fails = a_fails + 1;
            st_fail = 1;
         end

      end

      % if step was successful (solves succeeded, and error acceptable)
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

      % if step solves or error test failed
      else

         % if already at minimum step, just return with failure
         if (h <= hmin)
            error('Cannot achieve desired accuracy.\n  Consider reducing hmin or increasing rtol.\n');
         end

         % otherwise, reset guess, reduce time step, retry solve
         Ynew = Y0;
         h = h * h_reduce;

      end  % end logic tests for step success/failure

   end  % end while loop attempting to solve steps to next output time

   % store updated solution in output array
   Y(:,tstep) = Ynew;

end  % time step loop

% end solve_DIRK function
end



function [y,yerr,cfail,lits] = DIRKstep_embedded(fcn, Jfcn, y0, t0, h, B)
% usage: [y,yerr,cfail,lits] = DIRKstep_embedded(fcn, Jfcn, y0, t0, h, B)
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
%    cfail = convergence failure flag (0=success; 1=failure)
%    lits  = total linear iterations for step

   % extract ERK method information from B
   [Brows, Bcols] = size(B);
   s = Bcols - 1;        % number of stages
   c = B(1:s,1);         % stage time fraction array
   b = (B(s+1,2:s+1))';  % solution weights (convert to column)
   A = B(1:s,2:s+1);     % RK coefficients
   d = (B(s+2,2:s+1))';  % embedding coefficients

   % initialize storage for RHS vectors, outputs
   k = zeros(length(y0),s);
   lits = 0;
   cfail = 0;

   % set the solver parameters
   newt_maxit = 20;           % max number of Newton iterations
   newt_ftol  = 1e-10;        % Newton solver residual tolerance
   newt_stol  = 1e-10;        % Newton solver solution tolerance

   % set function names for Newton solver residual/Jacobian
   Fun = @F_DIRK;
   Jac = @A_DIRK;

   % set Fdata values for this step
   Fdata.frhs = fcn;    % ODE RHS function name
   Fdata.Jrhs = Jfcn;   % ODE RHS Jacobian function name
   Fdata.B    = B;      % Butcher table
   Fdata.s    = s;      % number of stages
   Fdata.h    = h;      % current step size
   Fdata.yold = y0;     % solution from previous step
   Fdata.t    = t0;     % time of last successful step

   % loop over stages
   for stage = 1:s

      % set Newton initial guess as previous stage solution
      z = y0;

      % set current stage index into Fdata structure
      Fdata.stage = stage;

      % construct RHS comprised of old time data
      %    zi = y_n + h*sum_{j=1}^s (a(i,j)*fj)
      % <=>
      %    zi - h*(a(i,i)*fi) = y_n + h*sum_{j=1}^{i-1} (a(i,j)*fj)
      % =>
      %    rhs = y_n + h*sum_{j=1}^{i-1} (a(i,j)*fj)
      Fdata.rhs = y0;
      for j = 1:stage-1
         Fdata.rhs = Fdata.rhs + h*A(stage,j)*k(:,j);
      end

      % call Newton solver to compute new stage solution
      [z,lin,ierr] = newton(Fun, Jac, z, Fdata, ...
                            newt_ftol, newt_stol, newt_maxit);

      % increment total linear solver statistics
      lits = lits + lin;

      % if Newton method failed, set relevant flags/statistics
      % and break out of stage loop
      if (ierr ~= 0)
         cfail = 1;
         return;
      end

      % construct new stage RHS
      k(:,stage) = fcn(t0+h*c(stage),z);

   end

   % compute new solution and error estimate
   %    ynew = yold + h*sum(b(j)*fj)
   y = y0 + h*k*b;
   yerr = h*k*(b-d);

% end of function
end



function [y,cfail,lits] = DIRKstep_basic(fcn, Jfcn, y0, t0, h, B)
% usage: [y,cfail,lits] = DIRKstep_basic(fcn, Jfcn, y0, t0, h, B)
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
%    cfail = convergence failure flag (0=success; 1=failure)
%    lits  = total linear iterations for step

   % extract ERK method information from B
   [Brows, Bcols] = size(B);
   s = Bcols - 1;        % number of stages
   c = B(1:s,1);         % stage time fraction array
   b = (B(s+1,2:s+1))';  % solution weights (convert to column)
   A = B(1:s,2:s+1);     % RK coefficients

   % initialize storage for RHS vectors, outputs
   k = zeros(length(y0),s);
   lits = 0;
   cfail = 0;

   % set the solver parameters
   newt_maxit = 20;           % max number of Newton iterations
   newt_ftol  = 1e-10;        % Newton solver residual tolerance
   newt_stol  = 1e-10;        % Newton solver solution tolerance

   % set function names for Newton solver residual/Jacobian
   Fun = @F_DIRK;
   Jac = @A_DIRK;

   % set Fdata values for this step
   Fdata.frhs = fcn;    % ODE RHS function name
   Fdata.Jrhs = Jfcn;   % ODE RHS Jacobian function name
   Fdata.B    = B;      % Butcher table
   Fdata.s    = s;      % number of stages
   Fdata.h    = h;      % current step size
   Fdata.yold = y0;     % solution from previous step
   Fdata.t    = t0;     % time of last successful step

   % loop over stages
   for stage = 1:s

      % set Newton initial guess as previous stage solution
      z = y0;

      % set current stage index into Fdata structure
      Fdata.stage = stage;

      % construct RHS comprised of old time data
      %    zi = y_n + h*sum_{j=1}^s (a(i,j)*fj)
      % <=>
      %    zi - h*(a(i,i)*fi) = y_n + h*sum_{j=1}^{i-1} (a(i,j)*fj)
      % =>
      %    rhs = y_n + h*sum_{j=1}^{i-1} (a(i,j)*fj)
      Fdata.rhs = y0;
      for j = 1:stage-1
         Fdata.rhs = Fdata.rhs + h*A(stage,j)*k(:,j);
      end

      % call Newton solver to compute new stage solution
      [z,lin,ierr] = newton(Fun, Jac, z, Fdata, ...
                            newt_ftol, newt_stol, newt_maxit);

      % increment total linear solver statistics
      lits = lits + lin;

      % if Newton method failed, set relevant flags/statistics
      % and break out of stage loop
      if (ierr ~= 0)
         cfail = 1;
         return;
      end

      % construct new stage RHS
      k(:,stage) = fcn(t0+h*c(stage),z);

   end

   % compute new solution and error estimate
   %    ynew = yold + h*sum(b(j)*fj)
   y = y0 + h*k*b;

% end of function
end



function [y,yerr,cfail,lits] = DIRKstep_Richardson(fcn, Jfcn, y0, t0, h, B)
% usage: [y,yerr,cfail,lits] = DIRKstep_Richardson(fcn, Jfcn, y0, t0, h, B)
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
%    cfail = convergence failure flag (0=success; 1=failure)
%    lits  = total linear iterations for step

   % extract ERK method information from B
   [Brows, Bcols] = size(B);
   s = Bcols - 1;        % number of stages
   c = B(1:s,1);         % stage time fraction array
   b = (B(s+1,2:s+1))';  % solution weights (convert to column)
   A = B(1:s,2:s+1);     % RK coefficients

   % initialize storage for RHS vectors, outputs
   k = zeros(length(y0),s);
   lits = 0;
   cfail = 0;

   % set the solver parameters
   newt_maxit = 20;           % max number of Newton iterations
   newt_ftol  = 1e-10;        % Newton solver residual tolerance
   newt_stol  = 1e-10;        % Newton solver solution tolerance

   % set function names for Newton solver residual/Jacobian
   Fun = @F_DIRK;
   Jac = @A_DIRK;

   % set Fdata values for this step
   Fdata.frhs = fcn;    % ODE RHS function name
   Fdata.Jrhs = Jfcn;   % ODE RHS Jacobian function name
   Fdata.B    = B;      % Butcher table
   Fdata.s    = s;      % number of stages
   Fdata.h    = h;      % current step size
   Fdata.yold = y0;     % solution from previous step
   Fdata.t    = t0;     % time of last successful step

   % First compute solution with a single step
   for stage = 1:s

      % set Newton initial guess as previous stage solution
      z = y0;

      % set current stage index into Fdata structure
      Fdata.stage = stage;

      % construct RHS comprised of old time data
      %    zi = y_n + h*sum_{j=1}^s (a(i,j)*fj)
      % <=>
      %    zi - h*(a(i,i)*fi) = y_n + h*sum_{j=1}^{i-1} (a(i,j)*fj)
      % =>
      %    rhs = y_n + h*sum_{j=1}^{i-1} (a(i,j)*fj)
      Fdata.rhs = y0;
      for j = 1:stage-1
         Fdata.rhs = Fdata.rhs + h*A(stage,j)*k(:,j);
      end

      % call Newton solver to compute new stage solution
      [z,lin,ierr] = newton(Fun, Jac, z, Fdata, ...
                            newt_ftol, newt_stol, newt_maxit);

      % increment total linear solver statistics
      lits = lits + lin;

      % if Newton method failed, set relevant flags/statistics
      % and break out of stage loop
      if (ierr ~= 0)
         cfail = 1;
         return;
      end

      % construct new stage RHS
      k(:,stage) = fcn(t0+h*c(stage),z);

   end

   % compute full-step solution
   %    ynew = yold + h*sum(b(j)*fj)
   y1 = y0 + h*k*b;


   % Second compute solution with two half steps
   Fdata.h = h/2;
   for stage = 1:s
      z = y0;   % consider 'smarter' approach for constructing
                % initial guess using results from full-step solution
      Fdata.stage = stage;
      Fdata.rhs = y0;
      for j = 1:stage-1
         Fdata.rhs = Fdata.rhs + h/2*A(stage,j)*k(:,j);
      end
      [z,lin,ierr] = newton(Fun, Jac, z, Fdata, ...
                            newt_ftol, newt_stol, newt_maxit);
      lits = lits + lin;
      if (ierr ~= 0)
         cfail = 1;
         return;
      end
      k(:,stage) = fcn(t0+h/2*c(stage),z);
   end
   y2 = y0 + h/2*k*b;
   Fdata.yold = y2;
   Fdata.t    = t0+h/2;
   for stage = 1:s
      z = y2;
      Fdata.stage = stage;
      Fdata.rhs = y2;
      for j = 1:stage-1
         Fdata.rhs = Fdata.rhs + h/2*A(stage,j)*k(:,j);
      end
      [z,lin,ierr] = newton(Fun, Jac, z, Fdata, ...
                            newt_ftol, newt_stol, newt_maxit);
      lits = lits + lin;
      if (ierr ~= 0)
         cfail = 1;
         return;
      end
      k(:,stage) = fcn(t0+h/2*(1+c(stage)),z);
   end
   y2 = y2 + h/2*k*b;


   % Compute Richardson extrapolant and error estimate
   y = 2*y2-y1;
   yerr = y-y2;

% end of function
end
