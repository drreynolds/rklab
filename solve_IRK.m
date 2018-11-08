function [tvals,Y,nsteps,lits,h] = solve_IRK(fcn,Jfcn,tvals,Y0,B,rtol,atol,hmin,hmax,hinit)
% usage: [tvals,Y,nsteps,lits,h] = solve_IRK(fcn,Jfcn,tvals,Y0,B,rtol,atol,hmin,hmax,hinit)
%
% Fixed time step implicit Runge-Kutta solver for the vector-valued
% ODE problem
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
%                      p b ]
%              Here, c is a vector of stage time fractions (s-by-1),
%                    A is a matrix of Butcher coefficients (s-by-s),
%                    p is an integer denoting the method order of accuracy,
%                    b is a vector of solution weights (1-by-s),
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
%     nsteps = number of internal time steps taken by method.
%              Note: if adaptivity is enabled, this includes the
%              two half steps used in the Richardson error estimate
%              and extrapolation.
%     lits   = number of linear solves required by method
%     h      = last internal step size
%
% Note1: to run in fixed-step mode, call with hmin=hmax as the desired
% time step size, and set the tolerances to large positive numbers.
%
% Note2: if adaptivity is requested, we use Richardson
% extrapolation.  Specifically, for every time step, we
% solve using both a single step of size h and two steps of size
% h/2.  The difference of these provides an estimate on the local
% error,  Moreover, an appropriately-chosen linear combination of
% these provides a solution that is one order of accuracy higher
% than the method itself -- this is the solution that is stored and
% returned to the user.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2012
% All Rights Reserved

% get number of stages and internal time fractions for IRK method
[Brows, Bcols] = size(B);
s = Bcols - 1;
c = B(1:s,1);

% check whether time step adaptivity is desired
adaptive = 0;
if (abs(hmax-hmin)/abs(hmax) > sqrt(eps))
   p = B(s+1,1);     % order of accuracy for method
   c1 = -1/(2^p-1);  % Richardson extrapolation factor for h step
   c2 = 1 - c1;      % Richardson extrapolation factor for h/2 step
   adaptive = 1;
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
h_reduce = 0.1;          % failed step reduction factor
h_safety = 0.9;          % adaptivity safety factor
h_growth = 10;           % adaptivity growth bound
ONEMSM   = 1-sqrt(eps);  % coefficients to account for
ONEPSM   = 1+sqrt(eps);  %   floating-point roundoff
ERRTOL   = 1.1;          % upper bound on allowed step error
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

      % reset solve failure flag
      st_fail = 0;

      % call stepper routine to take the step and compute error
      % estimate (if applicable); increment internal time steps counter
      if (adaptive)
         [Ynew,Yerr,cfail,lin] = IRKstep_Richardson(fcn, Jfcn, Y0, t, h, B);
         nsteps = nsteps + 3;
      else
         [Ynew,cfail,lin] = IRKstep_basic(fcn, Jfcn, Y0, t, h, B);
         nsteps = nsteps + 1;
      end

      % increment linear iteration counter
      lits = lits + lin;

      % check for nonlinear convergence/divergence
      if (cfail ~= 0)
         st_fail = 1;
         c_fails = c_fails + 1;
      end

      % if solves succeeded and time step adaptivity enabled, check step accuracy
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
            error('Cannot achieve desired accuracy.\n Consider reducing hmin or increasing rtol.\n');
         end

         % otherwise, reset guess, reduce time step, retry solve
         Ynew = Y0;
         h = h * h_reduce;

      end  % end logic tests for step success/failure

   end  % end while loop attempting to solve steps to next output time

   % store updated solution in output array
   Y(:,tstep) = Ynew;

end  % time step loop

% end solve_IRK function
end


%------------------------- Utility routines -------------------------%


function [y,cfail,lits] = IRKstep_basic(fcn, Jfcn, y0, t0, h, B)
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

   % extract IRK method information from B
   [Brows, Bcols] = size(B);
   s = Bcols - 1;        % number of stages
   c = B(1:s,1);         % stage time fraction array
   b = (B(s+1,2:s+1))';  % solution weights (convert to column)
   A = B(1:s,2:s+1);     % RK coefficients

   % initialize outputs
   m = length(y0);
   lits = 0;
   cfail = 0;

   % set the solver parameters
   newt_maxit = 20;           % max number of Newton iterations
   newt_ftol  = 1e-10;        % Newton solver residual tolerance
   newt_stol  = 1e-10;        % Newton solver solution tolerance

   % set function names for Newton solver residual/Jacobian
   Fun = @F_IRK;
   Jac = @A_IRK;

   % set Fdata values for this step
   Fdata.frhs = fcn;    % ODE RHS function name
   Fdata.Jrhs = Jfcn;   % ODE RHS Jacobian function name
   Fdata.B    = B;      % Butcher table
   Fdata.s    = s;      % number of stages
   Fdata.h    = h;      % current step size
   Fdata.yold = y0;     % solution from previous step
   Fdata.t    = t0;     % time of last successful step


   % solve with time step h

   % set Newton initial guesses as previous step solution
   z = zeros(s*m,1);
   for i = 0:s-1
      z(i*m+1:(i+1)*m) = y0;
   end

   % call Newton solver to update solution in time
   [z,lin,ierr] = newton(Fun, Jac, z, Fdata, newt_ftol, newt_stol, newt_maxit);

   % increment total linear solver statistics
   lits = lits + lin;

   % if Newton method failed, set relevant flags/statistics
   if (ierr ~= 0)
      cfail = 1;
      return;
   end

   % compute solution with this h
   y = Y_IRK(z,Fdata);

% end of function
end



function [y,yerr,cfail,lits] = IRKstep_Richardson(fcn, Jfcn, y0, t0, h, B)
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

   % extract IRK method information from B
   [Brows, Bcols] = size(B);
   s = Bcols - 1;        % number of stages
   c = B(1:s,1);         % stage time fraction array
   b = (B(s+1,2:s+1))';  % solution weights (convert to column)
   A = B(1:s,2:s+1);     % RK coefficients
   p = B(Bcols,1);

   % initialize outputs
   m = length(y0);
   lits = 0;
   cfail = 0;

   % set the solver parameters
   newt_maxit = 20;           % max number of Newton iterations
   newt_ftol  = 1e-10;        % Newton solver residual tolerance
   newt_stol  = 1e-10;        % Newton solver solution tolerance

   % set function names for Newton solver residual/Jacobian
   Fun = @F_IRK;
   Jac = @A_IRK;

   % set Fdata values for this step
   Fdata.frhs = fcn;    % ODE RHS function name
   Fdata.Jrhs = Jfcn;   % ODE RHS Jacobian function name
   Fdata.B    = B;      % Butcher table
   Fdata.s    = s;      % number of stages
   Fdata.h    = h;      % current step size
   Fdata.yold = y0;     % solution from previous step
   Fdata.t    = t0;     % time of last successful step


   % First compute solution with a single step

   % set Newton initial guesses as previous step solution
   z = zeros(s*m,1);
   for i = 0:s-1
      z(i*m+1:(i+1)*m) = y0;
   end

   % call Newton solver to update solution in time
   [z,lin,ierr] = newton(Fun, Jac, z, Fdata, newt_ftol, newt_stol, newt_maxit);

   % increment total linear solver statistics
   lits = lits + lin;

   % if Newton method failed, set relevant flags/statistics
   if (ierr ~= 0)
      cfail = 1;
      return;
   end

   % compute solution with this h
   y1 = Y_IRK(z,Fdata);


   % Second compute solution with two half steps
   Fdata.h    = h/2;
   Fdata.yold = y0;
   Fdata.t    = t0;
   z = zeros(s*m,1);
   for i = 0:s-1
      ti = c(i+1)*0.5;
      z(i*m+1:(i+1)*m) = (1-ti)*y0 + ti*y1;
   end
   [z,lin,ierr] = newton(Fun, Jac, z, Fdata, newt_ftol, newt_stol, newt_maxit);
   lits = lits + lin;
   if (ierr ~= 0)
      cfail = 1;
      return;
   end
   y2 = Y_IRK(z,Fdata);
   Fdata.yold = y2;
   Fdata.t    = t0+h/2;
   z = zeros(s*m,1);
   for i = 0:s-1
      ti = 0.5 + c(i+1)*0.5;
      z(i*m+1:(i+1)*m) = (1-ti)*y0 + ti*y1;
   end
   [z,lin,ierr] = newton(Fun, Jac, z, Fdata, newt_ftol, newt_stol, newt_maxit);
   lits = lits + lin;
   if (ierr ~= 0)
      cfail = 1;
      return;
   end
   y2 = Y_IRK(z,Fdata);

   % Compute Richardson extrapolant and error estimate
   y = (2^p)/(2^p-1)*y2 - 1/(2^p-1)*y1;
   yerr = 1/(2^p-1)*(y1-y2);

% end of function
end



function y = Y_IRK(z, Fdata)
% Inputs:
%    z     = stage solutions [z1, ..., zs]
%    Fdata = structure containing extra problem information
%
% Outputs:
%    y     = step solution built from the z values

   % extract method information from Fdata
   B = Fdata.B;
   [Brows, Bcols] = size(B);
   s = Bcols - 1;
   c = B(1:s,1);
   b = (B(s+1,2:s+1))';
   A = B(1:s,2:s+1);

   % get some problem information
   zlen = length(z);
   nvar = floor(zlen/s);
   if (nvar*s ~= zlen)
      error('Y_IRK error: input has incorrect length (must be a multiple of s)');
   end

   % reshape our z arguments into separate vectors for each stage
   z = reshape(z,nvar,s);

   % call f at our stages
   f = zeros(nvar,s);
   for is=1:s
      t = Fdata.t + Fdata.h*c(is);
      f(:,is) = Fdata.frhs(t, z(:,is));
   end

   % form the solution
   %    ynew = yold + h*sum(b(j)*fj)
   y = Fdata.yold + Fdata.h*f*b;

   % end of function
end



function F = F_IRK(z, Fdata)
% Inputs:  z = current guesses for [z1, ..., zs]
%          Fdata = structure containing extra information for evaluating F.
% Outputs: F = residual at current guess
%
% This function computes the (non)linear residuals for each intermediate
% stage solution, through calling the user-supplied (in Fdata) ODE
% right-hand side function.

   % extract IRK method information from Fdata
   B = Fdata.B;
   [Brows, Bcols] = size(B);
   s = Bcols - 1;
   c = B(1:s,1);
   b = (B(s+1,2:s+1))';
   A = B(1:s,2:s+1);

   % get some problem information
   zlen = length(z);
   nvar = floor(zlen/s);
   if (nvar*s ~= zlen)
      error('F_IRK error: input argument has incorrect length (must be a multiple of s)');
   end

   % reshape our z arguments
   z = reshape(z,nvar,s);

   % call f at our guesses
   f = zeros(nvar,s);
   for is=1:s
      t = Fdata.t + Fdata.h*c(is);
      f(:,is) = Fdata.frhs(t, z(:,is));
   end

   % form the IRK residuals
   %    Fs = zs - y_n - h*sum(a(s,j)*fj)
   F = zeros(nvar,s);
   for is=1:s
      F(:,is) = z(:,is) - Fdata.yold;
      for j=1:s
         F(:,is) = F(:,is) - Fdata.h*A(is,j)*f(:,j);
      end
   end

   % reshape our output
   F = reshape(F, nvar*s, 1);

% end of function
end



function Amat = A_IRK(z, Fdata)
% Inputs:  z = current guesses for [z1, ..., zs]
%          Fdata = structure containing extra information for evaluating F.
% Outputs: Amat = Jacobian at current guess
%
% This function computes the Jacobian of each intermediate stage residual
% for a multi-stage IRK method, through calling the user-supplied (in Fdata)
% ODE Jacobian function.

   % extract IRK method information from Fdata
   B = Fdata.B;
   [Brows, Bcols] = size(B);
   s = Bcols - 1;
   c = B(1:s,1);
   b = (B(s+1,2:s+1))';
   A = B(1:s,2:s+1);

   % get some problem information
   zlen = length(z);
   nvar = floor(zlen/s);
   if (nvar*s ~= zlen)
      error('A_IRK error: input argument has incorrect length (must be a multiple of s)');
   end

   % reshape our z arguments
   z = reshape(z,nvar,s);

   % call J at each of our guesses
   J = zeros(nvar,nvar,s);
   for is=1:s
      t = Fdata.t + Fdata.h*c(is);
      J(:,:,is) = Fdata.Jrhs(t, z(:,is));
   end

   % form the IRK Jacobian
   Amat = zeros(nvar*s);
   for j=1:s
      for i=1:s
         Amat(nvar*(i-1)+1:nvar*i,nvar*(j-1)+1:nvar*j) = A(i,j)*J(:,:,j);
      end
   end
   Amat = eye(nvar*s) - Fdata.h*Amat;

% end of function
end
