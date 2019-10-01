function [tvals, Y, nsteps, lits] = solve_ARK(fe,fi,Ji,tvals,Y0,Be,Bi,rtol,atol,hmin,hmax,hinit)
% usage: [tvals, Y, nsteps, lits] = solve_ARK(fe,fi,Ji,tvals,Y0,Be,Bi,rtol,atol,hmin,hmax,hinit)
%
% Adaptive time step additive Runge-Kutta solver for the
% vector-valued ODE problem
%     y' = fe(t,Y) + fi(t,Y), t in tvals, y in R^m,
%     Y(t0) = [y1(t0), y2(t0), ..., ym(t0)]'.
%
% Inputs:
%     fe     = function handle for fe(t,Y)
%     fi     = function handle for fi(t,Y)
%     Ji     = function handle for Jacobian of fi, J(t,Y)
%     tvals  = [t0, t1, t2, ..., tN]
%     Y0     = initial value array (column vector of length m)
%     Be,Bi  = Butcher table matrices for ARK coefficients, of the form
%                 Be = [ce Ae;       Bi = [ci Ai;
%                       q  be;             q  bi;
%                       p  be2 ]           p  bi2 ]
%              Here, ce,ci are vectors of stage time fractions (s-by-1),
%                    Ae,Ai are matrices of Butcher coefficients (s-by-s),
%                    q is an integer denoting the method order of accuracy,
%                    be,bi are vectors of solution weights (1-by-s),
%                    p is an integer denoting the embedding order of accuracy,
%                    be2,bi2 are vectors of embedding weights (1-by-s),
%              The [p, be2] and [p, bi2] rows are optional.  If
%              both of those are not provided the method will default to
%              taking fixed step sizes of size hmin.
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
%
% Note: to run in fixed-step mode, call with hmin=hmax as the desired
% time step size, and set the tolerances to large positive numbers.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% March 2017
% All Rights Reserved

% check for compatible Be,Bi tables
if (size(Be) ~= size(Bi))
   error('solve_ARK error: Be and Bi must have the same size')
end
s = size(Be,2) - 1;          % number of stages
if (Be(s+1,1) ~= Bi(s+1,1))
   error('solve_ARK error: Be and Bi must have the same method order')
end
if (size(Be,1) > size(Be,2))
   if (Be(s+1,2) ~= Bi(s+1,2))
      error('solve_ARK error: Be and Bi must have the same embedding order')
   end
end

% determine whether adaptivity is desired
adaptive = 0;
if (abs(hmax-hmin)/abs(hmax) > sqrt(eps))
   adaptive = 1;
end

% if adaptivity enabled, determine approach for error estimation,
% and set the lower-order of accuracy accordingly
[Brows, Bcols] = size(Be);
embedded = 0;
p = 0;
if (hmax > hmin)      % check whether adaptivity is desired
   if (Brows > Bcols)
      if ( (max(abs(Be(Bcols+1,2:Bcols))) > eps) && ...
           (max(abs(Bi(Bcols+1,2:Bcols))) > eps) )    % check for embedding coeffs
         embedded = 1;
         p = Be(Bcols+1,1);
      end
   end
end
if (embedded == 0)
   p = Be(Bcols,1);
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
   while ((t-tvals(tstep))*h < 0)

      % bound internal time step
      h = max([h, hmin]);            % enforce minimum time step size
      h = min([h, hmax]);            % maximum time step size
      h = min([h, tvals(tstep)-t]);  % stop at output time

      % reset stage failure flag
      st_fail = 0;

      % call stepper routine to take the step and compute error
      % estimate (if applicable); increment internal time steps counter
      if (adaptive)
         if (embedded)
            [Ynew,Yerr,cfail,lin] = ARKstep_embedded(fe, fi, Ji, Y0, t, h, Be, Bi);
            nsteps = nsteps + 1;
         else
            [Ynew,Yerr,cfail,lin] = ARKstep_Richardson(fe, fi, Ji, Y0, t, h, Be, Bi);
            nsteps = nsteps + 3;
         end
      else
         [Ynew,cfail,lin] = ARKstep_basic(fe, fi, Ji, Y0, t, h, Be, Bi);
         nsteps = nsteps + 1;
      end

      % increment linear iteration counter
      lits = lits + lin;

      % check for nonlinear convergence/divergence
      if (cfail ~= 0)
         st_fail = 1;
         c_fails = c_fails + 1;
      end

      % if stages succeeded and time step adaptivity enabled, check step accuracy
      if ((st_fail == 0) & adaptive)

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

% end solve_ARK function
end


%------------------------- Utility routines -------------------------%


function [y,yerr,cfail,lits] = ARKstep_embedded(fe, fi, Ji, y0, t0, h, Be, Bi)
% Inputs:
%    fe = function handle for fe(t,Y)
%    fi = function handle for fi(t,Y)
%    Ji = function handle for Jacobian of fi, J(t,Y)
%    y0 = solution at beginning of time step
%    t0 = 'time' at beginning of time step
%    h  = step size to take
%    Be = explicit Butcher table to use
%    Bi = implicit Butcher table to use
%
% Outputs:
%    y     = new solution at t0+h
%    yerr  = error vector
%    cfail = convergence failure flag (0=success; 1=failure)
%    lits  = total linear iterations for step

   % extract ERK method information from Be
   [Brows, Bcols] = size(Be);
   s = Bcols - 1;          % number of stages
   ce = Be(1:s,1);         % stage time fraction array
   be = (Be(s+1,2:s+1))';  % solution weights (convert to column)
   Ae = Be(1:s,2:s+1);     % RK coefficients
   de = (Be(s+2,2:s+1))';  % embedding coefficients

   % extract DIRK method information from Bi
   [Brows, Bcols] = size(Bi);
   ci = Bi(1:s,1);         % stage time fraction array
   bi = (Bi(s+1,2:s+1))';  % solution weights (convert to column)
   Ai = Bi(1:s,2:s+1);     % RK coefficients
   di = (Be(s+2,2:s+1))';  % embedding coefficients

   % initialize storage for RHS vectors, outputs
   ke = zeros(length(y0),s);
   ki = zeros(length(y0),s);
   lits = 0;
   cfail = 0;

   % set the solver parameters
   newt_maxit = 20;           % max number of Newton iterations
   newt_ftol  = 1e-10;        % Newton solver residual tolerance
   newt_stol  = 1e-10;        % Newton solver solution tolerance

   % set function names for Newton solver residual/Jacobian
   Fun = @F_ARK;
   Jac = @A_ARK;

   % set Fdata values for this step
   Fdata.fe   = fe;    % ODE RHS function names
   Fdata.fi   = fi;
   Fdata.Ji   = Ji;    % ODE RHS Jacobian function name
   Fdata.Be   = Be;    % Butcher tables
   Fdata.Bi   = Bi;
   Fdata.s    = s;     % number of stages
   Fdata.h    = h;     % current step size
   Fdata.yold = y0;    % solution from previous step
   Fdata.t    = t0;    % time of last successful step

   % loop over stages
   for stage = 1:s

      % set Newton initial guess as previous stage solution
      z = y0;

      % set current stage index into Fdata structure
      Fdata.stage = stage;

      % construct RHS comprised of old time data
      %    zi = y_n + h*ai(i,i)*fi(i) + h*sum_{j=1}^{i-1} [ae(i,j)*fe(j) + ai(i,j)*fi(j)]
      % <=>
      %    zi - h*(ai(i,i)*fi) = y_n + h*sum_{j=1}^{i-1} [ae(i,j)*fe(j) + ai(i,j)*fi(j)]
      % =>
      %    rhs = y_n + h*sum_{j=1}^{i-1} [ae(i,j)*fe(j) + ai(i,j)*fi(j)]
      Fdata.rhs = y0;
      for j = 1:stage-1
         Fdata.rhs = Fdata.rhs + h*Ae(stage,j)*ke(:,j) + h*Ai(stage,j)*ki(:,j);
      end

      % call Newton solver to compute new stage solution
      [z,lin,ierr] = newton(Fun, Jac, z, Fdata, newt_ftol, newt_stol, newt_maxit);

      % increment total linear solver statistics
      lits = lits + lin;

      % if Newton method failed, set relevant flags/statistics
      % and break out of stage loop
      if (ierr ~= 0)
         cfail = 1;
         return;
      end

      % construct new stage RHS
      ke(:,stage) = fe(t0+h*ce(stage),z);
      ki(:,stage) = fi(t0+h*ci(stage),z);

   end

   % compute new solution and error estimate
   %    ynew = yold + h*sum(be(j)*fe(j)) + h*sum(bi(j)*fi(j))
   y = y0 + h*ke*be + h*ki*bi;
   yerr = h*ke*(be-de) + h*ki*(bi-di);

% end of function
end





function [y,cfail,lits] = ARKstep_basic(fe, fi, Ji, y0, t0, h, Be, Bi)
% Inputs:
%    fe = function handle for fe(t,Y)
%    fi = function handle for fi(t,Y)
%    Ji = function handle for Jacobian of fi, J(t,Y)
%    y0 = solution at beginning of time step
%    t0 = 'time' at beginning of time step
%    h  = step size to take
%    Be = explicit Butcher table to use
%    Bi = implicit Butcher table to use
%
% Outputs:
%    y     = new solution at t0+h
%    cfail = convergence failure flag (0=success; 1=failure)
%    lits  = total linear iterations for step

   % extract ERK method information from Be
   [Brows, Bcols] = size(Be);
   s = Bcols - 1;          % number of stages
   ce = Be(1:s,1);         % stage time fraction array
   be = (Be(s+1,2:s+1))';  % solution weights (convert to column)
   Ae = Be(1:s,2:s+1);     % RK coefficients
   de = (Be(s+2,2:s+1))';  % embedding coefficients

   % extract DIRK method information from Bi
   [Brows, Bcols] = size(Bi);
   ci = Bi(1:s,1);         % stage time fraction array
   bi = (Bi(s+1,2:s+1))';  % solution weights (convert to column)
   Ai = Bi(1:s,2:s+1);     % RK coefficients
   di = (Be(s+2,2:s+1))';  % embedding coefficients

   % initialize storage for RHS vectors, outputs
   ke = zeros(length(y0),s);
   ki = zeros(length(y0),s);
   lits = 0;
   cfail = 0;

   % set the solver parameters
   newt_maxit = 20;           % max number of Newton iterations
   newt_ftol  = 1e-10;        % Newton solver residual tolerance
   newt_stol  = 1e-10;        % Newton solver solution tolerance

   % set function names for Newton solver residual/Jacobian
   Fun = @F_ARK;
   Jac = @A_ARK;

   % set Fdata values for this step
   Fdata.fe   = fe;    % ODE RHS function names
   Fdata.fi   = fi;
   Fdata.Ji   = Ji;    % ODE RHS Jacobian function name
   Fdata.Be   = Be;    % Butcher tables
   Fdata.Bi   = Bi;
   Fdata.s    = s;     % number of stages
   Fdata.h    = h;     % current step size
   Fdata.yold = y0;    % solution from previous step
   Fdata.t    = t0;    % time of last successful step

   % loop over stages
   for stage = 1:s

      % set Newton initial guess as previous stage solution
      z = y0;

      % set current stage index into Fdata structure
      Fdata.stage = stage;

      % construct RHS comprised of old time data
      %    zi = y_n + h*ai(i,i)*fi(i) + h*sum_{j=1}^{i-1} [ae(i,j)*fe(j) + ai(i,j)*fi(j)]
      % <=>
      %    zi - h*(ai(i,i)*fi) = y_n + h*sum_{j=1}^{i-1} [ae(i,j)*fe(j) + ai(i,j)*fi(j)]
      % =>
      %    rhs = y_n + h*sum_{j=1}^{i-1} [ae(i,j)*fe(j) + ai(i,j)*fi(j)]
      Fdata.rhs = y0;
      for j = 1:stage-1
         Fdata.rhs = Fdata.rhs + h*Ae(stage,j)*ke(:,j) + h*Ai(stage,j)*ki(:,j);
      end

      % call Newton solver to compute new stage solution
      [z,lin,ierr] = newton(Fun, Jac, z, Fdata, newt_ftol, newt_stol, newt_maxit);

      % increment total linear solver statistics
      lits = lits + lin;

      % if Newton method failed, set relevant flags/statistics
      % and break out of stage loop
      if (ierr ~= 0)
         cfail = 1;
         return;
      end

      % construct new stage RHS
      ke(:,stage) = fe(t0+h*ce(stage),z);
      ki(:,stage) = fi(t0+h*ci(stage),z);

   end

   % compute new solution,  ynew = yold + h*sum(be(j)*fe(j)) + h*sum(bi(j)*fi(j))
   y = y0 + h*ke*be + h*ki*bi;

% end of function
end



function [y,yerr,cfail,lits] = ARKstep_Richardson(fe, fi, Ji, y0, t0, h, Be, Bi)
% Inputs:
%    fe = function handle for fe(t,Y)
%    fi = function handle for fi(t,Y)
%    Ji = function handle for Jacobian of fi, J(t,Y)
%    y0 = solution at beginning of time step
%    t0 = 'time' at beginning of time step
%    h  = step size to take
%    Be = explicit Butcher table to use
%    Bi = implicit Butcher table to use
%
% Outputs:
%    y     = new solution at t0+h
%    yerr  = error vector
%    cfail = convergence failure flag (0=success; 1=failure)
%    lits  = total linear iterations for step

   % extract ERK method information from Be
   [Brows, Bcols] = size(Be);
   s = Bcols - 1;          % number of stages
   ce = Be(1:s,1);         % stage time fraction array
   be = (Be(s+1,2:s+1))';  % solution weights (convert to column)
   Ae = Be(1:s,2:s+1);     % RK coefficients
   de = (Be(s+2,2:s+1))';  % embedding coefficients
   p = Be(Bcols,1);

   % extract DIRK method information from Bi
   [Brows, Bcols] = size(Bi);
   ci = Bi(1:s,1);         % stage time fraction array
   bi = (Bi(s+1,2:s+1))';  % solution weights (convert to column)
   Ai = Bi(1:s,2:s+1);     % RK coefficients
   di = (Be(s+2,2:s+1))';  % embedding coefficients

   % initialize storage for RHS vectors, outputs
   ke = zeros(length(y0),s);
   ki = zeros(length(y0),s);
   lits = 0;
   cfail = 0;

   % set the solver parameters
   newt_maxit = 20;           % max number of Newton iterations
   newt_ftol  = 1e-10;        % Newton solver residual tolerance
   newt_stol  = 1e-10;        % Newton solver solution tolerance

   % set function names for Newton solver residual/Jacobian
   Fun = @F_ARK;
   Jac = @A_ARK;

   % set Fdata values for this step
   Fdata.fe   = fe;    % ODE RHS function names
   Fdata.fi   = fi;
   Fdata.Ji   = Ji;    % ODE RHS Jacobian function name
   Fdata.Be   = Be;    % Butcher tables
   Fdata.Bi   = Bi;
   Fdata.s    = s;     % number of stages
   Fdata.h    = h;     % current step size
   Fdata.yold = y0;    % solution from previous step
   Fdata.t    = t0;    % time of last successful step

   % First compute solution with a single step
   for stage = 1:s

      % set Newton initial guess as previous stage solution
      z = y0;

      % set current stage index into Fdata structure
      Fdata.stage = stage;

      % construct RHS comprised of old time data
      %    zi = y_n + h*ai(i,i)*fi(i) + h*sum_{j=1}^{i-1} [ae(i,j)*fe(j) + ai(i,j)*fi(j)]
      % <=>
      %    zi - h*(ai(i,i)*fi) = y_n + h*sum_{j=1}^{i-1} [ae(i,j)*fe(j) + ai(i,j)*fi(j)]
      % =>
      %    rhs = y_n + h*sum_{j=1}^{i-1} [ae(i,j)*fe(j) + ai(i,j)*fi(j)]
      Fdata.rhs = y0;
      for j = 1:stage-1
         Fdata.rhs = Fdata.rhs + h*Ae(stage,j)*ke(:,j) + h*Ai(stage,j)*ki(:,j);
      end

      % call Newton solver to compute new stage solution
      [z,lin,ierr] = newton(Fun, Jac, z, Fdata, newt_ftol, newt_stol, newt_maxit);

      % increment total linear solver statistics
      lits = lits + lin;

      % if Newton method failed, set relevant flags/statistics
      % and break out of stage loop
      if (ierr ~= 0)
         cfail = 1;
         return;
      end

      % construct new stage RHS
      ke(:,stage) = fe(t0+h*ce(stage),z);
      ki(:,stage) = fi(t0+h*ci(stage),z);

   end

   % compute full step solution,  ynew = yold + h*sum(be(j)*fe(j)) + h*sum(bi(j)*fi(j))
   y1 = y0 + h*ke*be + h*ki*bi;


   % Second compute solution with two half steps
   Fdata.h = h/2;
   for stage = 1:s
      z = y0;   % consider 'smarter' approach for constructing
                % initial guess using results from full-step solution
      Fdata.stage = stage;
      Fdata.rhs = y0;
      for j = 1:stage-1
         Fdata.rhs = Fdata.rhs + h/2*Ae(stage,j)*ke(:,j) + h/2*Ai(stage,j)*ki(:,j);
      end
      [z,lin,ierr] = newton(Fun, Jac, z, Fdata, newt_ftol, newt_stol, newt_maxit);
      lits = lits + lin;
      if (ierr ~= 0)
         cfail = 1;
         return;
      end
      ke(:,stage) = fe(t0+h*ce(stage),z);
      ki(:,stage) = fi(t0+h*ci(stage),z);
   end
   y2 = y0 + h/2*ke*be + h/2*ki*bi;
   Fdata.yold = y2;
   Fdata.t    = t0+h/2;
   for stage = 1:s
      z = y2;
      Fdata.stage = stage;
      Fdata.rhs = y0;
      for j = 1:stage-1
         Fdata.rhs = Fdata.rhs + h/2*Ae(stage,j)*ke(:,j) + h/2*Ai(stage,j)*ki(:,j);
      end
      [z,lin,ierr] = newton(Fun, Jac, z, Fdata, newt_ftol, newt_stol, newt_maxit);
      lits = lits + lin;
      if (ierr ~= 0)
         cfail = 1;
         return;
      end
      ke(:,stage) = fe(t0+h*ce(stage),z);
      ki(:,stage) = fi(t0+h*ci(stage),z);
   end
   y2 = y2 + h/2*ke*be + h/2*ki*bi;


   % Compute Richardson extrapolant and error estimate
   y = (2^p)/(2^p-1)*y2 - 1/(2^p-1)*y1;
   yerr = 1/(2^p-1)*(y1-y2);

% end of function
end



function F = F_ARK(z, Fdata)
% Inputs:  z = current guess for stage solution
%          Fdata = structure containing extra information for evaluating F.
% Outputs: F = residual at current guess
%
% This function computes the (non)linear residuals for an intermediate
% stage solution, through calling the user-supplied (in Fdata) ODE
% right-hand side function.

   % extract ARK method information from Fdata
   Bi = Fdata.Bi;
   [Brows, Bcols] = size(Bi);
   s  = Bcols - 1;
   ci = Bi(1:s,1);
   Ai = Bi(1:s,2:s+1);
   h  = Fdata.h;
   st = Fdata.stage;
   t  = Fdata.t + Fdata.h*ci(st);

   % form the ARK residual
   %    F = z - rhs - h*(ai(stage,stage)*fstage)
   F = z - Fdata.rhs - h*Ai(st,st)*Fdata.fi(t, z);

% end of function
end



function Amat = A_ARK(z, Fdata)
% Inputs:  z = current guess for stage solution
%          Fdata = structure containing extra information for evaluating F.
% Outputs: Amat = Jacobian at current guess
%
% This function computes the Jacobian of each intermediate stage residual
% for a multi-stage ARK method, through calling the user-supplied (in
% Fdata) ODE Jacobian function.

   % extract ARK method information from Fdata
   Bi = Fdata.Bi;
   [Brows, Bcols] = size(Bi);
   s   = Bcols - 1;
   ci  = Bi(1:s,1);
   bi  = (Bi(s+1,2:s+1))';
   Ai  = Bi(1:s,2:s+1);
   st  = Fdata.stage;
   t   = Fdata.t + Fdata.h*ci(st);

   % form the ARK Jacobian
   Amat = eye(length(z)) - Fdata.h*Ai(st,st)*Fdata.Ji(t, z);

% end of function
end
