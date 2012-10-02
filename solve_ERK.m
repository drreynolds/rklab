function [tvals,Y,nsteps] = solve_ERK(fcn,StabFn,tvals,Y0,B,rtol,atol,hmin,hmax)
% usage: [tvals,Y,nsteps] = solve_ERK(fcn,StabFn,tvals,Y0,B,rtol,atol,hmin,hmax)
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
%
% Outputs: 
%     tvals  = the same as the input array tvals
%     y      = [y(t0), y(t1), y(t2), ..., y(tN)], where each
%               y(t*) is a column vector of length m.
%     nsteps = number of internal time steps taken by method
%
% Note: to run in fixed-step mode, call with hmin=hmax as the desired 
% time step size, and set the tolerances to large positive numbers.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2012
% All Rights Reserved

   
% extract ERK method information from B
[Brows, Bcols] = size(B);
s = Bcols - 1;        % number of stages
c = B(1:s,1);         % stage time fraction array
b = (B(s+1,2:s+1))';  % solution weights (convert to column)
A = B(1:s,2:s+1);     % RK coefficients

% initialize as non-embedded, until proven otherwise
embedded = 0;
p = 0;
if (Brows > Bcols)
   if (max(abs(B(s+2,2:s+1))) > eps)
      embedded = 1;
      b2 = (B(s+2,2:s+1))';
      p = B(s+2,1);
   end
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

% create Fdata structure for evaluating solution
Fdata.fname = fcn;    % ODE RHS function name
Fdata.B     = B;      % Butcher table 

% set initial time step size
h = hmin;

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
      Fdata.h = h;
      Fdata.yold = Y0;

      % set Fdata values for this step
      Fdata.t    = t;    % time of last successful step

      % initialize data storage for multiple stages
      z = zeros(m,s);

      % reset step failure flag
      st_fail = 0;
      
      % loop over stages
      for stage=1:s
	 
	 % construct stage solution
	 %    zi = y_n + h*sum_{j=1}^{i-1} (A(i,j)*f(zj))
	 z(:,stage) = Y0;
	 for j=1:stage-1
	    z(:,stage) = z(:,stage) + h*A(stage,j)*feval(fcn,t+h*c(j),z(:,j));
	 end
	 
      end

      % increment number of internal time steps taken
      nsteps = nsteps + 1;

      % compute new solution (and embedding if available)
      [Ynew,Y2] = Y_ERK(z,Fdata);
      
      % if time step adaptivity enabled, check step accuracy
      if (embedded)

	 % estimate error in current step
	 err_step = max(norm((Ynew - Y2)./(rtol*Ynew + atol),inf), eps);
	 
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
      
         % for embedded methods, use error estimate to adapt the time step
	 if (embedded) 
            
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
	 hstab = h_stable * feval(StabFn, t, Ynew);

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
            fprintf('Cannot achieve desired accuracy.\n');
            fprintf('Consider reducing hmin or increasing rtol.\n');
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





function [y,y2] = Y_ERK(z, Fdata)
% usage: [y,y2] = Y_ERK(z, Fdata)
%
% Inputs:
%    z     = stage solutions [z1, ..., zs]
%    Fdata = structure containing extra problem information
%
% Outputs: 
%    y     = step solution built from the z values
%    y2    = embedded solution (if embedding included in Butcher 
%               table; otherwise the same as y)

% extract method information from Fdata
B = Fdata.B;
[Brows, Bcols] = size(B);
s = Bcols - 1;
c = B(1:s,1);
b = (B(s+1,2:s+1))';

% check to see if we have coefficients for embedding
if (Brows > Bcols)
   b2 = (B(s+2,2:s+1))';
else
   b2 = b;
end

% get some problem information
[zrows,zcols] = size(z);
nvar = zrows;
if (zcols ~= s)
   error('Y_ERK error: z has incorrect number of stages');
end

% call RHS at our stages
f = zeros(nvar,s);
for is=1:s
   t = Fdata.t + Fdata.h*c(is);
   f(:,is) = feval(Fdata.fname, t, z(:,is));
end

% form the solutions
%    ynew = yold + h*sum(b(j)*fj)
y  = Fdata.yold + Fdata.h*f*b;
y2 = Fdata.yold + Fdata.h*f*b2;

% end of function
end
