function [tvals,Y,nsteps,lits] = solve_IRK(fcn,Jfcn,tvals,Y0,B,rtol,atol,hmin,hmax)
% usage: [tvals,Y,nsteps,lits] = solve_IRK(fcn,Jfcn,tvals,Y0,B,rtol,atol,hmin,hmax)
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
%     lits   = number of linear solves required by method
%
% Note: currently rtol and atol are ignored, and the solver only
% operates in fixed-step mode, using steps of size hmin.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2012
% All Rights Reserved

   
% get number of stages for IRK method
[Brows, Bcols] = size(B);
s = Bcols - 1;

% extract order of accuracy for method
p = B(s+1,1);

% initialize output arrays
N = length(tvals);
m = length(Y0);
Y = zeros(m,N);
Y(:,1) = Y0;

% set the solver parameters
newt_maxit = 20;           % max number of Newton iterations
newt_tol   = 1e-10;        % Newton solver tolerance
newt_alpha = 1;            % Newton damping parameter
ONEMSM     = 1-sqrt(eps);  % coefficients to account for
ONEPSM     = 1+sqrt(eps);  %   floating-point roundoff

% initialize temporary variables
t = tvals(1);
Ynew = Y0;

% create Fdata structure for Newton solver and step solutions
Fdata.fname = fcn;    % ODE RHS function name
Fdata.Jname = Jfcn;   % ODE RHS Jacobian function name
Fdata.B     = B;      % Butcher table 
Fdata.s     = s;      % number of stages

% initialize work counters
nsteps = 0;
lits   = 0;

% iterate over time steps
for tstep = 2:length(tvals)

   % loop over internal time steps to get to desired output time
   while (t < tvals(tstep)*ONEMSM)
      
      % set internal time step
      h = max([h, hmin]);            % enforce minimum time step size
      h = min([h, hmax]);            % maximum time step size
      h = min([h, tvals(tstep)-t]);  % stop at output time

      
      % set Fdata values for this step
      Fdata.h    = h;      % current step size
      Fdata.yold = Ynew;   % solution from previous step
      Fdata.t    = t;      % time of last successful step
      
      % set Newton initial guesses as previous step solution
      z = zeros(s*m,1);
      for i = 0:s-1
	 z(i*m+1:(i+1)*m) = Ynew;
      end

      % call Newton solver to update solution in time
      [z,lin,ierr] = newton_damped('F_IRK', 'A_IRK', z, Fdata, ...
                                   newt_tol, newt_maxit, newt_alpha);

      % increment total linear solver statistics
      lits = lits + lin;
      
      % compute new solution
      Ynew = Y_IRK(z,Fdata);
      
      % check newton error flag, if failure return with error
      if (ierr ~= 0) 
	 error('Newton failure: consider reducing hmin and retry.');
      end
	 
      % update time, time step counter
      t = t + h;
      nsteps = nsteps + s;
   
   end

   % store updated solution in output array
   Y(:,tstep) = Ynew;
   
end


% end function
end





function y = Y_IRK(z, Fdata)
% usage: [y,y2] = Y_IRK(z, Fdata)
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
   f(:,is) = feval(Fdata.fname, t, z(:,is));
end

% form the solution
%    ynew = yold + h*sum(b(j)*fj)
y = Fdata.yold + Fdata.h*f*b;

% end of function
end

