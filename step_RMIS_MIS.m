function [Y,Yerr,m,hf] = step_RMIS_MIS(fs,ff,Jf,t0,Y0,Bo,Bi,hs,hf,rtol,atol)
% usage: [Y,Yerr,m,hf] = step_RMIS_MIS(fs,ff,Jf,t0,Y0,Bo,Bi,hs,hf,rtol,atol)
%
% This routine performs a single step of the relaxed multirate
% infinitesimal step (RMIS) method for the vector-valued ODE problem
%     y' = fs(t,Y) + ff(t,Y), t >= t0, y in R^n,
%     Y0 = [y1(t0), y2(t0), ..., yn(t0)]'.
% We perform the solve using a variation of the approach used for
% MIS methods by Knoth & Wolke (1998), i.e., we do NOT consider the
% problem in full GARK form.  We do this to compute both an MIS and
% an RMIS solution, with the difference constituting the temporal
% error estimate.
%
% Inputs:
%     fs     = function handle for (slow) ODE RHS
%     ff     = function handle for (fast) ODE RHS
%     Jf     = function handle for Jacobian of ff (only required if an
%              implicit inner method is supplied)
%     t0     = value of independent variable at start of step
%     Y0     = solution vector at start of step (column vector of length n)
%     Bo     = ERK Butcher table for 'outer' (slow) method
%                 Bo = [co Ao;
%                       qo bo;
%                       po do ]
%              Here, co is a vector of stage time fractions (so-by-1),
%                    Ao is a matrix of Butcher coefficients (so-by-so),
%                    qo is an integer denoting the method order of accuracy
%                    bo is a vector of solution weights (1-by-so),
%                    po is an integer denoting the embedding order of accuracy,
%                    do is a vector of embedding weights (1-by-so),
%              The [po, do] row is optional, and is unused in this
%              routine.  Also, the qo value is unused here.
%              Note: the MIS method assumes that the outer method stage
%              times are sorted, i.e. co(i+1) > co(i)  for all i.
%     Bi     = Butcher table for a single step of an 'inner' (fast) method
%              (can be ERK, DIRK or IRK)
%                 Bi = [ci Ai;
%                       qi bi;
%                       pi di ]
%              All components have the same role as with Bi; we
%              assume that Bi encodes a method with si stages
%     hs     = step size to use for slow time scale
%     hf     = initial step size to use for fast time scale
%     rtol   = desired relative error of solution at the fast scale (scalar)
%     atol   = desired absolute error of solution at the fast scale (vector or scalar)
%
% Outputs:
%     Y      = updated RMIS solution, [y1(t0+hs), y2(t0+hs), ..., yn(t0+hs)].
%     Yerr   = multirate temporal error estimate.
%     m      = actual number of substeps used for inner method.
%     hf     = final step size used at fast time scale
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% July 2018
% All Rights Reserved

% dummy stability function, tolerances
estab = @(t,y) inf;

% extract ERK method information from Bo, and ensure that it is legal
[Brows, Bcols] = size(Bo);
so = Bcols - 1;           % number of stages
co = Bo(1:so,1);          % stage time fraction array
bo = (Bo(so+1,2:so+1))';  % solution weights (convert to column)
Ao = Bo(1:so,2:so+1);     % ERK coefficients
if (max(max(abs(triu(Ao)))) > 0)
   error('Error: Bo does not specify an explicit RK method table')
end

% extract RK method information from Bi
[Brows, Bcols] = size(Bi);
si = Bcols - 1;           % number of stages
ci = Bi(1:si,1);          % stage time fraction array
bi = (Bi(si+1,2:si+1))';  % solution weights (convert to column)
Ai = Bi(1:si,2:si+1);     % RK coefficients
innerRK = 0;
if (max(max(abs(triu(Ai)))) > 0)        % implicit components exist
   innerRK = 1;
   if (max(max(abs(triu(Ai,1)))) > 0)   % method is IRK
      innerRK = 2;
   end
end

% initialize outputs
m = 0;
n = length(Y0);
Y = reshape(Y0,n,1);

% initialize temporary variables
Fs = zeros(n,so);
fscale = zeros(so,1);
Ys = Y;

% first outer stage
Fs(:,1) = fs(t0,Y);
tcur = t0;

% add contributions from first outer stage to overall solution.
% Note: since inner method has explicit first stage, then its
% contribution to overall solution may be obtained by evaluating
% fast RHS at the same time as the slow
Y = Y + hs*bo(1)*(Fs(:,1) + ff(t0,Y));

% iterate over remaining outer stages
for i=2:so

   % determine 'inner' ODE for this stage
   %   RHS function
   for j=1:i-1
      fscale(j) = (Ao(i,j)-Ao(i-1,j))/(co(i)-co(i-1));
   end
   fi = @(t,y) ff(t,y) + Fs*fscale;
   %   time interval
   tspan = [tcur, tcur + (co(i)-co(i-1))*hs];
   %   step size bounds for fast scale
   hf_max = (co(i)-co(i-1))*hs;
   hf_min = hf_max * 1e-5;

   % call inner RK method solver to perform substepping
   if (innerRK == 2)         % IRK inner method
      [tvals, V, mi, ~, hf] = solve_IRK(fi, Jf, tspan, Ys, Bi, rtol, atol, hf_min, hf_max, hf);
   elseif (innerRK == 1)     % DIRK inner method
      [tvals, V, mi, ~, hf] = solve_DIRK(fi, Jf, tspan, Ys, Bi, rtol, atol, hf_min, hf_max, hf);
   else                      % ERK inner method
      [tvals, V, mi, hf] = solve_ERK(fi, estab, tspan, Ys, Bi, rtol, atol, hf_min, hf_max, hf);
   end
   m = m + mi;

   % update slow 'solution' as result from fast solve
   Ys = V(:,end);
   tcur = t0 + c(i)*hs;
   Fs(:,i) = fs(tcur,Ys);

   % update overall solution; note that since inner method has
   % explicit first stage, then its contribution to overall
   % solution may be obtained by evaluating fast RHS at the same
   % time as the slow
   Y = Y + hs*bo(i)*(Fs(:,i) + ff(tcur,Ys));

end


% all slow stages completed, if any of the time interval remains,
% finish that off here
if (c(so) < 1)

   % determine 'inner' ODE for this stage
   %   RHS function
   for j=1:so
      fscale(j) = (bo(j)-Ao(so,j))/(1-co(so));
   end
   fi = @(t,y) ff(t,y) + Fs*fscale;
   %   time interval
   tspan = [tcur, hs];
   %   num internal time steps
   ni = ceil((1-co(so))*hs/hf);
   %   step size
   hf = (1-co(so))*hs/ni;

   % call inner RK method solver to perform substepping
   if (innerRK == 2)         % IRK inner method
      [tvals, V, mi, ~, hf] = solve_IRK(fi, Jf, tspan, Ys, Bi, rtol, atol, hf_min, hf_max, hf);
   elseif (innerRK == 1)     % DIRK inner method
      [tvals, V, mi, ~, hf] = solve_DIRK(fi, Jf, tspan, Ys, Bi, rtol, atol, hf_min, hf_max, hf);
   else                      % ERK inner method
      [tvals, V, mi, hf] = solve_ERK(fi, estab, tspan, Ys, Bi, rtol, atol, hf_min, hf_max, hf);
   end
   m = m + mi;

   % update slow 'solution' as result from fast solve
   Ys = V(:,end);

end

Yerr = Y-Ys;


% end of function
