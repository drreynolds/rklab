function [qse, qsi, qsa] = check_ark_stage_order(ce,ci,Ae,Ai,tol)
% Usage: [qse, qsi, qsa] = check_ark_stage_order(ce,ci,Ae,Ai,tol)
% 
% Applies the Butcher simplifying assumption C(q) as prescribed in
%
%J.C. Butcher, "Implicit Runge-Kutta processes", Math Comp
%    18 (1964), pp 50-64.
%
% to determine the stage order(s) of an additive Runge Kutta
% method.  Specifically, the stage order of a RK is based on the
% Butcher table components c and A.  We perform the following
% tests: 
%
%   A. Compute the stage order of the ERK method using (ce,Ae).
%      This output is returned as 'qse'
%
%   B. Compute the stage order of the DIRK method using (ci,Ai).
%      This output is returned as 'qsi'
%
%   C. Compute the stage order of the coupling via (ce,Ai) and
%      (ci,Ae); the minimum of these two values is is returned 
%      as 'qsa'
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2017, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

% check for legal input tolerance
if (tol <= 0.0)
   error('illegal input tolerance, must be > 0')
end

% verify compatability of Butcher table inputs
if ( (length(ce) ~= length(ci)) || ...
     (size(Ae,1) ~= size(Ai,1)) || ...
     (size(Ae,2) ~= size(Ai,2)) || ...
     (length(ce) ~= size(Ae,1)) )
   error('incompatible Butcher table inputs')
end
   
% determine number of stages
s = length(ce);


% perform test for ERK method
qse = 0;
for k=1:1000
   alltrue = 1;
   for i=1:s
      LHS = Ae(i,:)*(ce.^(k-1));
      RHS = ce(i)^k/k;
      if (abs(RHS-LHS)>tol)
         alltrue=0;
         break;
      end
   end
   if (alltrue == 1)
      qse = qse+1;
   else
      break;
   end
end


% perform test for DIRK method
qsi = 0;
for k=1:1000
   alltrue = 1;
   for i=1:s
      LHS = Ai(i,:)*(ci.^(k-1));
      RHS = ci(i)^k/k;
      if (abs(RHS-LHS)>tol)
         alltrue=0;
         break;
      end
   end
   if (alltrue == 1)
      qsi = qsi+1;
   else
      break;
   end
end


% perform tests for ARK method
qsa = 0;
for k=1:1000
   alltrue = 1;
   for i=1:s
      LHS = Ae(i,:)*(ci.^(k-1));
      RHS = ci(i)^k/k;
      if (abs(RHS-LHS)>tol)
         alltrue=0;
         break;
      end
      LHS = Ai(i,:)*(ce.^(k-1));
      RHS = ce(i)^k/k;
      if (abs(RHS-LHS)>tol)
         alltrue=0;
         break;
      end
   end
   if (alltrue == 1)
      qsa = qsa+1;
   else
      break;
   end
end

% end of function
