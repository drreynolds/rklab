function [q,p,Bs,As,Ls,BsE,AsE,LsE] = check_butcher(B)
% Usage: [q,p,Bs,As,Ls,BsE,AsE,LsE] = check_butcher(B)
% 
% Checks the Butcher table B to determine the analytical order of
% accuracy for the method (q) and its embedding (p), whether the
% method and embedding are B stable, and estimates of whether the
% method and embedding are A and/or L stable.  It is assumed that B
% has block structure 
%     B = [c, A; 0, b]
% for a standard Runge-Kutta method, or 
%     B = [c, A; 0, b; 0, b2]
% if the method has an embedded error indicator.  
%
% If the method has no embedding, then we set p=0.
%
% This function uses Butcher's "simplifying assumptions"
% (i.e. sufficient order condition equations) to check the order of
% accuracy for the method, from the article
%
%    J.C. Butcher, "Implicit Runge-Kutta processes", Math Comp
%    18 (1964), pp 50-64.
%
% For B-stability, we check for positive semi-definiteness of 
%    M(i,j) = [b(i)*A(i,j) + b(j)*A(j,i) - b(i)*b(j)]
% as outlined in 
%    K. Burrage and J.C. Butcher, "Stability criteria for implicit
%    Runge-Kutta methods", SIAM J. Numer. Anal. 16 (1979), pp
%    46-57.
%
% For A-stability we perform stability tests on the stability
% function R(z) for 1000 random values in the left half-plane.
%
% For L-stability we check that R(z) decreases monotonically in the
% left half-plane as z increases in magnitude from -1 to -1e8.
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

% set tolerance on 'equality'
tol = 1e-8;

% extract components of Butcher table
[m,n] = size(B);
if (m == n)        % no embedding
   d = 0;
elseif (m == n+1)  % has an embedding
   d = B(m,2:n)';
else   % illegal input
   error('illegal Butcher table input')
end
s = n-1;
c = B(1:s,1);
b = B(s+1,2:n)';
A = B(1:s,2:n);

% determine P, Q and R for the Butcher simplifying assumptions
%   B(P):
P = 0;
for i=1:1000
   LHS = b'*(c.^(i-1));
   RHS = 1/i;
   if (abs(RHS-LHS)>tol)
      break;  
   end
   P = P+1;
end

%   C(Q):
Q = 0;
for k=1:1000
   alltrue = 1;
   for i=1:s
      LHS = A(i,:)*(c.^(k-1));
      RHS = c(i)^k/k;
      if (abs(RHS-LHS)>tol)
         alltrue=0;
         break;
      end
   end
   if (alltrue == 1)
      Q = Q+1;
   else
      break;
   end
end

%   D(R):
R = 0;
for k=1:1000
   alltrue = 1;
   for j=1:s
      LHS = 0;
      for i=1:s
         LHS = LHS + A(i,j)*b(i)*c(i)^(k-1);
      end
      RHS = b(j)/k*(1-c(j)^k);
      if (abs(RHS-LHS)>tol)
         alltrue = 0;
         break;  
      end
   end
   if (alltrue == 1)
      R = R+1;
   else
      break;
   end
end

% determine q
q = 0;
for i=1:P
   if ((q > Q+R+1) || (q > 2*Q+2)),  
      break;  
   end
   q = q+1;
end


% if there's an embedding, determine the order
if (length(d) > 1) 
   P = 0;
   for i=1:1000
      LHS = d'*(c.^(i-1));
      RHS = 1/i;
      if (abs(RHS-LHS)>tol),  break;  end
      P = P+1;
   end
   R = 0;
   for k=1:1000
      alltrue = 1;
      for j=1:s
         LHS = 0;
         for i=1:s
            LHS = LHS + A(i,j)*d(i)*c(i)^(k-1);
         end
         RHS = d(j)/k*(1-c(j)^k);
         if (abs(RHS-LHS)>tol)
            alltrue = 0;
            break;  
         end
      end
      if (alltrue == 1)
         R = R+1;
      else
         break;
      end
   end
   p = 0;
   for i=1:P
      if ((p > Q+R+1) || (p > 2*Q+2)),  break;  end
      p = p+1;
   end
  
else
   p = 0;
end



% determine B stability
M = zeros(s,s);
for j=1:s
   for i=1:s
      M(i,j) = b(i)*A(i,j) + b(j)*A(j,i) - b(i)*b(j);
   end
end
lam = eig(M);
Bs = 1;
for i=1:s
   if (lam(i) < -tol)
      Bs = 0;
      break;
   end
end


% estimate A stability:
%   Hairer & Wanner, Solving ODEs II: an RK method is A-stable iff
%      (i)  |R(iy)| < 1  for all y in R
%      (ii) R(z) is analytic in the left half-plane
%   generate stability function coefficients (reverse ordering)
[alpha,beta] = stab_function(A,b);
alpha_dbl = double(alpha);
beta_dbl = double(beta);
%   remove zero coefficients for highest-order terms
for i=1:length(alpha_dbl)-1
   if (abs(alpha_dbl(end)) < eps)
      alpha_dbl = alpha_dbl(1:end-1);
   else
      break;
   end
end
for i=1:length(beta_dbl)-1
   if (abs(beta_dbl(end)) < eps)
      beta_dbl = beta_dbl(1:end-1);
   else
      break;
   end
end
   
%   check analytic in left half-plane
beta2 = fliplr(beta_dbl);
alpha2 = fliplr(alpha_dbl);
rt = roots(beta2);
As = 1;
for i=1:length(rt)
    if (real(rt(i)) < -tol)
        As=0;
        break
    end
end
%   check along imaginary axis
ztests = [0, sqrt(-1)*logspace(-5,4,10000)];
for j = 1:length(ztests)
   z = ztests(j);
   if (abs(polyval(alpha2,z)/polyval(beta2,z))>1+tol)
      As = 0;
      break;
   end
end


% verify that z = -0.01 is in the stability region; if not, report
% A-stability as "-1"
z = -0.01;
if (abs(polyval(alpha2,z)/polyval(beta2,z))>1)
    As = -1;
end


% check L stability:
%   if degree of denominator is greater than degree of numerator in
%   ational stability polynomial, then it will be L-stable
numdeg = length(alpha)-1;
for i=0:length(alpha)-1
   if (abs(alpha(end-i)) < tol^2)
      numdeg = numdeg-1;
   end
end
dendeg = length(beta)-1;
for i=0:length(beta)-1
   if (abs(beta(end-i)) < tol^2)
      dendeg = dendeg-1;
   end
end
if (dendeg > numdeg)
   Ls = 1;
else
   Ls = 0;
end


% determine B stability of embedding
BsE = 0;
if (length(d) > 1) 
   M = zeros(s,s);
   for j=1:s
      for i=1:s
         M(i,j) = d(i)*A(i,j) + d(j)*A(j,i) - d(i)*d(j);
      end
   end
   lam = eig(M);
   BsE = 1;
   for i=1:s
      if (lam(i) < -tol)
         BsE = 0;
         break;
      end
   end
end

% estimate A stability of embedding:
%   Hairer & Wanner, Solving ODEs II: an RK method is A-stable iff
%      (i)  |R(iy)| < 1  for all y in R
%      (ii) R(z) is analytic in the left half-plane
%   generate stability function coefficients (reverse ordering)
AsE = 0;
if (length(d) > 1)
   [alpha,beta] = stab_function(A,d);
   alpha_dbl = double(alpha);
   beta_dbl = double(beta);
   %   remove zero coefficients for highest-order terms
   for i=1:length(alpha_dbl)-1
      if (abs(alpha_dbl(end)) < eps)
         alpha_dbl = alpha_dbl(1:end-1);
      else
         break;
      end
   end
   for i=1:length(beta_dbl)-1
      if (abs(beta_dbl(end)) < eps)
         beta_dbl = beta_dbl(1:end-1);
      else
         break;
      end
   end
   
   %   check analytic in left half-plane
   beta2 = fliplr(beta_dbl);
   alpha2 = fliplr(alpha_dbl);
   rt = roots(beta2);
   AsE = 1;
   for i=1:length(rt)
      if (real(rt(i)) < -tol)
         AsE=0;
         break
      end
   end
   %   check along imaginary axis
   ztests = [0, sqrt(-1)*logspace(-5,4,10000)];
   for j = 1:length(ztests)
      z = ztests(j);
      if (abs(polyval(alpha2,z)/polyval(beta2,z))>1+tol)
         AsE = 0;
         break;
      end
   end

   % verify that z = -0.01 is in the stability region; if not, report
   % A-stability as "-1"
   z = -0.01;
   if (abs(polyval(alpha2,z)/polyval(beta2,z))>1)
      AsE = -1;
   end
end


% check L stability of embedding:
%   if degree of denominator is greater than degree of numerator in
%   ational stability polynomial, then it will be L-stable
LsE = 0;
if (length(d) > 1)
   numdeg = length(alpha)-1;
   for i=0:length(alpha)-1
      if (abs(alpha(end-i)) < tol^2)
         numdeg = numdeg-1;
      end
   end
   dendeg = length(beta)-1;
   for i=0:length(beta)-1
      if (abs(beta(end-i)) < tol^2)
         dendeg = dendeg-1;
      end
   end
   if (dendeg > numdeg)
      LsE = 1;
   end
end


% end of function
