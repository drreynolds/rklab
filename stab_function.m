function [alpha,beta] = stab_function(A,b)
% usage: [alpha,beta] = stab_function(A,b)
%
% Computes the coefficients of the polynomials that comprise the
% numerator and denominator of the rational RK stability function
% defined by Butcher table matrix A and solution coefficients b.
%
% The coefficients alpha and beta are listed in order of increasing
% degree, i.e.
%
%   R(z) = num(z) / den(z),
%   num(z) = alpha(1) + alpha(2)*z + ...
%   den(z) = beta(1) + beta(2)*z + ...
%
% Daniel R. Reynolds
% SMU Mathematics
% December 2016

% check A and b for compatibility
[rows,cols] = size(A);
s = length(b);
if ((s ~= rows) || (s ~= cols))
    error('stab_function: incompatible Butcher table inputs')
end

% convert b to a row vector
b = reshape(b,1,s);

% create stability function
syms z;
e = ones(s,1);
R = 1 + z*b*((eye(s)-z*A)\e);

% simplify the rational function
Rsimp = simplifyFraction(R, 'Expand', true);

% extract numerator and denominator polynomials
[n,d] = numden(Rsimp);

% extract coefficients of each polynomial 
% (sorted from lowest to highest degree)
alpha = fliplr(coeffs(n,z,'All'));
beta = fliplr(coeffs(d,z,'All'));

% scale both so that constant terms are at most 1
if (abs(alpha(1)) > abs(beta(1)))
    beta = beta / alpha(1);
    alpha = alpha / alpha(1);
else
    alpha = alpha / beta(1);
    beta = beta / beta(1);
end

% remove any trailing zeros for non-symbolic coefficients
if (~isa(alpha(1),'sym'))
   for i=1:length(alpha)-1
      if (abs(alpha(end)) < 1e-18)
         alpha = alpha(1:end-1);
      else
         break;
      end
   end
   for i=1:length(beta)-1
      if (abs(beta(end)) < 1e-18)
         beta = beta(1:end-1);
      else
         break;
      end   
   end
end
