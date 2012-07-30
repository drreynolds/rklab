function F = F_IRK(z, Fdata)
% usage: F = F_IRK(z, Fdata)
%
% Inputs:  z = current guesses for [z1, ..., zs]
%          Fdata = structure containing extra information for evaluating F.
% Outputs: F = residual at current guess
%
% This function computes the (non)linear residuals for each intermediate
% stage solution, through calling the user-supplied (in Fdata) ODE
% right-hand side function.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2012
% All Rights Reserved

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
   f(:,is) = feval(Fdata.fname, t, z(:,is));
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
