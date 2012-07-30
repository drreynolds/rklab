function y = Y_IRK(z, Fdata)
% usage: y = Y_IRK(z, Fdata)
%
% Inputs:  z = current guesses for [z1, ..., zs]
%          Fdata = structure containing extra information for evaluating F.
% Outputs: y = glued-together time-evolved solution
%
% This routine takes as input the intermediate-time states (z) for a
% multi-stage IRK method, and pieces them together to form the time-evolved
% solution y(t_{n+1}). 
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
   error('Y_IRK error: input argument has incorrect length (must be a multiple of s)');
end

% reshape our z arguments
z = reshape(z,nvar,s);

% call f at our guesses
f = zeros(nvar,s);
for is=1:s
   t = Fdata.t + Fdata.h*c(is);
   f(:,is) = feval(Fdata.fname, t, z(:,is));
end

% form the solution
%    ynew = yold + h*sum(b(j)*fj)
y = Fdata.yold + Fdata.h*f*b;

% end of function
