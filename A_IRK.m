function Amat = A_IRK(z, Fdata)
% usage: Amat = A_IRK(z, Fdata)
%
% Inputs:  z = current guesses for [z1, ..., zs]
%          Fdata = structure containing extra information for evaluating F.
% Outputs: Amat = Jacobian at current guess
%
% This function computes the Jacobian of each intermediate stage residual
% for a multi-stage IRK method, through calling the user-supplied (in Fdata)
% ODE Jacobian function. 
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
   error('A_IRK error: input argument has incorrect length (must be a multiple of s)');
end

% reshape our z arguments
z = reshape(z,nvar,s);

% call J at each of our guesses
J = zeros(nvar,nvar,s);
for is=1:s
   t = Fdata.t + Fdata.h*c(is);
   J(:,:,is) = feval(Fdata.Jname, t, z(:,is));
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
