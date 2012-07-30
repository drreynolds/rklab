function Amat = A_DIRK(z, Fdata)
% usage: Amat = A_DIRK(z, Fdata)
%
% Inputs:  z = current guess for stage solution
%          Fdata = structure containing extra information for evaluating F.
% Outputs: Amat = Jacobian at current guess
%
% This function computes the Jacobian of each intermediate stage residual
% for a multi-stage DIRK method, through calling the user-supplied (in
% Fdata) ODE Jacobian function. 
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2012
% All Rights Reserved

% extract DIRK method information from Fdata
B = Fdata.B;
[Brows, Bcols] = size(B);
s  = Bcols - 1;
c  = B(1:s,1);
b  = (B(s+1,2:s+1))';
A  = B(1:s,2:s+1);
st = Fdata.stage;
t  = Fdata.t + Fdata.h*c(st);

% form the DIRK Jacobian
Amat = eye(length(z)) - Fdata.h*A(st,st)*feval(Fdata.Jname, t, z);

% end of function
