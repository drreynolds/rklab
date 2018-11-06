% driver for brusselator PDE test problem:
%      u' = d1*u_xx + a - (b+1)*u + u^2*v,
%      v' = d2*v_xx + b*u - u^2*v,
% with parameters a=0.6, b=2, d1=d2=0.25, over the spatial
% interval [0,1] and the time interval [0,80].  Initial
% conditions are u(x,0) = a + 0.1*sin(pi*x),
% v(x,0) = b/a + 0.1*sin(pi*x), and boundary conditions are
% homogeneous Dirichlet.  We use a mesh with 100 spatial zones.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% March 2017
% All Rights Reserved
clear

% set problem parameters
global a b d1 d2 m dx
a = 0.6;
b = 2;
d1 = 0.025;
d2 = 0.025;
m = 100;
dx = 1/(m-1);
fn = @f_p3;
fe = @fe_p3;
fi = @fi_p3;
Jn = @J_p3;
Ji = @Ji_p3;
Es = @EStab_p3;
Tf = 10;
tout = linspace(0,Tf,100);
hmin = 1e-6;
hmax = 1.0;
rtol = 1e-4;
atol = 1e-14*ones(2*m,1);

% initial conditions
xspan = linspace(0,1,m)';
u0 = a + 0.1*sin(pi*xspan);
v0 = b/a + 0.1*sin(pi*xspan);
Y0 = [u0; v0];


% plot "true" solution at end of run
opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s(fn, tout, Y0, opts);
figure()
plot(xspan,Ytrue(end,1:m),xspan,Ytrue(end,m+1:2*m))
xlabel('t','FontSize',12), ylabel('y','FontSize',12)
title('1D Brusselator PDE test','FontSize',14)
set(gca,'FontSize',12)
print('-dpng','brusselator1D')


% run with an embedded diagonally-implicit RK method
mname = 'Kvaerno(7,4,5)-ESDIRK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with DIRK integrator: %s (order = %i)\n',mname,B(s+1,1))
[t,Y,ns,nl,~] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmin);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);


% run with a non-embedded diagonally-implicit RK method
mname = 'Cooper6-ESDIRK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with DIRK integrator: %s (order = %i)\n',mname,B(s+1,1))
[t,Y,ns,nl,~] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmin);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);


% run with a fully-implicit RK method
mname = 'RadauIIA-3-5-IRK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with IRK integrator: %s (order = %i)\n',mname,B(s+1,1))
[t,Y,ns,nl,~] = solve_IRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmin);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);


% run with an embedded ARK method
mname1 = 'ARK5(4)8L[2]SA-ERK';
Be = butcher(mname1);  s = numel(Be(1,:))-1;
mname2 = 'ARK5(4)8L[2]SA-ESDIRK';
Bi = butcher(mname2);  s = numel(Bi(1,:))-1;
fprintf('\nRunning with ARK integrator: %s/%s (order = %i)\n',...
        mname1,mname2,Be(s+1,1))
[t,Y,ns,nl] = solve_ARK(fe, fi, Ji, tout, Y0, Be, Bi, rtol, atol, hmin, hmax, hmin);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);


% $$$ % run with an embedded explicit RK method
% $$$ mname = 'Merson-5-4-ERK';
% $$$ B = butcher(mname);  s = numel(B(1,:))-1;
% $$$ fprintf('\nRunning with ERK integrator: %s (order = %i)\n',mname,B(s+1,1))
% $$$ [t,Y,ns,~] = solve_ERK(fn, Es, tout, Y0, B, rtol, atol, hmin, hmax, hmin);
% $$$ err_max = max(max(abs(Y'-Ytrue)));
% $$$ err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
% $$$ fprintf('Accuracy/Work Results:\n')
% $$$ fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
% $$$ fprintf('   steps = %i (stages = %i)\n',ns,ns*s);
% $$$
% $$$
% $$$ % run with a non-embedded explicit RK method
% $$$ mname = 'ERK-4-4';
% $$$ B = butcher(mname);  s = numel(B(1,:))-1;
% $$$ fprintf('\nRunning with ERK integrator: %s (order = %i)\n',mname,B(s+1,1))
% $$$ [t,Y,ns,~] = solve_ERK(fn, Es, tout, Y0, B, rtol, atol, hmin, hmax, hmin);
% $$$ err_max = max(max(abs(Y'-Ytrue)));
% $$$ err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
% $$$ fprintf('Accuracy/Work Results:\n')
% $$$ fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
% $$$ fprintf('   steps = %i (stages = %i)\n',ns,ns*s);

% end of script



%------------------------- Utility routines -------------------------%



function dy = f_p3(t, y)

   % extract solution components
   global a b d1 d2 m dx
   u = y(1:m);
   v = y(m+1:2*m);

   % initialize RHS terms
   du = zeros(m,1);
   dv = zeros(m,1);

   % enforce stationary boundary conditions
   du(1) = 0;  du(m) = 0;  dv(1) = 0;  dv(m) = 0;

   % diffusion components
   du(2:m-1) = d1/dx/dx*(u(3:m)+u(1:m-2)-2*u(2:m-1));
   dv(2:m-1) = d2/dx/dx*(v(3:m)+v(1:m-2)-2*v(2:m-1));

   % reaction components
   du(2:m-1) = du(2:m-1) + a - (b+1)*u(2:m-1) + u(2:m-1).*u(2:m-1).*v(2:m-1);
   dv(2:m-1) = dv(2:m-1) + b*u(2:m-1) - u(2:m-1).*u(2:m-1).*v(2:m-1);

   % combine together into output
   dy = [du; dv];

   % end function
end



function dy = fe_p3(t, y)

   % extract solution components
   global a b d1 d2 m dx
   u = y(1:m);
   v = y(m+1:2*m);

   % initialize RHS terms
   du = zeros(m,1);
   dv = zeros(m,1);

   % enforce stationary boundary conditions
   du(1) = 0;  du(m) = 0;  dv(1) = 0;  dv(m) = 0;

   % reaction components
   du(2:m-1) = a - (b+1)*u(2:m-1) + u(2:m-1).*u(2:m-1).*v(2:m-1);
   dv(2:m-1) = b*u(2:m-1) - u(2:m-1).*u(2:m-1).*v(2:m-1);

   % combine together into output
   dy = [du; dv];

   % end function
end



function dy = fi_p3(t, y)

   % extract solution components
   global a b d1 d2 m dx
   u = y(1:m);
   v = y(m+1:2*m);

   % initialize RHS terms
   du = zeros(m,1);
   dv = zeros(m,1);

   % enforce stationary boundary conditions
   du(1) = 0;  du(m) = 0;  dv(1) = 0;  dv(m) = 0;

   % diffusion components
   du(2:m-1) = d1/dx/dx*(u(3:m)+u(1:m-2)-2*u(2:m-1));
   dv(2:m-1) = d2/dx/dx*(v(3:m)+v(1:m-2)-2*v(2:m-1));

   % combine together into output
   dy = [du; dv];

   % end function
end



function J = J_p3(t, y)

   % extract solution components
   global a b d1 d2 m dx
   u = y(1:m);
   v = y(m+1:2*m);

   % initialize Jacobian blocks terms
   Juu = sparse([],[],[],m,m,3*m);
   Juv = sparse([],[],[],m,m,m);
   Jvu = sparse([],[],[],m,m,m);
   Jvv = sparse([],[],[],m,m,3*m);

   % diffusion components
   for j=2:m-1
      Juu(j,j-1) = Juu(j,j-1) + d1/dx/dx;
      Juu(j,j)   = Juu(j,j) - 2*d1/dx/dx;
      Juu(j,j+1) = Juu(j,j+1) + d1/dx/dx;
   end
   for j=2:m-1
      Jvv(j,j-1) = Jvv(j,j-1) + d2/dx/dx;
      Jvv(j,j)   = Jvv(j,j) - 2*d2/dx/dx;
      Jvv(j,j+1) = Jvv(j,j+1) + d2/dx/dx;
   end

   % reaction components
   for j=2:m-1
      Juu(j,j) = Juu(j,j) - (b+1) + 2*u(j)*v(j);
   end
   for j=2:m-1
      Juv(j,j) = Juv(j,j) + u(j)^2;
   end
   for j=2:m-1
      Jvv(j,j) = Jvv(j,j) - u(j)^2;
   end
   for j=2:m-1
      Jvu(j,j) = Jvu(j,j) + b - 2*u(j)*v(j);
   end

   % combine together into output
   J = [Juu, Juv; Jvu, Jvv];

   % end function
end



function J = Ji_p3(t, y)

   % extract solution components
   global a b d1 d2 m dx
   u = y(1:m);
   v = y(m+1:2*m);

   % initialize Jacobian blocks terms
   Juu = sparse([],[],[],m,m,3*m);
   Juv = sparse([],[],[],m,m,m);
   Jvu = sparse([],[],[],m,m,m);
   Jvv = sparse([],[],[],m,m,3*m);

   % diffusion components
   for j=2:m-1
      Juu(j,j-1) = Juu(j,j-1) + d1/dx/dx;
      Juu(j,j)   = Juu(j,j) - 2*d1/dx/dx;
      Juu(j,j+1) = Juu(j,j+1) + d1/dx/dx;
   end
   for j=2:m-1
      Jvv(j,j-1) = Jvv(j,j-1) + d2/dx/dx;
      Jvv(j,j)   = Jvv(j,j) - 2*d2/dx/dx;
      Jvv(j,j+1) = Jvv(j,j+1) + d2/dx/dx;
   end

   % combine together into output
   J = [Juu, Juv; Jvu, Jvv];

   % end function
end


function dt = EStab_p3(t, y)

   % compute the Jacobian ODE RHS
   J = J_p3(t,y);

   % determine the largest eigenvalue magnitude
   lam = max(abs(eigs(J)));

   % assume explicit stability region includes Euler stability region
   % (this assumes that the eigenvalue is in fact negative).
   dt = 1/lam;

   % end function
end
