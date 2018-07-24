% driver for Van der Pol ODE test problem:
%    u' = v
%    v' = (v - v*u^2/ep - u
% where u(0) = 2,  v(0) = 0, and ep = 0.2, integrated over
% the time interval [0,12].
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2012
% All Rights Reserved
clear

% set problem parameters
ep = 0.2;
fn = @(t,y) [y(2); y(2)*(1-y(1)^2)/ep - y(1)];
Jn = @(t,y) [0, 1; -1-2*y(1)*y(2)/ep, (1-y(1)^2)/ep];
global Pdata;
Pdata.ep = ep;
Es = @EStab_p1;
Tf = 12;
tout = linspace(0,Tf,100);
hmin = 1e-6;
hmax = 1.0;
rtol = 1e-3;
atol = 1e-14*ones(2,1);
u0 = 2;
v0 = 0;
Y0 = [u0; v0];

% plot "true" solution
opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s(fn, tout, Y0, opts);
figure()
plot(tout,Ytrue)
xlabel('t','FontSize',12), ylabel('y','FontSize',12)
title('Van der Pol test','FontSize',14)
set(gca,'FontSize',12)
print('-dpng','vanderPol')


% run with an embedded explicit RK method
mname = 'Bogacki-Shampine-ERK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with ERK integrator: %s (order = %i)\n',mname,B(s+1,1))
[t,Y,ns,~] = solve_ERK(fn, Es, tout, Y0, B, rtol, atol, hmin, hmax, hmin);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i)\n',ns,ns*s);


% run with a non-embedded explicit RK method
mname = 'ERK-4-4';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with ERK integrator: %s (order = %i)\n',mname,B(s+1,1))
[t,Y,ns,~] = solve_ERK(fn, Es, tout, Y0, B, rtol, atol, hmin, hmax, hmin);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i)\n',ns,ns*s);


% run with an embedded diagonally-implicit RK method
mname = 'TRBDF2-ESDIRK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with DIRK integrator: %s (order = %i)\n',mname,B(s+1,1))
[t,Y,ns,nl,~] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmin);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);


% run with a non-embedded diagonally-implicit RK method
mname = 'Cooper4-ESDIRK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with DIRK integrator: %s (order = %i)\n',mname,B(s+1,1))
[t,Y,ns,nl,~] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmin);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);


% run with a fully-implicit RK method
mname = 'RadauIIA-2-3-IRK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with IRK integrator: %s (order = %i)\n',mname,B(s+1,1))
[t,Y,ns,nl,~] = solve_IRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, hmin);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);


% end of script
