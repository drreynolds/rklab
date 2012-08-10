% driver for stiff brusselator test problem:
%      u' = a - (w+1)*u + u^2*v,
%      v' = w*u - u^2*v,
%      w' = (b-w)/ep - u*w,
% where u(0) = 1.2, v(0) = 3.1 and w(0) = 3, with prameters a=1,
% b=3.5 and ep=5e-6.  We evaluate over the time interval [0,10].  
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2012
% All Rights Reserved
clear

% set problem parameters
fn = 'f_p2';
Jn = 'J_p2';
Es = 'EStab_p2';
Tf = 10;
tout = linspace(0,Tf,100);
hmin = 1e-7;
hmax = 1.0;
rtol = 1e-3;
atol = 1e-14*ones(3,1);
global Pdata;
Pdata.a = 1; 
Pdata.b = 3.5; 
%Pdata.ep = 5e-6;
Pdata.ep = 1e-3;
u0 = 1.2;
v0 = 3.1;
w0 = 3;
Y0 = [u0; v0; w0];

% plot "true" solution 
opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s(fn, tout, Y0, opts);
figure()
plot(tout,Ytrue)
xlabel('t','FontSize',12), ylabel('y','FontSize',12)
title('Brusselator ODE test','FontSize',14)
set(gca,'FontSize',12)
print('-dpng','brusselator')


% run with a diagonally-implicit RK method
mname = 'Cash(5,3,4)-SDIRK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with DIRK integrator: %s (order = %i)\n',mname,B(s+1,1))
[t,Y,ns,nl] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);


% run with a fully-implicit RK method
mname = 'LobattoIIIC-3-4-IRK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with IRK integrator: %s (order = %i)\n',mname,B(s+1,1))
[t,Y,ns,nl] = solve_IRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);


% run with an explicit RK method
mname = 'Fehlberg-ERK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with ERK integrator: %s (order = %i)\n',mname,B(s+1,1))
[t,Y,ns] = solve_ERK(fn, Es, tout, Y0, B, rtol, atol, hmin, hmax);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i)\n',ns,ns*s);



% end of script
