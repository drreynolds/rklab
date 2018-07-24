% multirate driver for Van der Pol ODE test problem:
%    [u'] = [  v ]  + [     0     ] = fslow(u,v) + ffast(u,v)
%    [v']   [ -u ]    [v(1-u^2)/ep]
% where u(0) = 2,  v(0) = 0, and ep = 0.2, integrated over
% the time interval [0,12].
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% July 2018
% All Rights Reserved
clear

% set problem parameters
ep = 0.2;
fs  = @(t,y) [y(2); -y(1)];
ff  = @(t,y) [0; y(2)*(1 - y(1)^2)/ep];
fn  = @(t,y) fs(t,y) + ff(t,y);
Jff = @(t,y) [0, 0; -2*y(1)*y(2)/ep, (1-y(1)^2)/ep];
Tf = 12;
tout = linspace(0,Tf,100);
hmin = 1e-6;
hmax = 1.0;
ntests = 6;
rtol = 1e-3;
atol = 1e-14*ones(2,1);
u0 = 2;
v0 = 0;
Y0 = [u0; v0];

% construct reference solution
opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s(fn, tout, Y0, opts);

% perform convergence tests with fixed-step multirate solvers
solvers = {@solve_MIS_KW3, @solve_MIS_38, @solve_RMIS_KW3, @solve_RMIS_38};
snames = {'MIS-KW3', 'MIS-3/8', 'RMIS-KW3', 'RMIS-3/8'};
orders = [3, 3, 3, 4];
for isol=1:length(solvers)
   solver = solvers{isol};
   fprintf('\nRunning convergence test with %s integrator (theoretical order %i)\n',...
        snames{isol},orders(isol))
   for i=1:ntests
      hs = Tf/100/2^(i-1);
      hf = hs/100;
      [t,Y,ns,nf] = solver(fs, ff, tout, Y0, hs, hf);
      err_max(i) = max(max(abs(Y'-Ytrue)));
      err_rms(i) = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
      fprintf('  Accuracy/Work Results (hs = %.3e, hf = %.3e):\n', hs, hf)
      if (i>1)
         fprintf('    maxerr = %.2e,  rmserr = %.2e (order = %.2e)\n',...
                 err_max(i), err_rms(i), log(err_rms(i)/err_rms(i-1))/log(0.5));
      else
         fprintf('    maxerr = %.2e,  rmserr = %.2e\n', err_max(i), err_rms(i));
      end
      fprintf('    slow steps = %i, fast steps = %i\n', ns, nf);
   end
end

% run with a adaptive RMIS+MIS method
fprintf('\nRunning with RMIS+MIS adaptive integrator\n')
[t,Y,ns,nf] = solve_RMIS_MIS_38(fs,ff,tout,Y0,rtol,atol,hmin,hmax,hmin,hmin);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n', err_max, err_rms);
fprintf('   slow steps = %i, fast steps = %i\n', ns, nf);


% end of script
