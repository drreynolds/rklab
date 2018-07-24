% multirate driver for stiff brusselator test problem:
%      [u]'   [ a - (w+1)*u + u^2*v ]   [     0    ]
%      [v]' = [     w*u - u^2*v     ] + [     0    ] = fslow(u,v,w) + ffast(u,v,w)
%      [w]'   [         -u*w        ]   [ (b-w)/ep ]
% where u(0) = 1.2, v(0) = 3.1 and w(0) = 3, with prameters a=1,
% b=3.5 and ep=5e-6.  We evaluate over the time interval [0,10].
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% July 2018
% All Rights Reserved
clear

% set problem parameters
a = 1;
b = 3.5;
ep = 1e-3;
fs  = @(t,y) [a - (y(3)+1)*y(1) + y(1)*y(1)*y(2); y(3)*y(1) - y(1)*y(1)*y(2); -y(3)*y(1)];
ff  = @(t,y) [0; 0; (b-y(3))/ep];
fn  = @(t,y) fs(t,y) + ff(t,y);
Jff = @(t,y) [0, 0, 0; 0, 0, 0; 0, 0, -1/ep ];
Tf = 10;
tout = linspace(0,Tf,100);
hmin = 1e-7;
hmax = 1.0;
rtol = 1e-3;
atol = 1e-14*ones(3,1);
u0 = 1.2;
v0 = 3.1;
w0 = 3;
Y0 = [u0; v0; w0];
ntests = 7;
test_adaptive = false;

% construct reference solution
opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s(fn, tout, Y0, opts);
figure()
plot(tout,Ytrue)
xlabel('t','FontSize',12), ylabel('y','FontSize',12)
title('Brusselator ODE test','FontSize',14)
set(gca,'FontSize',12)
print('-dpng','brusselator')


% perform convergence tests with fixed-step multirate solvers
solvers = {@solve_MIS_KW3, @solve_MIS_38, @solve_RMIS_KW3, @solve_RMIS_38};
snames = {'MIS-KW3', 'MIS-3/8', 'RMIS-KW3', 'RMIS-3/8'};
orders = [3, 3, 3, 4];
for isol=1:length(solvers)
   solver = solvers{isol};
   fprintf('\nRunning convergence test with %s integrator (theoretical order %i)\n',...
        snames{isol},orders(isol))
   for i=1:ntests
      hs = Tf/100/2^(i+1);
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

if (test_adaptive)
   % run with a adaptive RMIS+MIS method
   fprintf('\nRunning with RMIS+MIS adaptive integrator\n')
   [t,Y,ns,nf] = solve_RMIS_MIS_38(fs,ff,tout,Y0,rtol,atol,hmin,hmax,hmin,hmin);
   err_max = max(max(abs(Y'-Ytrue)));
   err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
   fprintf('Accuracy/Work Results:\n')
   fprintf('   maxerr = %.5e,  rmserr = %.5e\n', err_max, err_rms);
   fprintf('   slow steps = %i, fast steps = %i\n', ns, nf);
end

% end of script
