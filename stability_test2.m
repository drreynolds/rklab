% stability_test2 -- driver for simple ODE system
%      u' = lamda*u + cos(t) - lamda*sin(t)
%      v' = gamma*v + 100*cos(100*t) - gamma*sin(100*t)
% for t in the interval [0.0, 10.0], with initial conditions:
% u(0)=v(0)=0, with parameters 
%     lamda = -1
%     gamma = -100
% Based on these parameters, as well as the linear stability
% regions for both the KW3 and 3/8-Rule ERK tables, we should 
% expect that the methods have maximum stable step sizes:
%    KW3: 
%         h*lambda = -2.5  =>  h = 2.5
%         hfast*gamma = -2.5  =>  hfast = 0.025
%    3/8-Rule: 
%         h*lambda = -2.75  =>  h = 2.75
%         hfast*gamma = -2.75  =>  hfast = 0.0275
% The tests are set up to start with those values, and subsequently
% shrink the value of "h".  This will eventually require that hfast
% also begin to shrink (once h gets too small).
%
% Daniel Reynolds
% Mathematics @ SMU
% August 2018               

clear, close all

% Set problem parameters
lamda = -1;
gamma = -100;
Ti    = 0;
Tf    = 10;
Y0    = [0;0];
fs    = @(t,y) [lamda*y(1) + cos(t) - lamda*sin(t); 0];    % slow/fast RHS functions
ff    = @(t,y) [0; gamma*y(2) + 100*cos(100*t) - gamma*sin(100*t)];
n     = 5;       % this guarantees hmax = 2.5, which should be stable for both methods
tout  = linspace(Ti,Tf,n);           % intermediate times for solution
hmax  = tout(2)-tout(1);
hfast = 0.025;
ytrue = @(t) [sin(t); sin(100*t)];
Ytrue = ytrue(tout);
ntest = 10;

% loop over MIS-KW3 tests, systematically decreasing h at each test
fprintf('\nMIS-KW3 integrator:\n')
for k=0:ntest-1
   
   % Set up time steps, m value
   h = hmax/(2^k);
   m = ceil(h/hfast);

   % MIS-KW3 test
   fprintf('  h = %g,  hfast (goal) = %g,',h,hfast)
   [t,Yout,nfslow,nffast] = solve_MIS_KW3(fs, ff, tout, Y0, h, hfast);
   err_max = max(max(abs(Yout - Ytrue)));
   fprintf('  nffast = %5i,  nfslow = %4i, maxerr = %.3e\n',nffast,nfslow,err_max);
   
end

   
% loop over RMIS-3/8 tests, systematically decreasing h at each test
fprintf('\nRMIS-3/8 integrator:\n')
for k=0:ntest-1
   
   % Set up time steps, m value
   h = hmax/(2^k);
   m = ceil(h/hfast);

   % RMIS-3/8 test
   fprintf('  h = %g,  hfast (goal) = %g,',h,hfast)
   [t,Yout,nfslow,nffast] = solve_RMIS_38(fs, ff, tout, Y0, h, hfast);
   err_max = max(max(abs(Yout - Ytrue)));
   fprintf('  nffast = %5i,  nfslow = %4i, maxerr = %.3e\n',nffast,nfslow,err_max);
   
end
