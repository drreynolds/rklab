function [q,p,qs,lq,lp,tol,Bs,As,Ls,BsE,AsE,LsE] = check_rk(B,reportL,doPlot,box,mname,fname)
% Usage: [q,p,qs,lq,lp,tol,Bs,As,Ls,BsE,AsE,LsE] = check_rk(B,reportL,doPlot,box,mname,fname)
%  or
% Usage: [q,p,qs,lq,lp,tol,Bs,As,Ls,BsE,AsE,LsE] = check_rk(B,reportL)
%  or
% Usage: [q,p,qs,lq,lp,tol,Bs,As,Ls,BsE,AsE,LsE] = check_rk(B)
%
% Checks the Butcher table B to determine:
% * the analytical order of accuracy (up to 6th) -> q
% * the embedding order of accuracy (up to 6th) -> p
% * the stage order of the method -> qs
% * the linear order of the method & embedding -> lq, lp
% * estimate B-stability of the method & embedding -> Bs and BsE
% * A-stability of the method & embedding -> As, AsE
% * L-stability of the method & embedding -> Ls, LsE
%
% When determining q, p, qs, lq and lq, we first declare 'success'
% at a very loose tolerance of 1e-8, and then tighten this until
% any of these decrease.  The tightest tested tolerance where these
% outputs show full order is returned as 'tol'.
%
% It is assumed that B has block structure
%     B = [c, A; 0, b]
% for a standard Runge-Kutta method, or
%     B = [c, A; 0, b; 0, b2]
% if the method has an embedded error indicator.
%
% If the method has no embedding, then we set p=lp=BsE=AsE=LsE=0.
%
% We check for "linear order of accuracy" on the autonomous ODE
%    y' = A*y,
% and for actual order of accuracy, using analytical order
% conditions from Sandu & Gunther, SINUM 53, 2015.
%
% We compute the 'stage order' of the method through checking
% the Butcher simplifying assumption C(q) as prescribed in the article
% J.C. Butcher, "Implicit Runge-Kutta processes", Math Comp 18 (1964), pp 50-64.
%
% For B-stability, we check for positive semi-definiteness of
%    M(i,j) = [b(i)*A(i,j) + b(j)*A(j,i) - b(i)*b(j)]
% as outlined in
%    K. Burrage and J.C. Butcher, "Stability criteria for implicit
%    Runge-Kutta methods", SIAM J. Numer. Anal. 16 (1979), pp
%    46-57.
%
% For A-stability we do ...
%
% For L-stability we do ...
%
% Inputs:
%   B -- RK table
%   reportL -- integer flag denoting print level:
%              >1 -> all results
%              =1 -> final results
%            else -> no results
%   doPlot -- boolean flag denoting whether to create/save plot
%   box -- [xl,xr,yl,yr] sub-region of the complex plane to
%          use in plot (ignored if doPlot is false)
%   mname -- string containing method name to insert into plot title
%   fname -- string containing method name to insert into plot filenames
%
% Outputs are described above.
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2018, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

% handle different call structures
if ((nargin ~= 6) && (nargin ~= 2) && (nargin ~= 1))
   error('check_rk error: must be called with exactly 1, 2 or 6 arguments');
end
if (nargin < 6)   % set defaults for doPlot,box,mname,fname
   doPlot = false;
   mname = 'Butcher table';
end
if (nargin == 1)
   reportL = 0;
end

% set tolerance for assessing stability
StabTol = 1e-8;

% set initial tolerance on method/embedding order
tol = 1e-8;

% extract components of Butcher table
[m,n] = size(B);
if (m == n)        % no embedding
   embedded = 0;
   d = 0;
elseif (m == n+1)  % has an embedding
   embedded = 1;
   d = B(m,2:n)';
else   % illegal input
   error('illegal Butcher table input')
end
s = n-1;
c = B(1:s,1);
b = B(s+1,2:n)';
A = B(1:s,2:n);

% check stage order of method
qs = stage_order(c,A,tol);

% determine whether method satisfies the 'stiffly accurate' property
SA = 1;
if (norm(b'-A(end,:),inf) > tol)
   SA = 0;
end

% assess order & stability of method
[q,lq] = table_order(c,A,b,tol,reportL);
[As,Bs,Ls] = stability(A,b,StabTol);

% report on method
if (reportL>0)
   fprintf('  %s stage order = %i\n', mname, qs);
   fprintf('    method:    q = %i,  lq = %i,  As = %i,  Bs = %i,  Ls = %i,  SA = %i\n', ...
           q, lq, As, Bs, Ls, SA);
end

% if an embedding exists, assess embedding properties
if (embedded)
   [p,lp] = table_order(c,A,d,tol,reportL);
   [AsE,BsE,LsE] = stability(A,d,StabTol);
   if (reportL>0)
      fprintf('    embedding: p = %i,  lp = %i,  As = %i,  Bs = %i,  Ls = %i\n', ...
              p, lp, AsE, BsE, LsE);
   end
else
   p = 0;
   lp = 0;
   AsE = 0;
   BsE = 0;
   LsE = 0;
end


% hone in on tightest tolerance where method and embedding order
% are retained
if (q > 0)
   for testtol = 0.1.^(8:40)
      [q_,lq_] = table_order(c,A,b,testtol,0);
      qs_ = stage_order(c,A,tol);
      if ((q_ < q) || (lq_ < lq) || (qs_ < qs))
         break
      end
      if (embedded)
         [p_,lp_] = table_order(c,A,d,testtol,0);
         if ((p_ < p) || (lp_ < lp))
            break
         end
      end
      tol = testtol;
   end
else
  tol = 1;
end


% generate plot of stability region
if (doPlot)
   figure()
   xl = box(1:2);  yl = box(3:4);
   xax = plot(linspace(xl(1),xl(2),10),zeros(1,10),'k:'); hold on
   yax = plot(zeros(1,10),linspace(yl(1),yl(2),10),'k:');
   [X,Y] = stab_region(A,b,box);
   plot(X,Y,'r-')
   if (embedded)
      [X,Y] = stab_region(A,d,box);
      plot(X,Y,'b--')
   end
   set(get(get(xax,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
   set(get(get(yax,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
   axis(box)
   xlabel('Re(z)')
   ylabel('Im(z)')
   if (embedded)
      title(sprintf('%s stability regions, order %i',mname,q))
      legend('method','embedding')
   else
      title(sprintf('%s stability region, order %i',mname,q))
   end
   print(sprintf('%s_stab_region.png', fname), '-dpng');
   print(sprintf('%s_stab_region.eps', fname), '-depsc');
   savefig(sprintf('%s_stab_region.fig', fname));
end

% end check_rk function
end



%---------------------- Utility routines ----------------------%


function [q,lq] = table_order(c,A,b,tol,reportL)
%
% Utility routine to assess analytical order of accuracy for
% (c,A,b) RK table

   % initialize failure flags, order of accuracy
   Ofail = false;
   Lfail = false;
   q = -1;
   lq = -1;

   % get number of stages
   s = length(b);

   % convert b to column vector for these tests
   b = reshape(b,s,1);

   % create vector of ones
   e = ones(s,1);

   % check row sum condition
   for i=1:s
      tst = c(i) - sum(A(i,:));
      if (abs(tst) > tol)
         Ofail = true;
         Lfail = true;
         if (reportL>1)
            fprintf('    Method fails row sum condition, i = %i, tst = %g\n', i, tst);
         end
      end
   end
   if (~Ofail)
      q = 0;
      lq = 0;
   end

   % check for first order
   if (~Ofail || ~Lfail)
      tst = double(b'*e - sym(1));
      if (abs(tst) > tol)
         Ofail = true;
         Lfail = true;
         if (reportL>1)
            fprintf('    Method fails 1st order condition (tst = %g)\n', tst);
         end
      end
      if (reportL>1)
         if (~Ofail), fprintf('  Method passes order 1 condition\n'); end
      end
      if (~Ofail)
         q = 1;
         lq = 1;
      end
   end

   % check for second order
   if (~Ofail || ~Lfail)
      tst = double(b'*c - sym(0.5));
      if (abs(tst) > tol)
         Ofail = true;
         Lfail = true;
         if (reportL>1)
            fprintf('    Method fails 2nd order condition (tst = %g)\n', tst);
         end
      end
      if (reportL>1)
         if (~Ofail),  fprintf('  Method passes order 2 condition\n'); end
      end
      if (~Ofail)
         q = 2;
         lq = 2;
      end
   end

   % check for third order
   if (~Ofail || ~Lfail)

      tst = double(b'*(c.*c) - sym(1)/sym(3));
      if (abs(tst) > tol)
         Ofail = true;
         if (reportL>1)
            fprintf('    Method fails 3rd order condition (tst = %g)\n', tst);
         end
      end

      tst = double(b'*(A*c) - sym(1)/sym(6));
      if (abs(tst) > tol)
         Ofail = true;
         Lfail = true;
         if (reportL>1)
            fprintf('    Method fails linear 3rd order condition (tst = %g)\n', tst);
         end
      end

      if (reportL>1)
         if (~Ofail),  fprintf('  Method passes order 3 conditions\n'); end
         if (~Lfail && Ofail),  fprintf('  Method passes linear order 3 condition\n'); end
      end
      if (~Ofail)
         q = 3;
      end
      if (~Lfail)
         lq = 3;
      end
   end

   % check for fourth order
   if (~Ofail || ~Lfail)

      if (~Ofail)
         tst = double(b'*(c.^3) - sym(1)/sym(4));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 4th order condition A (tst = %g)\n', tst);
            end
         end

         tst = double((b.*c)'*(A*c) - sym(1)/sym(8));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 4th order condition B (tst = %g)\n', tst);
            end
         end

         tst = double(b'*A*(c.*c) - sym(1)/sym(12));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 4th order condition C (tst = %g)\n', tst);
            end
         end
      end

      tst = double(b'*A*A*c - sym(1)/sym(24));
      if (abs(tst) > tol)
         Ofail = true;
         Lfail = true;
         if (reportL>1)
            fprintf('    Method fails linear 4th order condition (tst = %g)\n', tst);
         end
      end

      if (reportL>1)
         if (~Ofail),  fprintf('  Method passes order 4 conditions\n'); end
         if (~Lfail && Ofail),  fprintf('  Method passes linear order 4 condition\n'); end
      end
      if (~Ofail)
         q = 4;
      end
      if (~Lfail)
         lq = 4;
      end
   end

   % check for fifth order
   if (~Ofail || ~Lfail)

      if (~Ofail)

         tst = double(b'*(c.^4) - sym(1)/sym(5));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 5th order condition A (tst = %g)\n', tst);
            end
         end

         tst = double((b.*c.*c)'*(A*c) - sym(1)/sym(10));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 5th order condition B (tst = %g)\n', tst);
            end
         end

         tst = double(b'*((A*c).^2) - sym(1)/sym(20));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 5th order condition C (tst = %g)\n', tst);
            end
         end

         tst = double((b.*c)'*A*(c.*c) - sym(1)/sym(15));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 5th order condition D (tst = %g)\n', tst);
            end
         end

         tst = double(b'*A*(c.^3) - sym(1)/sym(20));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 5th order condition E (tst = %g)\n', tst);
            end
         end

         tst = double((b.*c)'*A*A*c - sym(1)/sym(30));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 5th order condition F (tst = %g)\n', tst);
            end
         end

         tst = double(b'*A*(c.*(A*c)) - sym(1)/sym(40));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 5th order condition G (tst = %g)\n', tst);
            end
         end

         tst = double(b'*A*A*(c.*c) - sym(1)/sym(60));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 5th order condition H (tst = %g)\n', tst);
            end
         end

      end

      tst = double(b'*A*A*A*c - sym(1)/sym(120));
      if (abs(tst) > tol)
         Ofail = true;
         Lfail = true;
         if (reportL>1)
            fprintf('    Method fails linear 5th order condition I (tst = %g)\n', tst);
         end
      end

      if (reportL>1)
         if (~Ofail),  fprintf('  Method passes order 5 conditions\n'); end
         if (~Lfail && Ofail),  fprintf('  Method passes linear order 5 condition\n'); end
      end
      if (~Ofail)
         q = 5;
      end
      if (~Lfail)
         lq = 5;
      end
   end

   % check for sixth order
   if (~Ofail || ~Lfail)

      if (~Ofail)

         tst = double(b'*(c.^5) - sym(1)/sym(6));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 6th order condition A (tst = %g)\n', tst);
            end
         end

         tst = double((b.*c.^3)'*(A*c) - sym(1)/sym(12));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 6th order condition B (tst = %g)\n', tst);
            end
         end

         tst = double(b'*(c.*(A*c).^2) - sym(1)/sym(24));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 6th order condition C (tst = %g)\n', tst);
            end
         end

         tst = double((b.*(c.*c))'*A*(c.*c) - sym(1)/sym(18));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 6th order condition D (tst = %g)\n', tst);
            end
         end

         tst = double((b.*(c.*c))'*A*A*c - sym(1)/sym(36));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 6th order condition E (tst = %g)\n', tst);
            end
         end

         tst = double(b'*((A*A*c).*(A*c)) - sym(1)/sym(72));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 6th order condition F (tst = %g)\n', tst);
            end
         end

         tst = double(b'*(c.*(A*(c.^3))) - sym(1)/sym(24));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 6th order condition G (tst = %g)\n', tst);
            end
         end

         tst = double(b'*(c.*(A*(c.*(A*c)))) - sym(1)/sym(48));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 6th order condition H (tst = %g)\n', tst);
            end
         end

         tst = double( b'*(c.*(A*A*(c.*c))) - sym(1)/sym(72));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 6th order condition I (tst = %g)\n', tst);
            end
         end

         tst = double( b'*(c.*(A*A*A*c)) - sym(1)/sym(144));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 6th order condition J (tst = %g)\n', tst);
            end
         end

         tst = double( b'*A*(c.^4) - sym(1)/sym(30));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 6th order condition K (tst = %g)\n', tst);
            end
         end

         tst = double( b'*A*(c.*c.*(A*c)) - sym(1)/sym(60));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 6th order condition L (tst = %g)\n', tst);
            end
         end

         tst = double( b'*A*((A*c).^2) - sym(1)/sym(120));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 6th order condition M (tst = %g)\n', tst);
            end
         end

         tst = double( b'*A*(c.*(A*(c.*c))) - sym(1)/sym(90));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 6th order condition N (tst = %g)\n', tst);
            end
         end

         tst = double( b'*A*(c.*(A*A*c)) - sym(1)/sym(180));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 6th order condition O (tst = %g)\n', tst);
            end
         end

         tst = double( b'*A*A*(c.^3) - sym(1)/sym(120));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 6th order condition P (tst = %g)\n', tst);
            end
         end

         tst = double( b'*A*A*(c.*(A*c)) - sym(1)/sym(240));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 6th order condition Q (tst = %g)\n', tst);
            end
         end

         tst = double( b'*A*A*A*(c.*c) - sym(1)/sym(360));
         if (abs(tst) > tol)
            Ofail = true;
            if (reportL>1)
               fprintf('    Method fails 6th order condition R (tst = %g)\n', tst);
            end
         end

      end

      tst = double(b'*A*A*A*A*c - sym(1)/sym(720));
      if (abs(tst) > tol)
         Ofail = true;
         Lfail = true;
         if (reportL>1)
            fprintf('    Method fails linear 6th order condition S (tst = %g)\n', tst);
         end
      end

      if (reportL>1)
         if (~Ofail),  fprintf('  Method passes order 6 conditions\n'); end
         if (~Lfail && Ofail),  fprintf('  Method passes linear order 6 condition\n'); end
      end
      if (~Ofail)
         q = 6;
      end
      if (~Lfail)
         lq = 6;
      end
   end

   % check for linear seventh order
   if (~Lfail)
      tst = double(b'*A*A*A*A*A*c - sym(1)/sym(factorial(7)));
      if (abs(tst) > tol)
         Lfail = true;
         if (reportL>1)
            fprintf('    Method fails linear 7th order condition (tst = %g)\n', tst);
         end
      end

      if (reportL>1)
         if (~Lfail),  fprintf('  Method passes linear order 7 condition\n'); end
      end
      if (~Lfail)
         lq = 7;
      end
   end

   % check for linear eighth order
   if (~Lfail)
      tst = double(b'*A*A*A*A*A*A*c - sym(1)/sym(factorial(8)));
      if (abs(tst) > tol)
         Lfail = true;
         if (reportL>1)
            fprintf('    Method fails linear 8th order condition (tst = %g)\n', tst);
         end
      end

      if (reportL>1)
         if (~Lfail),  fprintf('  Method passes linear order 8 condition\n'); end
      end
      if (~Lfail)
         lq = 8;
      end
   end


   %--- if the method seems to have order 6, continue with Butcher's simplifying assumptions ---%
   if (q == 6)

      if (reportL>1)
         fprintf('  Using Butcher''s simplifying assumptions to continue tests\n')
      end

      % determine P, Q and R for the Butcher simplifying assumptions
      %   B(P):
      P = 0;
      for i=1:1000
         LHS = b'*(c.^(i-1));
         RHS = sym(1)/(i);
         tst = double(RHS-LHS);
         if (abs(tst)>tol)
            break;
         end
         P = P+1;
      end

      %   C(Q):
      Q = 0;
      for k=1:1000
         alltrue = 1;
         for i=1:s
            LHS = A(i,:)*(c.^(k-1));
            RHS = sym(c(i)^k)/sym(k);
            tst = double(RHS-LHS);
            if (abs(tst)>tol)
               alltrue=0;
               break;
            end
         end
         if (alltrue == 1)
            Q = Q+1;
         else
            break;
         end
      end

      %   D(R):
      R = 0;
      for k=1:1000
         alltrue = 1;
         for j=1:s
            LHS = 0;
            for i=1:s
               LHS = LHS + A(i,j)*b(i)*c(i)^(k-1);
            end
            RHS = b(j)/sym(k)*(sym(1)-c(j)^k);
            tst = double(RHS-LHS);
            if (abs(tst)>tol)
               alltrue = 0;
               break;
            end
         end
         if (alltrue == 1)
            R = R+1;
         else
            break;
         end
      end

      % determine q
      q2 = 0;
      for i=1:P
         if ((q2 > Q+R+1) || (q2 > 2*Q+2))
            break;
         end
         q2 = q2+1;
      end
      q = max(q,q2);
      if (reportL>1)
         fprintf('    passes order %i conditions\n',q);
      end

   end


% end table_order function
end


function [qs] = stage_order(c,A,tol)
%
% Applies the Butcher simplifying assumption C(q) as prescribed in
%
% J.C. Butcher, "Implicit Runge-Kutta processes", Math Comp
%    18 (1964), pp 50-64.
%
% to determine the stage order(s) of a Runge Kutta method.
% Specifically, the stage order of a RK is based on the
% Butcher table components c and A.

   % check for legal input tolerance
   if (tol <= 0.0)
      error('illegal input tolerance, must be > 0')
   end

   % verify compatability of Butcher table inputs
   if ( (length(c) ~= size(A,1)) )
      error('incompatible Butcher table inputs')
   end

   % determine number of stages
   s = length(c);

   % perform test for RK method
   qs = 0;
   for k=1:s
      alltrue = 1;
      for i=1:s
         LHS = A(i,:)*(c.^(k-1));
         RHS = c(i)^k/sym(k);
         tst = double(RHS-LHS);
         if (abs(tst)>tol)
            alltrue=0;
            break;
         end
      end
      if (alltrue == 1)
         qs = qs+1;
      else
         break;
      end
   end

% end of function
end



function [As,Bs,Ls] = stability(A,b,tol)
%
% estimates A, B and L stability of a RK method

   % estimate B stability
   s = length(b);
   M = zeros(s,s);
   for j=1:s
      for i=1:s
         M(i,j) = b(i)*A(i,j) + b(j)*A(j,i) - b(i)*b(j);
      end
   end
   lam = eig(M);
   Bs = 1;
   for i=1:s
      if (lam(i) < -tol)
         Bs = 0;
         break;
      end
   end


   % estimate A stability:
   %   Hairer & Wanner, Solving ODEs II: an RK method is A-stable iff
   %      (i)  |R(iy)| < 1  for all y in R
   %      (ii) R(z) is analytic in the left half-plane

   % generate stability function coefficients (reverse ordering)
   [alpha,beta] = stab_function(A,b);
   alpha_dbl = double(alpha);
   beta_dbl = double(beta);

   % remove zero coefficients for highest-order terms
   for i=1:length(alpha_dbl)-1
      if (abs(alpha_dbl(end)) < eps)
         alpha_dbl = alpha_dbl(1:end-1);
      else
         break;
      end
   end
   for i=1:length(beta_dbl)-1
      if (abs(beta_dbl(end)) < eps)
         beta_dbl = beta_dbl(1:end-1);
      else
         break;
      end
   end

   % check analytic in left half-plane
   beta2 = fliplr(beta_dbl);
   alpha2 = fliplr(alpha_dbl);
   rt = roots(beta2);
   As = 1;
   for i=1:length(rt)
      if (real(rt(i)) < -tol)
         As=0;
         break
      end
   end

   % check along imaginary axis
   ztests = [0, sqrt(-1)*logspace(-5,4,10000)];
   for j = 1:length(ztests)
      z = ztests(j);
      if (abs(polyval(alpha2,z)/polyval(beta2,z))>1+tol)
         As = 0;
         break;
      end
   end


   % check L stability:
   %   if degree of denominator is greater than degree of numerator in
   %   ational stability polynomial, then it will be L-stable
   numdeg = length(alpha)-1;
   for i=0:length(alpha)-1
      if (abs(alpha(end-i)) < tol^2)
         numdeg = numdeg-1;
      end
   end
   dendeg = length(beta)-1;
   for i=0:length(beta)-1
      if (abs(beta(end-i)) < tol^2)
         dendeg = dendeg-1;
      end
   end
   if (dendeg > numdeg)
      Ls = 1;
   else
      Ls = 0;
   end

end
