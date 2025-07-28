function [xk, fk, gradfk_norm, k, xseq, btseq] = ...
innewton_general(x0, f, gradf, Hessf, ...
kmax, tolgrad, c1, rho, btmax, FDgrad, FDHess, h, fterms, pcg_maxit)
% INPUTS:
% x0 = n-dimensional column vector;
% f = function handle that describes a function R^n->R;
% gradf = function handle that describes the gradient of f (not necessarily
% used);
% Hessf = function handle that describes the Hessian of f (not necessarily
% used);
% kmax = maximum number of iterations permitted;
% tolgrad = value used as stopping criterion w.r.t. the norm of the
% gradient;
% c1 = the factor of the Armijo condition that must be a scalar in (0,1);
% rho = fixed factor, lesser than 1, used for reducing alpha0;
% btmax = maximum number of steps for updating alpha during the
% backtracking strategy;
% FDgrad = 'fw' (FD Forward approx. for gradf), 'c' (FD Centered approx.
% for gradf), any other string (usage of input Hessf)
% FDHess = 'fw' (FD approx. for Hessf), 'Jfw' (Jacobian FD Forward
% approx. of Hessf), 'Jc' (Jacobian FD Centered approx. of Hessf), 'MF'
% (Matrix Free implementation for solving Hessf(xk)pk=-gradf(xk)), any
% other string (usage of input Hessf);
% h = approximation step for FD (if used);
% fterms = f. handle "@(gradfk, k) ..." that returns the forcing term
% eta_k at each iteration
% pcg_maxit = maximum number of iterations for the pcg solver.
%%
OUTPUTS:
% xk = the last x computed by the function;
% fk = the value f(xk);
% gradfk_norm = value of the norm of gradf(xk)
% k = index of the last iteration performed
% xseq = n-by-k matrix where the columns are the xk computed during the
% iterations
% btseq = 1-by-k vector where elements are the number of backtracking
% iterations at each optimization step.
%
switch FDgrad
case 'fw'
% Calculating the gradient with forward finite differences
gradf = @(x) findiff_grad(f, x, h, 'fw');
case 'c'
%Calculating the gradient with central finite differences
gradf = @(x) findiff_grad(f, x, h, 'c');
otherwise
% Use the exact gradient as in inpunt gradf
end

if isequal(FDgrad, 'fw') || isequal(FDgrad, 'c')
switch FDHess
case 'fw'
% OVERWRITE Hessf WITH A F. HANDLE THAT USES findiff_Hess
Hessf = @(x) findiff_Hess(f, x, sqrt(h));
case 'MF'
% DEFINE a f. handle for the product of Hessf * p USING THE
% GRADIENT
Hessf_pk = @(x, p) (gradf(x + h * p) - gradf(x)) / h;
otherwise
% Use the Hessf as in input
end
else
switch FDHess
case 'fw'
% OVERWRITE Hessf WITH A F. HANDLE THAT USES findiff_Hess
Hessf = @(x) findiff_Hess(f, x, sqrt(h));
case 'Jfw'
% OVERWRITE Hessf WITH A F. HANDLE THAT USES findiff_J
% (with option 'fw')
Hessf = @(x) findiff_J(gradf, x, h, 'fw');
case 'Jc'
% OVERWRITE Hessf WITH A F. HANDLE THAT USES findiff_J
% (with option 'c')
Hessf = @(x) findiff_J(gradf, x, h, 'c');
case 'MF'
% DEFINE a f. handle for the product of Hessf * p USING THE
% GRADIENT
Hessf_pk = @(x, p) (gradf(x + h * p) - gradf(x)) / h;
otherwise
% Use the Hessf as in input
end
end
% Function handle for the armijo condition
farmijo = @(fk, alpha, gradfk, pk) ...
fk + c1 * alpha * gradfk' * pk;
n = length(x0);
% Initializations
xseq = zeros(n, kmax);
btseq = zeros(1, kmax);
xk = x0;
fk = f(xk);
k = 0;
gradfk = gradf(xk);
gradfk_norm = norm(gradfk);
while k < kmax && gradfk_norm >= tolgrad
% "INEXACTLY" compute the descent direction as solution of
% Hessf(xk) p = - gradf(xk)
% TOLERANCE VARYING W.R.T. FORCING TERMS:
eta_k = fterms(gradfk, k);
% We will use directly eta_k as tolerance in the pcg because
% this function looks at the RELATIVE RESIDUAL and not the RESIDUAL
switch FDHess
case 'MF'
pk = pcg(Hessfk_pk, -gradfk, eta_k, pcg_maxit);
otherwise
% DIRECT METHOD
% pk = -Hessf(xk)\gradf(xk)
% ITERATIVE METHOD
pk = pcg(Hessf(xk), -gradfk, eta_k, pcg_maxit);
end
% Reset the value of alpha
alpha = 1;
% Compute the candidate new xk
xnew = xk + alpha * pk;
% Compute the value of f in the candidate new xk
fnew = f(xnew);
bt = 0;
% Backtracking strategy:
% 2nd condition is the Armijo condition not satisfied
while bt < btmax && fnew > farmijo(fk, alpha, xk, pk)
% Reduce the value of alpha
alpha = rho * alpha;
% Update xnew and fnew w.r.t. the reduced alpha
xnew = xk + alpha * pk;
fnew = f(xnew);
% Increase the counter by one
bt = bt + 1;
end
% Update xk, fk, gradfk_norm
xk = xnew;
fk = fnew;
gradfk = gradf(xk);
gradfk_norm = norm(gradfk);
% Increase the step by one
k = k + 1;
% Store current xk in xseq
xseq(:, k) = xk;
% Store bt iterations in btseq
btseq(k) = bt;
end
% "Cut" xseq and btseq to the correct size
xseq = xseq(:, 1:k);
btseq = btseq(1:k);
end
