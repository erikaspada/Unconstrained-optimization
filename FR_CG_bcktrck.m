function [xk, fk, gradfk_norm, k, xseq, btseq] = FR_CG_bcktrck(x0, f, gradf,alpha0, ...
kmax, tolgrad, c1, rho, btmax)
% INPUTS:
% x0 = n-dimensional column vector;
% f = function handle that describes a function R^n->R;
% gradf = function handle that describes the gradient of f;
% kmax = maximum number of iterations permitted;
% tolgrad = value used as stopping criterion w.r.t. the norm of the
% gradient;
% c1 = the factor of the Armijo condition that must be a scalar in (0,1);
% rho = fixed factor, lesser than 1, used for reducing alpha0;
% btmax = maximum number of steps for updating alpha during the
% backtracking strategy.
% Function handle for the armijo condition
farmijo = @(fk, alpha, gradfk, pk) ...
fk + c1 * alpha * gradfk' * pk;
% Initializations
xseq = zeros(length(x0), kmax);
btseq = zeros(1, kmax);
xk = x0;
fk = f(xk);
k = 0;
gradfk = gradf(xk);
gradfk_norm = norm(gradfk);
pk = -gradf(xk);
while k < kmax && gradfk_norm >= tolgrad
gradf_initial = gradf(xk);
% Reset the value of alpha
alpha = alpha0;
% Compute the candidate new xk
xnew = xk + alpha * pk;
% Compute the value of f in the candidate new xk
fnew = f(xnew);
bt = 0;
% Backtracking strategy:
% 2nd condition is the Armijo condition not satisfied
while bt < btmax && fnew > farmijo(fk, alpha, gradfk, pk)
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
% Find b_k+1 and the descent direction
beta = (gradfk_norm^2)/(norm(gradf_initial)^2);
pk = -gradfk + beta*pk;
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
