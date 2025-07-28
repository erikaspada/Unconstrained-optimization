function [xk, fk, k] = nelder_mead(x0, f, tol)
% [xk, fk] = nelder_mead(x0, f, tol)
% Function that performs the Nelder Mead method
% INPUTS:
% x0 = n-dimensional column vector;
% f = function handle that describes a function R^n->R;
% tol = tolerance for stopping criterion;
% OUTPUTS:
% xk = the last x computed by the function;
% fk = the value f(xk);
%Definition of parameters
rho = 1; %rho > 0 (tipically 1);
chi = 2; %chi > 1 (tipically 2);
gamma = 0.5; %0 < gamma < 1 (tipically 0.5);
sigma = 0.5; %0 < sigma < 1 (tipically 0.5);
k=0;
% Create initial simplex
n = length(x0);
X = zeros(n, n+1);
X(:,1) = x0;
for i = 2:n+1
X(:,i) = x0;
X(i-1,i) = x0(i-1) + 1;
end
err=1;
while err>tol
k=k+1;
% Sort simplex by function value (ordering phase)
f_vect = zeros(1,n+1);
for i = 1:n+1
f_vect(i) = f(X(:,i));
end
[f_vect,idx] = sort(f_vect);
X = X(:,idx);
% Compute centroid
x_bar = mean(X(:,1:n),2);
% Reflection phase
x_r = x_bar + rho*(x_bar - X(:,n+1));
f_r = f(x_r);
if f_vect(1) <= f_r && f_r < f_vect(n)
X(:,n+1) = x_r;
continue; %Returns to the beginning of the while (new simplex)
end
% Expansion phase
if f_r < f_vect(1)
x_e = x_bar + chi*(x_r - x_bar);
f_e = f(x_e);
if f_e < f_r
X(:,n+1) = x_e;
else
X(:,n+1) = x_r;
end
continue;
end
% Contraction phase
if f_r >= f_vect(n)
x_c = x_bar - gamma*(x_bar - X(:,n+1));
f_c = f(x_c);
if f_c < f_vect(n+1)
X(:,n+1) = x_c;
else
% Shrinking phase
for i = 2:n+1
X(:,i) = X(:,1) + sigma*(X(:,i) - X(:,1));
end
end
end
% Stopping criterion
err=norm(f_vect(1) - f_vect(n+1));
end
% Output
xk = X(:,1);
fk = f_vect(1);
end
