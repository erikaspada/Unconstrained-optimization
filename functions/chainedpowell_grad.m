function grad = chainedpowell_grad(x)
n = length(x);
grad = zeros(n,1);
for i = 1:(n-2)/2
    grad(2*i-1) = grad(2*i-1) + 2*(x(2*i-1) + 10*x(2*i)) + 40*(x(2*i-1) - x(2*i + 2))^3;
    grad(2*i) = grad(2*i) + 20*(x(2*i-1) + 10*x(2*i)) + 4*(x(2*i) - 2*x(2*i+1))^3;

    grad(2*i+1) = grad(2*i+1) + 10*(x(2*i+1) - x(2*i+2)) -8*(x(2*i) - 2*x(2*i+1))^3;
    grad(2*i+2) = grad(2*i+2) -10*(x(2*i+1) - x(2*i+2)) - 40*(x(2*i-1) - x(2*i + 2))^3;
end
end
