function grad = rosenbrock_grad(x)
n = length(x);
grad = zeros(n,1);
for i=2:n
    grad(i-1) = grad(i-1)+ 400*x(i-1)*(x(i-1)^2-x(i))+2*(x(i-1)-1);
    grad(i) = -200*(x(i-1)^2-x(i));
end
end
