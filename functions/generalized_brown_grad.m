function grad = generalized_brown_grad(x)
n = length(x);
grad = zeros(n,1);
q = 0;
for i = 1:n/2
    q = q+(x(2*i-1)-3);
    grad(2*i-1) = (x(2*i -1)-3)/500 -1 + 20*exp(20*(x(2*i -1)-x(2*i)))+2*q;
    grad(2*i) = 1-20*exp(20*(x(2*i -1)-x(2*i)));
end
end
