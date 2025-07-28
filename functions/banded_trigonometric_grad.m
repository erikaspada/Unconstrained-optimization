function grad = banded_trigonometric_grad(x)
n = length(x);
grad = zeros(n,1);
for i = 1:n
    if i-1 == 0
        grad(i) = grad(i) + i*sin((x(i)));
        grad(i+1) = -i*cos((x(i+1)));
    elseif i+1 == n+1
        grad(i-1) = grad(i-1) + i*cos((x(i-1)));
        grad(i) = grad(i) + i*sin((x(i)));
    else
        grad(i-1) = grad(i-1) + i*cos((x(i-1)));
        grad(i) = grad(i) + i*sin((x(i)));
        grad(i+1) = -i*cos((x(i+1)));
    end
end
end
