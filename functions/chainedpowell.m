function f = chainedpowell(x)
n = length(x);
f = 0;
for i = 1:(n-2)/2
    f = f + (x(2*i-1) +10*x(2*i))^2 + 5*(x(2*i+1)-x(2*i+2))^2 + (x(2*i)-2*x(2*i+1))^4 + ...
        10*(x(2*i-1)- x(2*i + 2))^4;
end
end
