function f = rosenbrock(x)
f = 0;
n = length(x);
for i=2:n
    f = f + 100*(x(i-1)^2-x(i))^2+(x(i-1)-1)^2;
end
end
