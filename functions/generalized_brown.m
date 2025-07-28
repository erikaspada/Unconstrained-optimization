function f = generalized_brown(x)
n = length(x);
f = 0;
k=0;
for i = 1:n/2
    f = f + ((x(2*i-1)-3)^2)/1000 -(x(2*i-1)-x(2*i))+ exp(20*(x(2*i-1)-x(2*i)));
end
for j = 1:n/2
    k= k + x(2*j-1)-3;
end
f= f+ k^2;
end
