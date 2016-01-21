N = 100;
r = zeros(10000,1);

for pp = 0.0001:0.0001:1;
   r(round(pp*10000),1) = N*pp*(1-pp)^(N-1); 
end
sum(r)/10000
figure, plot(r);

fun = @(x) N.*x.*(1-x).^(N-1);
q = integral(fun,0,1)