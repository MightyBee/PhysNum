function [a,erra, yFit] = fit(x,y)
X=[ones(length(x),1) x];
coef=X\y;
a=coef(2);
b=coef(1);

n=size(x,1);
k=0;
l=0;
for i = 1:n
  k = k + (y(i)-a*x(i)-b)^2;
  l = l + (x(i)-mean(x(1:n)))^2;
end
erra = sqrt(k/((n-2)*l));
yFit=X*coef;
end

