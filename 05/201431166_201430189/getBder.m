function b = getBder(t,i)
t0 = 0;
tf = 10;
u = (t-t0)/(tf-t0);
if i == 0
    b = -5*((1-u)^4);
elseif i == 5
    b = 5*(u^4);
else
    b = nchoosek(5,i)*((1-u)^(4-i))*(u^(i-1))*(i-5*u);
end
b = b/(tf-t0);    
end