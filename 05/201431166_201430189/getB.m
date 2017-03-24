function b = getB(t,i)
t0 = 0;
tf = 10;
u = (t-t0)/(tf-t0);
b = nchoosek(5,i)*((1-u)^(5-i))*(u^i);
end