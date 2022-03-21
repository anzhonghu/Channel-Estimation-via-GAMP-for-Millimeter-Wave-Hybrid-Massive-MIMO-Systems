 function x=laplace(a,b,m,n)
u=rand(m,n)-0.5;
x=a-b*sign(u).*log(1-2*abs(u));
end   