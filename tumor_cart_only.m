function dx = tumor_cart_only(t, x, param)
%param={a,b,n,s,p,g,m,d}
global a b %n d p g m
E=x(1);
T=x(2);
dx(1,1)=param(2)*E*T/(param(3)+T)-param(4)*E*T-param(5)*E;
dx(2,1)=a*T*(1-b*T)-param(1)*E*T;