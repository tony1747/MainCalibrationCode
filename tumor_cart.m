function dx = tumor_cart(t, x, param)
%param={a,b,n,s,p,g,m,d}
E=x(1);
T=x(2);
dx(1,1)=param(4)*E*T/(param(5)+T)-param(6)*E*T-param(7)*E;
dx(2,1)=param(1)*T*(1-param(2)*T)-param(3)*E*T;