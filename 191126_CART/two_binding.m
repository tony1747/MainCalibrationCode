function dx = two_binding(t, x, param)
%param={a,b,p,d,g,alpha,beta,gamma,delta,epsilon}
global a b %n d p g m
E=x(1);
T=x(2);
dx(1,1)=param(1)*E*T/(param(3)+T)-param(5)/(E+param(4))*E*T-param(2)*E-param(6)/(E+param(4))*T*E*E;
dx(2,1)=a*T*(1-b*T)-param(7)/(E+param(4))*E*T-param(8)/(E+param(4))*T*E*E;