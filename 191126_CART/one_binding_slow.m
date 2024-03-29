function dx=one_binding_slow(t, x, param)
global a b d %n d p g m
%param={k_1,alpha(=k_n1+k_3),beta(=k_n1+k_2),gamma(=k_n1+k_3+k_2),p,g}
%or calibrate k_n1param={k_1,k_2,k_3,k_n1,p,g}
%param={p,g,k_1_1,k_1_2,k_1_3,k_n1}

E=x(2);
T=x(1);
I=x(3);
dx(2,1)=param(1)*E*T/(param(2)+T)-param(3)*E*T-d*E+(param(6)+param(5))*I;
dx(1,1)=a*T*(1-b*T)-param(3)*E*T+(param(6)+param(4))*I;
dx(3,1)=param(3)*E*T-(param(6)+param(5)+param(4))*I;