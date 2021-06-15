function dx = two_binding_slow(t, x, param)
global a b d %n d p g m
%param={p,g,k_1,alpha(=k_n1+k_3),beta(=k_n1+k_2),gamma(=k_n1+k_3+k_2)}
%param={p,g,k_1_1,k_1_2,k_1_3,k_2_1,k_2_2,k_2_3,k_2_4}
%k_1_1 = k_2_1
%k_1_3 = k_2_4
%k_1_2 = k_2_2+k_2_3 or k_2_2=k_2_3= k_1_2/2
E=x(2);
T=x(1);
I(1)=x(3);
I(2)=x(4);
p=param(1);
g=param(2);
k(1,1)=param(3);
k(1,2)=param(5);
k(1,3)=param(6)-param(5);
k(2,1)=param(3);
k(2,2)=param(5)/2;
k(2,3)=param(5)/2;
k(2,4)=param(6)-param(5);

dx(2,1)=p*E*T/(g+T)-k(1,1)*E*T-d*E+k(1,3)*I(1)-k(2,1)*I(1)*E+2*k(2,4)*I(2)+k(2,2)*I(2);
dx(1,1)=a*T*(1-b*T)-k(1,1)*E*T+k(1,2)*I_1+k(2,3)*I_2+k(2,2)*I(2);
dx(3,1)=k(1,1)*E*T-k(1,3)*I(1)-k(1,2)*I(1)-k(2,1)*I(1)*E;
dx(4,1)=k(2,1)*I(1)*E-k(2,4)*I(2)-k(2,3)*I(2)-k(2,2)*I(2);