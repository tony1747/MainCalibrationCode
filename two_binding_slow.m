function dx = two_binding_slow(t, x, param)
global a b d %n d p g m
%param={p,g,k_1,k_2,k_3,k_n1}
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
k_n(1:2)=param(6);
k(1,1)=param(3);
k(1,2)=param(4);
k(1,3)=param(5);
k(2,1)=param(3);
k(2,2)=param(4)/2;
k(2,3)=param(4)/2;
k(2,4)=param(5);

dx(2,1)=p*E*T/(g+T)-k(1,1)*E*T-d*E+(k(1,3)+k_n(1))*I(1)-k(2,1)*I(1)*E+2*k(2,4)*I(2)+k(2,3)*I(2)+k_n(2)*I(2);
dx(1,1)=a*T*(1-b*T)-k(1,1)*E*T+k_n(1)*I(1)+k(1,2)*I(1)+k(2,2)*I(2)+k(2,3)*I(2);
dx(3,1)=k(1,1)*E*T-k(1,3)*I(1)-k(1,2)*I(1)-k(2,1)*I(1)*E-k_n(1)*I(1)+k_n(2)*I(2);
dx(4,1)=k(2,1)*I(1)*E-k(2,4)*I(2)-k(2,3)*I(2)-k(2,2)*I(2)-k_n(2)*I(2);