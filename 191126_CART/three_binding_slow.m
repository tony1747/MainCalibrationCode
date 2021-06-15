function dx = three_binding_slow(t, x, param)
global a b d %n d p g m

%param={p,g,k_1_1,k_1_2,k_1_3,k_2_1,k_2_2,k_2_3,k_2_4,
%k_3_1,k_3_2,k_3_3,k_3_4,k_3_5}
%param={p,g,k_1,alpha(=k_n1+k_1_3),beta(=k_n1+k_1_2),gamma(=k_n1+k_1_3+k_1_2)}
E=x(2);
T=x(1);

I(1)=x(3);
I(2)=x(4);
I(3)=x(5);
p=param(1);
g=param(2);
k(1,1)=param(3);
k(1,2)=param(5);
k(1,3)=param(6)-param(5);
k(2,1)=param(3);
k(2,2)=param(5)/2;
k(2,3)=param(5)/2;
k(2,4)=param(6)-param(5);
k(3,1)=param(3);
k(3,2)=param(5)/3;
k(3,3)=param(5)/3;
k(3,4)=param(5)/3;
k(3,5)=param(6)-param(5);
temp_1=zeros(1,3);
temp_3=zeros(1,3);
for j=1:3 
    temp_var_1=0;
    temp_var_2=0;
    for i=1:j    
    temp_var_1=temp_var_1+i*k(j,i+2); %sum( k(j, 3:(j+2)) )
    temp_var_2=temp_var_2+k(j,i+1);
    end
    temp_1(j)=temp_var_1*I(j);
    temp_3(j)=temp_var_2*I(j);
end
temp_1=sum(temp_1);
temp_2=sum(k(2:3,1).*I(1:2))*E;
temp_3=sum(temp_3);
temp_4=0;
for i=2:4
    temp_4=temp_4+k(2,i)+k(3,1)*E;
end
temp_4=I(2)*temp_4;
temp_5=sum(k(3,2:5))*I(3);

dx(2,1)=p*E*T/(g+T)-k(1,1)*E*T-d*E+temp_1-temp_2;
dx(1,1)=a*T*(1-b*T)-k(1,1)*E*T+temp_3;
dx(3,1)=k(1,1)*E*T-k(1,3)*I(1)-k(1,2)*I(1)-k(2,1)*I(1)*E;
dx(4,1)=k(2,1)*I(1)*E-temp_4;
dx(5,1)=k(3,1)*I(2)*E-temp_5;




