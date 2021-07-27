function dx = three_binding_slow(t, x, param)
global a b d %n d p g m

%param={p,g,k_1_1,k_1_2,k_1_3,k_n1}
%param={p,g,k_1,alpha(=k_n1+k_1_3),beta(=k_n1+k_1_2),gamma(=k_n1+k_1_3+k_1_2)}
%up to five-binding
E=x(2);
T=x(1);

I(1)=x(3);
I(2)=x(4);
I(3)=x(5);
p=param(1);
g=param(2);
k_n(1:3)=param(6);
k(1,1)=param(3);
k(1,2)=param(4);
k(1,3)=param(5);
k(2,1)=param(3);
k(2,2)=param(4)/2;
k(2,3)=param(4)/2;
k(2,4)=param(5);
k(3,1)=param(3);
k(3,2)=param(4)/3;
k(3,3)=param(4)/3;
k(3,4)=param(4)/3;
k(3,5)=param(5);
sum_1=0;
sum_2=0;
for j=1:3 
    temp_1=0;
    temp_2=0;
    for i=1:j    
    temp_1=temp_1+i*k(j,i+2); %sum( k(j, 3:(j+2)) )
    temp_2=temp_2+k(j,i+1);
    end
    temp_1=(temp_1+k_n(j))*I(j);
    sum_1=sum_1+temp_1;
    temp_2=temp_2*I(j);
    sum_2=sum_2+temp_2;
end
temp_3=0;
for j=1:2
    temp_3=temp_3+k(j+1,1)*I(j);
end
temp_3=temp_3*E;
    



dx(2,1)=p*E*T/(g+T)-k(1,1)*E*T-d*E+sum_1-temp_3;
dx(1,1)=a*T*(1-b*T)-k(1,1)*E*T+sum_2+k_n(1)*I(1);
dx(3,1)=k(1,1)*E*T-k(1,3)*I(1)-k(1,2)*I(1)-k(2,1)*I(1)*E+k_n(2)*I(2);
dx(4,1)=k(2,1)*I(1)*E-(sum(k(2,2:4))+k(3,1)*E+k_n(2))*I(2)+k_n(3)*I(3);
dx(5,1)=k(3,1)*I(2)*E-(sum(k(3,2:5))+k_n(3))*I(3);




