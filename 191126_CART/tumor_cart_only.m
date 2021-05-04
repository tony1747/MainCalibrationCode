function dx = tumor_cart_only(t, x, param)
%param={a,b,n,s,p,g,m,d}
global a b p m n d g  
E=x(1);
T=x(2);
dx(1,1)=param(1)*E*T/(param(5)+T)-param(2)*E*T-param(4)*E;
% dx(1,1)=param(1)*E*T/(param(4)+T)-param(2)*E*T-d*E;
dx(2,1)=a*T*(1-b*T)-param(3)*E*T;


% dx(1,1)=param(1)*E*T/(param(2)+T)-m*E*T-d*E;
% dx(2,1)=a*T*(1-b*T)-n*E*T;

% dx(1,1)=p *E*T/(param(3)+T)-param(1)*E*T-d*E;
% dx(2,1)=a*T*(1-b*T)-param(2)*E*T; 

% dx(1,1)=p*E*T/(param(3)+T)-param(1)*E*T-d*E;
% dx(2,1)=a*T*(1-b*T)-param(2)*E*T; 


