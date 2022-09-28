function dx = tumor_cart_only(t, x, param)
%param={p,m,n,g}
global a b d  %n  p g m
E=x(2);
T=x(1);
dx(2,1)=param(1)*E*T/(param(4)+T)-param(2)*E*T-d*E;
dx(1,1)=a*T*(1-b*T)-param(3)*E*T;
end

% dx(1,1)=param(1)*E*T/(param(2)+T)-m*E*T-d*E;
% dx(2,1)=a*T*(1-b*T)-n*E*T;

% dx(1,1)=p *E*T/(param(3)+T)-param(1)*E*T-d*E;
% dx(2,1)=a*T*(1-b*T)-param(2)*E*T; 

% dx(1,1)=p*E*T/(param(3)+T)-param(1)*E*T-d*E;
% dx(2,1)=a*T*(1-b*T)-param(2)*E*T; 


