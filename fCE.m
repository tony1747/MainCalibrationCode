
function value = fCE( t, x, param) 
global a b %n d p g m


C = x(1); 
E = x(2); 
value(1,1) =  param(1)*C*(1-b*C) - param(3)*C*E; 
% value(2,1) = -d*E + param(1)*E*(C/(param(4)+C))-param(2)*E*C;
value(2,1) = -param(4)*E + param(1)*E*(C/(param(5)+C))-param(2)*E*C;
            

end 