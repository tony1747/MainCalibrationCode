
function dx = tumor_only(t, x, param) 
global b 
T = x(1); 
% dx(1,1) =  param(1)*T*(1-param(2)*T);
dx(1,1) =  param(1)*T*(1-b*T);

end 