
function dx = tumor_only(t, x, param) 

T = x(1); 
dx(1,1) =  param(1)*T*(1-param(2)*T);

end 