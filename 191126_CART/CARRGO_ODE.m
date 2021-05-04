function CARRGO_ODE() 

tspan = [0 2];
% 1dim 
y0 = 10;
%[t,y] = ode45(@(t,y) 0.2*y, tspan, y0);

% 2dim 
y0 = [150 20]';
[t,y] = ode45(@myfun, tspan, y0);

figure; hold on; 
plot(t,y,'-o')

end 


function dy = myfun( t, y ) 


%% define parameters 
rho = 1; 


% y(1), y(2) 
dy(1,1) = rho*(y(1))*(1-(y(1))/1000)-0.02*y(1)*(y(2)) ; 
dy(2,1) =  0.03*(y(1))*(y(2))-0.1*(y(2));  

end 