%draw the phase diagram for each fixed point of the Kuz model
%use the parameters from low density case with ratio 0.2
%(sqrt(d)+sqrt(mg))^2=1.233, (m/b+d)(gb+1)=1.358
global a b
a=0.3953;
b=1/6;
fp_y=[0,0,80.75,1.5842];
fp_x=[0,6,-2.65,0.1565];
p=2.8271;
m=0.0267;
n=1.8584 ;
d=0.4018;
g=8.5005;
time=[0:200];
params=[p,m,n,d,g];
CARTratio=0.2;
for i=1:4
    p_x=(fp_x(i)-1):0.1:(fp_x(i)+1);
    p_y=(fp_y(i)-1):0.1:(fp_y(i)+1);
    [x,y]=meshgrid(p_x,p_y);
    dx=p*x.*y./(g+y)-d*x-m*x.*y;
    dy=a*y.*(1-b*y)-n*x.*y;
   % u=gradient(dx);
   % v=gradient(dy);
    figure()
    quiver(x,y,dx,dy)  
    hold on
    [t,tumVol] = ode23(@(t,y)tumor_cart_only(t,y,params), time, [fp_x(i)+0.1,fp_y(i)-0.1]' );
    plot(tumVol(:,1),tumVol(:,2))
end 


%fitData.ydata=[0];



