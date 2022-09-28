%compute the equilibrium points 
%try fsolve
global a b d
[num,txt,raw]=xlsread('fitted params (PBT138, slow, k_n1 included).xlsx');
d=0.0412;

%------compute the equilibrium points by using the formula we derived in
%the paper------
%C=zeros(27,2);
%T=zeros(27,2);
%I=zeros(27,2);
%eigen=zeros(3,2,27);
%for i=1:27
%    cancer=txt(i+1,1);
%    CancerType=strcat(cancer{1}(1));
%    [a,b]=set_CancerGrowthParams( CancerType );
%    alpha=num(i,6)+num(i,5);
%    beta=num(i,6)+num(i,4);
%    gamma=num(i,6)+num(i,4)+num(i,5);
%    m=num(i,3)-alpha*num(i,3)/gamma;
%    n=num(i,3)-beta*num(i,3)/gamma;
%    C(i,1)=(num(i,1)-d-m*num(i,2)+sqrt((num(i,1)-d-m*num(i,2))^2-4*m*d*num(i,2)))/(2*m);
    
%    if isreal(C(i,1))==1 && C(i,1)>0
 %       C(i,2)=(num(i,1)-d-m*num(i,2)-sqrt((num(i,1)-d-m*num(i,2))^2-4*m*d*num(i,2)))/(2*m);
  %      T(i,1)=a*(1-b*C(i,1))/n;
  %      T(i,2)=a*(1-b*C(i,2))/n;
  %      J_1=Jacobian(a,b,num(i,1),num(i,2),d,num(i,3),alpha,beta,gamma,C(i,1),T(i,1));
  %      J_2=Jacobian(a,b,num(i,1),num(i,2),d,num(i,3),alpha,beta,gamma,C(i,2),T(i,2));
  %      eigen(:,1,i)=eig(J_1);
  %      eigen(:,2,i)=eig(J_2);
  %      I(i,1)=num(i,3)/gamma*C(i,1)*T(i,1);    %Volume for I may not be needed
  %      I(i,2)=num(i,3)/gamma*C(i,2)*T(i,2); 
  %      if num(i,1)<(k*alpha/(b*gamma)+d+k/b)*(g*b+1) && num(i,1)>(sqrt(d)+sqrt(m*g))^2
  %          disp('Condition 2 is satisfied')
  %      elseif num(i,1)>(k*alpha/(b*gamma)+d+k/b)*(g*b+1)
  %          disp('Condition 3 is satisfied')
  %      end
           
  %  else
  %      C(i,1)=0;       %for complex C, set it to be 0
  %      num(i,:)=0;     %since the fixed point is complex, we don't need to consider
                        %the corresponding parameter set. Let this set of
                        %parameters be 0
   % end
    
%end

%------compute the equilibrium points by using the built-in function
%jacobian(f,v)------
C=zeros(27,2);
T=zeros(27,2);
I=zeros(27,2);
eigen=zeros(3,2,27);

for i=1:27
    cancer=txt(i+1,1);
    CancerType=strcat(cancer{1}(1));
    [a,b]=set_CancerGrowthParams( CancerType );
    alpha=num(i,6)+num(i,5);
    beta=num(i,6)+num(i,4);
    gamma=num(i,6)+num(i,4)+num(i,5);
    m=num(i,3)-alpha*num(i,3)/gamma;
    num(i,1)=num(i,1)+(-num(i,3)*alpha/(b*gamma)+d+num(i,3)/b)*(num(i,2)*b+1);
    syms x y z %symbolic variables in order T, C, I
    matrix(1,1)=num(i,1)*x*y/(num(i,2)+y)-d*x-num(i,3)*x*y+alpha*z;
    matrix(2,1)=a*y*(1-b*y)-num(i,3)*x*y+beta*z;
    matrix(3,1)=num(i,3)*x*y-gamma*z;
   % matrix(x,y,z)=[eqn1 eqn2 eqn3];
    solns=solve([matrix(1)==0 matrix(2)==0 matrix(3)==0],[x y z]);
    
    solns.z=solns.z(solns.x~=0);
    solns.y=solns.y(solns.x~=0);
    solns.x=solns.x(solns.x~=0);
    
    C(i,1)=solns.y(1);
    C(i,2)=solns.y(2);
    
    T(i,1)=solns.x(1);
    %
    T(i,2)=solns.x(2);
    
    
    
    if isreal(C(i,:))==1 & C(i,:)>=0 & T(i,:)>=0
        I(i,1)=solns.z(1);
        I(i,2)=solns.z(2);
        J=jacobian(matrix,[x y z]);
        J_1=subs(J, [x y z], [T(i,1), C(i,1), I(i,1)]);
        J_2=subs(J, [x y z], [T(i,2), C(i,2), I(i,2)]);
        eigen(:,1,i)=eig(J_1);
        eigen(:,2,i)=eig(J_2);
        if num(i,1)<(-num(i,3)*alpha/(b*gamma)+d+num(i,3)/b)*(num(i,2)*b+1) && num(i,1)>(sqrt(d)+sqrt(m*num(i,2)))^2
            disp('Condition 2 is satisfied')
        elseif num(i,1)>(-num(i,3)*alpha/(b*gamma)+d+num(i,3)/b)*(num(i,2)*b+1)
            disp('Condition 3 is satisfied')
        else
            disp('neither condition 2 nor 3 is satisfied')
        end
        
    else
        C(i,:)=0;
        T(i,:)=0;
        num(i,:)=0;
    end
    
end


%get rid of complex fixed points and corresponding parameter sets


%C(~any(C,2),:)=[];
%T(~any(T,2),:)=[];
%I(~any(I,2),:)=[];
%num(~any(num,2),:)=[];
%eigen(:,:,~any(eigen,[1,2]))=[];

%save('eigenvalues','eigen')


function [a,b] = set_CancerGrowthParams( CancerType )

switch( CancerType )
    case 'L'
        a = 0.3953;
    case 'M'
        a = 0.3606;
    case 'H'
        a = 0.2187;
end
b = 1/6;
end
