function [FinalTumVol,TotalTumVol,CARTnum,receptor,n]=multiple_binding_solver(choice,n_bind)%for example, choice='H, 1:10, 1'
global a b d
%b=1/6;
d=0.0412;

[num,txt,raw]=xlsread('fitted params (PBT138, slow, k_n1 included).xlsx');

row=find(strcmp(txt,choice));
params=num(row-1,:);    %matrix of txt has 28 rows, while matrix of num has 27 rows because
%it only counts numeric values
receptor=choice(1);
[a,b] = set_CancerGrowthParams( receptor );
ratio=extractBetween(choice,' ',',','Boundaries','exclusive');%extract 1:10 out of the choice
if string(ratio{1})=='1:5'
    CARTnum=5;
elseif string(ratio{1})=='1:10'
    CARTnum=10;
elseif string(ratio{1})=='1:20'
    CARTnum=20;
end
CARTratio=1/CARTnum;
read_CARTdata()
n=str2double(choice(end));  %read_CARTdata provides 3 sets of data with same density and ratio,
%choose the set that we specified in the string
%'H, 1:10, 1'
data=selected_cancerdata(:,n);
switch( n_bind )
    case 1
        x(1)=data(1);
        x(2)=data(1)*CARTratio;
        x(3:3)=0; %3:5 for three-binding; 3:6 for four-binding ...
        [t,tumVol] = ode23(@(t,y)one_binding_slow(t,y,params), time, x, params );%change the ODE for different bindings
        t=t-t(1);
        diff=zeros(1,length(t));
        for i=1:(length(t)-1)
            diff(i)=t(i+1)-t(i);
        end
        TotalTumVol=0;
        for i=1:(length(t)-1)
            TotalTumVol=TotalTumVol+diff(i)*(tumVol(i,1)+tumVol(i+1,1))/2;
        end
        FinalTumVol=tumVol(end,1);
        figure;hold on;
        plot(t,tumVol(:,1),'-r')
        plot(t,tumVol(:,2),'-b')
        for i=3:3 %same as above, 3:5 for three-binding, 3:6 for four-binding...
            plot(t,tumVol(:,i),'y')
        end
        xlabel('Time','FontSize',14)
        ylabel('Tumor Size','FontSize',14)
        set(gca,'FontSize',14)
        legend('Cancer','CAR-T','Conjugate1','Conjugate2','Conjugate3')
        n_binding=strcat('1_binding_slow(with conjugates)',receptor,int2str(CARTnum),'_',int2str(n),'.jpg');
        %saveas(gcf,n_binding)
        hold off;
        plot(t,tumVol(:,1),'-r', t,tumVol(:,2),'-b')
        xlabel('Time','FontSize',14)
        ylabel('Tumor Size','FontSize',14)
        set(gca,'FontSize',14)
        legend('Cancer','CAR-T')
        n_binding=strcat('1_binding_slow(without conjugates)',receptor,int2str(CARTnum),'_',int2str(n),'.jpg');
        %saveas(gcf,n_binding)
        close all
    case 2
        x(1)=data(1);
        x(2)=data(1)*CARTratio;
        x(3:4)=0; %3:5 for three-binding; 3:6 for four-binding ...
        [t,tumVol] = ode23(@(t,y)two_binding_slow(t,y,params), time, x, params );%change the ODE for different bindings
        t=t-t(1);
        diff=zeros(1,length(t));
        for i=1:(length(t)-1)
            diff(i)=t(i+1)-t(i);
        end
        TotalTumVol=0;
        for i=1:(length(t)-1)
            TotalTumVol=TotalTumVol+diff(i)*(tumVol(i,1)+tumVol(i+1,1))/2;
        end
        FinalTumVol=tumVol(end,1);
        
        
        figure;hold on;
        plot(t,tumVol(:,1),'-r')
        plot(t,tumVol(:,2),'-b')
        for i=3:4 %same as above, 3:5 for three-binding, 3:6 for four-binding...
            plot(t,tumVol(:,i),'y')
        end
        xlabel('Time','FontSize',14)
        ylabel('Tumor Size','FontSize',14)
        set(gca,'FontSize',14)
        legend('Cancer','CAR-T','Conjugate1','Conjugate2','Conjugate3')
        n_binding=strcat('2_binding_slow(with conjugates)',receptor,int2str(CARTnum),'_',int2str(n),'.jpg');
        %saveas(gcf,n_binding)
        hold off;
        plot(t,tumVol(:,1),'-r', t,tumVol(:,2),'-b')
        xlabel('Time','FontSize',14)
        ylabel('Tumor Size','FontSize',14)
        set(gca,'FontSize',14)
        legend('Cancer','CAR-T')
        n_binding=strcat('2_binding_slow(without conjugates)',receptor,int2str(CARTnum),'_',int2str(n),'.jpg');
        %saveas(gcf,n_binding)
        close all
    case 3
        x(1)=data(1);
        x(2)=data(1)*CARTratio;
        x(3:5)=0; %3:5 for three-binding; 3:6 for four-binding ...
        [t,tumVol] = ode23(@(t,y)three_binding_slow(t,y,params), time, x, params );%change the ODE for different bindings
        t=t-t(1);
        diff=zeros(1,length(t));
        for i=1:(length(t)-1)
            diff(i)=t(i+1)-t(i);
        end
        TotalTumVol=0;
        for i=1:(length(t)-1)
            TotalTumVol=TotalTumVol+diff(i)*(tumVol(i,1)+tumVol(i+1,1))/2;
        end
        FinalTumVol=tumVol(end,1);
       
        figure;hold on;
        plot(t,tumVol(:,1),'-r')
        plot(t,tumVol(:,2),'-b')
        for i=3:5 %same as above, 3:5 for three-binding, 3:6 for four-binding...
            plot(t,tumVol(:,i),'y')
        end
        xlabel('Time','FontSize',14)
        ylabel('Tumor Size','FontSize',14)
        set(gca,'FontSize',14)
        legend('Cancer','CAR-T','Conjugate1','Conjugate2','Conjugate3')
        n_binding=strcat('3_binding_slow(with conjugates)',receptor,int2str(CARTnum),'_',int2str(n),'.jpg');
        %saveas(gcf,n_binding)
        hold off;
        plot(t,tumVol(:,1),'-r', t,tumVol(:,2),'-b')
        xlabel('Time','FontSize',14)
        ylabel('Tumor Size','FontSize',14)
        set(gca,'FontSize',14)
        legend('Cancer','CAR-T')
        n_binding=strcat('3_binding_slow(without conjugates)',receptor,int2str(CARTnum),'_',int2str(n),'.jpg');
        %saveas(gcf,n_binding)
        close all
    case 4
        x(1)=data(1);
        x(2)=data(1)*CARTratio;
        x(3:6)=0; %3:5 for three-binding; 3:6 for four-binding ...
        [t,tumVol] = ode23(@(t,y)four_binding_slow(t,y,params), time, x, params );%change the ODE for different bindings
        t=t-t(1);
        diff=zeros(1,length(t));
        for i=1:(length(t)-1)
            diff(i)=t(i+1)-t(i);
        end
        TotalTumVol=0;
        for i=1:(length(t)-1)
            TotalTumVol=TotalTumVol+diff(i)*(tumVol(i,1)+tumVol(i+1,1))/2;
        end
        FinalTumVol=tumVol(end,1);
       
        figure;hold on;
        plot(t,tumVol(:,1),'-r')
        plot(t,tumVol(:,2),'-b')
        for i=3:6 %same as above, 3:5 for three-binding, 3:6 for four-binding...
            plot(t,tumVol(:,i),'y')
        end
        xlabel('Time','FontSize',14)
        ylabel('Tumor Size','FontSize',14)
        set(gca,'FontSize',14)
        legend('Cancer','CAR-T','Conjugate1','Conjugate2','Conjugate3','Conjugate4')
        n_binding=strcat('4_binding_slow(with conjugates)',receptor,int2str(CARTnum),'_',int2str(n),'.jpg');
        %saveas(gcf,n_binding)
        hold off;
        plot(t,tumVol(:,1),'-r', t,tumVol(:,2),'-b')
        xlabel('Time','FontSize',14)
        ylabel('Tumor Size','FontSize',14)
        set(gca,'FontSize',14)
        legend('Cancer','CAR-T')
        n_binding=strcat('4_binding_slow(without conjugates)',receptor,int2str(CARTnum),'_',int2str(n),'.jpg');
        %saveas(gcf,n_binding)
        close all
    case 5
        x(1)=data(1);
        x(2)=data(1)*CARTratio;
        x(3:7)=0; %3:5 for three-binding; 3:6 for four-binding ...
        [t,tumVol] = ode23(@(t,y)five_binding_slow(t,y,params), time, x, params );%change the ODE for different bindings
        t=t-t(1);
        diff=zeros(1,length(t));
        for i=1:(length(t)-1)
            diff(i)=t(i+1)-t(i);
        end
        TotalTumVol=0;
        for i=1:(length(t)-1)
            TotalTumVol=TotalTumVol+diff(i)*(tumVol(i,1)+tumVol(i+1,1))/2;
        end
        FinalTumVol=tumVol(end,1);
       
        figure;hold on;
        plot(t,tumVol(:,1),'-r')
        plot(t,tumVol(:,2),'-b')
        for i=3:7 %same as above, 3:5 for three-binding, 3:6 for four-binding...
            plot(t,tumVol(:,i),'y')
        end
        xlabel('Time','FontSize',14)
        ylabel('Tumor Size','FontSize',14)
        set(gca,'FontSize',14)
        legend('Cancer','CAR-T','Conjugate1','Conjugate2','Conjugate3','Conjugate4','Conjugate5')
        n_binding=strcat('5_binding_slow(with conjugates)',receptor,int2str(CARTnum),'_',int2str(n),'.jpg');
        %saveas(gcf,n_binding)
        hold off;
        plot(t,tumVol(:,1),'-r', t,tumVol(:,2),'-b')
        xlabel('Time','FontSize',14)
        ylabel('Tumor Size','FontSize',14)
        set(gca,'FontSize',14)
        legend('Cancer','CAR-T')
        n_binding=strcat('5_binding_slow(without conjugates)',receptor,int2str(CARTnum),'_',int2str(n),'.jpg');
        %saveas(gcf,n_binding)
        close all
end
end



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