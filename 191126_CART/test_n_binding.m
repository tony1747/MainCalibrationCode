
[num,txt,raw]=xlsread('fitted params (PBT138, slow, k_n1 included).xlsx');
choices=txt(:,1);
choices(1,:)=[];
m=size(choices);
m=m(1);
TumVol=zeros(27,5,2);
for i=1:m
    for n_bind=1:5
    choice=choices(i);
    choice=char(choice);
    [FinalTumVol,TotalTumVol,CARTnum,receptor,n]=multiple_binding_solver(choice,n_bind);
    TumVol(i,n_bind,1)=FinalTumVol;
    TumVol(i,n_bind,2)=TotalTumVol;
    end
    figure;
    hold on
    bar(TumVol(i, 1:5, 1))
    plot( [1,2,3,4,5], TumVol(i, 1:5, 1) ,'x', 'markeredgecolor', 'k', 'markersize', 14 ) %Final Tumor volume
   
  %  xlim([0.7,5.3])
    %currentylim=ylim;
    ylabel('Final Tumor Volume')
    xlabel( 'binding number' )
   % xticks( [1:5] );
    yl=[min(TumVol(i, 1:5, 1)),max(TumVol(i, 1:5, 1))];
    ylim(yl)
    set(gca, 'FontSize',20);  box on;
    xticks([1,2,3,4,5])
    hold off
    TumVolPlot=strcat('Scenario_2/FinalTumVol',receptor,int2str(CARTnum),'_',int2str(n), '.jpg');
    saveas(gcf,TumVolPlot)
%     ylim([0,currentylim(2)*1.1])
   % TumVolPlot=strcat('./Scenario_2/FinalTumVol',receptor,int2str(CARTnum),'_',int2str(n), '_fromzero.jpg');
   % saveas(gcf,TumVolPlot)
    TumVolPlot=strcat('./Scenario_2/FinalTumVol',receptor,int2str(CARTnum),'_',int2str(n), '.fig');
    saveas(gcf,TumVolPlot)
   % plot( [1,2,3,4,5], TumVol(i, 1:5, 2) , 'x', 'markersize', 14  )
   % xlim([0.7,5.3])
    %currentylim=ylim;
%     ylim([0,currentylim(2)])
   % xticks( [1:5] );
   % set(gca, 'FontSize',20);  box on;  
   % xlabel( 'binding number' )
   % ylabel('Total Tumor Volume')
   % TumVolPlot=strcat('./Scenario_2/TotalTumVol',receptor,int2str(CARTnum),'_',int2str(n), '.jpg');
   % saveas(gcf,TumVolPlot)
end
