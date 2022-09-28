for receptor = ['L', 'M', 'H']
  for CARTnum = [0, 5, 10, 20]
    read_CARTdata() %% <- change  last part, data to cancerdata <- maybe you can also make something to select the cancer cell type outside. 
    if CARTnum ~= 0
        CARTratio=1/CARTnum;
    
    end
    
    for n=1:3
        if receptor=='M'
          selected_cancerdata(:,n) = selected_cancerdata(:,n)*3;
        end
        figure;
    
        plot(time-time(1),selected_cancerdata(:,n),'ok','MarkerSize',6,'MarkerFaceColor','k')
        %ax=gca;
        %ax.FontSize = 5;
        xlabel('Time','FontSize',14)
        ylabel('Tumor Size','FontSize',14)
        ylim([0,4])
        xlim([0,3.5])
        set(gca,'FontSize',32)
        cancerdata_plot=strcat('Run8_PBT138_data_plot',receptor,int2str(CARTnum),'_',int2str(n),'.jpg');
        saveas(gcf,cancerdata_plot)
    end
    
    %% here save parameters/posterior 
  end
end 


