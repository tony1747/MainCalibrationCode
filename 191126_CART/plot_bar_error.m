[num,txt,raw]=xlsread('fitted params (PBT138, slow, k_n1 included).xlsx');
slowerror=num(25:27,7);
fasterror=num(25:27,8);
figure; 

bar( mean( [fasterror, slowerror] ) )
hold on; 

plot( [1 1 1], fasterror, 'kx', 'markersize', 10, 'linewidth', 2 )
plot( [2 2 2], slowerror, 'kx', 'markersize', 10, 'linewidth', 2 )


set(gca,'FontSize',14)
ylabel( 'error' )
xticks([1 2])
xticklabels({'fast', 'slow'})
xlim( [0 3] )
plot=strcat('bar_plot', '.jpg');
saveas(gcf,plot)