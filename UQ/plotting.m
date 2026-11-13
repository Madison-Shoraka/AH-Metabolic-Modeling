
sI=[ 0.4138    0.0645    0.0935    0.4397
    0.2846    0.6028    0.5426    0.0948
    0.1040    0.2088    0.0871    0.3450];
tI=[0.6211    0.1537    0.2483    0.5332
    0.4477    0.7145    0.7363    0.1675
    0.2090    0.2480    0.1794    0.4257];

GSH=[sI(:,1),tI(:,1)];
MET=[sI(:,2),tI(:,2)];
urea=[sI(:,3),tI(:,3)];
FABs=[sI(:,4),tI(:,4)];
cat=categorical({'k_d','y_1','y_2'});



t=tiledlayout(2,2);
nexttile
hold on
title("GSH")
bar(cat,GSH,1)
fig = gcf;
fig.Color = 'w' ;                  % White
fontsize(fig, 20, "points")
ylim([0,1])
hold off
nexttile
hold on
title("MET")
bar(cat,MET,1)
fig = gcf;
fig.Color = 'w' ;                  % White
fontsize(fig, 20, "points")
ylim([0,1])
hold off
nexttile
hold on
bar(cat,urea,1)
title("Urea")
fig = gcf;
fig.Color = 'w' ;                  % White
fontsize(fig, 20, "points")
ylim([0,1])
hold off
nexttile
hold on
bar(cat,FABs,1)
title("FABs")
fig = gcf;
fig.Color = 'w' ;                  % White
fontsize(fig, 20, "points")
legend(["Sensitivity Index","Total-effect Index"])
ylim([0,1])
hold off

saveas(t,"UQ_fig_onehit.png",'png')



sI=[0.4107    0.0896    0.1854    0.4076
    0.1906    0.5254    0.2924    0.1008
    0.1001    0.2094    0.1687    0.4005];
tI=[0.6741    0.2624    0.5993    0.5039
    0.4122    0.6934    0.6777    0.1987
    0.2646    0.2695    0.5680    0.4797];

GSH=[sI(:,1),tI(:,1)];
MET=[sI(:,2),tI(:,2)];
urea=[sI(:,3),tI(:,3)];
FABs=[sI(:,4),tI(:,4)];
cat=categorical({'k_d','y_1','y_2'});



t=tiledlayout(2,2);
nexttile
hold on
title("GSH")
bar(cat,GSH,1)
fig = gcf;
fig.Color = 'w' ;                  % White
fontsize(fig, 20, "points")
ylim([0,1])
hold off
nexttile
hold on
title("MET")
bar(cat,MET,1)
fig = gcf;
fig.Color = 'w' ;                  % White
fontsize(fig, 20, "points")
ylim([0,1])
hold off
nexttile
hold on
bar(cat,urea,1)
title("Urea")
fig = gcf;
fig.Color = 'w' ;                  % White
fontsize(fig, 20, "points")
ylim([0,1])
hold off
nexttile
hold on
bar(cat,FABs,1)
title("FABs")
fig = gcf;
fig.Color = 'w' ;                  % White
fontsize(fig, 20, "points")
legend(["Sensitivity Index","Total-effect Index"])
ylim([0,1])
hold off

saveas(t,"UQ_fig_twohit.png",'png')
