%% standard parameters
y3=1; %Time delay (min) between alcohol metabolism and acetyl-CoA production
tau=10000;
dt=.005;
T=1640;
time=0:dt:T-dt;
nt=length(time);
sigma=5;
time_to_detox=60;
response_rate=.01;
l=[1 0 0 0 %metabolite = GSH_m
   0 1 0 0 %metbolite = SAM_c
   0 0 1 0 %metabolite = ammonia_m
   0 0 1 1 %metabolite = a-CoA_m
   0 0 0 1]; %metabolite=malonyl-CoA_m
type=[1 1 0 0 0]; %1=produced, 0=used
beta = 7; % sharpness parameter
alcohol_met=-.05*5.5;
Adj_mat=[0     7     4     2
     5     0     2     1
     6     5     0     4
     2     1     3     0];
velocities_IC= [2.0737    0.4002    0.6193    0.9068]; %GSH,MET,Urea,FABs


concentration_IC=[10 .08 0 0 0];
critical_concentration_vec=[7.5,.04,1,2,.3];
critical_concentrations=[1 0 0 0 %metabolite = GSH_m
    0 1 0 0 %metbolite = SAM_c
    0 0 1 0 %metabolite = ammonia_m
    0 0 1/2 1/2 %metabolite = a-CoA_m
    0 0 0 1]; %metabolite = malonyl-CoA_m 
critical_concentrations=critical_concentration_vec'.*critical_concentrations;
max_concentration=[concentration_IC(1) concentration_IC(2) 100 100 100];

mM_drink=100; %mM alcohol per drink
perturb=[200]; %minutes from start, can include multiple (24 hours later)

% failed decompensated
decay= 0.01;
y1=0.06;
y2=0.15;


[concentration,velocities,alcohol,ROS,total_oxidative_stress]=full_model(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
     
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(ROS)-0
disp("At 2 hours GSH concentration (mM)")
disp(concentration(64000,1))
disp("At 6 hours GSH concentration (mM)")
disp(concentration(112000,1))
disp("At 24 hours GSH concentration (mM)")
disp(concentration(end,1))


figure()
t=tiledlayout(2,2);
title(t,"Mouse Study")
ax1=nexttile;
hold on
plot(ax1,time, velocities(:,1)./velocities_IC(1),'LineWidth',3)
plot(ax1,time, (velocities(:,2))./velocities_IC(2),'LineWidth',3)
plot(ax1,time, (velocities(:,3))./velocities_IC(3),'LineWidth',3)
plot(ax1,time, (velocities(:,4))./velocities_IC(4),'LineWidth',3)
grid on;
hold off
xlabel(ax1,'time (min)');
legend(["GSH","MET","Urea","FABs"])
ylabel(ax1,[{'scaled velocities'},{'(mM/min)'}]);
title(ax1,'Scaled velocities (mM/min)');
ax2=nexttile;
hold on
plot(ax2,time, concentration(:,1)./concentration_IC(1),'LineWidth',3,'Color',[1, 0, 0])
ylim([0,1.1])
grid on;
hold off
xlabel(ax2,'time (min)');
ylabel(ax2,'Percent depletion');
title(ax2,'[GSH]_m depletion');
ax3=nexttile;
hold on
plot(ax3,time, alcohol,'LineWidth',3,'Color',[0.8500, 0.3250, 0.0980])
xlabel(ax3,'time (min)');
ylabel(ax3,'alcohol (mM)');
title(ax3,'Alcohol (mM) over time');
grid on;
hold off
ax4=nexttile;
hold on
plot(ax4,time, total_oxidative_stress,'LineWidth',3,'Color',[0.8500, 0.3250, 0.0980])
grid on;
hold off
xlabel(ax4,'time (min)');
ylabel(ax4,[{'Total oxidative stress'},{'(mM)'}]);
title(ax4,'Total oxidative stress over time');



figure()
hold on
plot(time, concentration(:,1),'LineWidth',3)
plot(time, (concentration(:,2)),'LineWidth',3)
plot(time, (concentration(:,3)),'LineWidth',3)
plot(time, (concentration(:,4)),'LineWidth',3)
plot(time, (concentration(:,5)),'LineWidth',3)
legend(["[GSH]_m","[SAMe]_c","[ammonia]_m","[a-CoA]_m","[malonyl-CoA]_m"])


