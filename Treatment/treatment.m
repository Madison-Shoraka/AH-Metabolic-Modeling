%% standard parameters
y3=1; %Time delay (min) between alcohol metabolism and acetyl-CoA production
tau=10000;
dt=.005;
T=25000;
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
alcohol_met=-.05;
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

mM_drink=17.4; %mM alcohol per drink
perturb=[200];

% failed decompensated
decay= 0.001;
y1=0.01;
y2=0.1;

%% no treatment
treatment_time=-1;


[concentration_normal,velocities_normal,alcohol,H2O2_normal,effect,total_oxidative_stress_normal]=treatment_cortocosteriods(treatment_time,tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
     
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_normal)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_normal)-0


%% treatment 1
time_IC=2880;
number_treatments=7;
treatment_time1=time_IC:1440:number_treatments*1440+time_IC;
[concentration_treated1,velocities_treated1,alcohol,H2O2_treated1,effect1,total_oxidative_stress_treatment1]=treatment_cortocosteriods(treatment_time1,tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
 
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_treated1)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_treated1)-0

%% treatment 2
time_IC=4320;
number_treatments=7;
treatment_time2=time_IC:1440:number_treatments*1440+time_IC;
[concentration_treated2,velocities_treated2,alcohol,H2O2_treated2,effect2,total_oxidative_stress_treatment2]=treatment_cortocosteriods(treatment_time2,tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
 
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_treated2)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_treated2)-0

%% treatment 3
time_IC=5760;
number_treatments=7;
treatment_time3=time_IC:1440:number_treatments*1440+time_IC;
[concentration_treated3,velocities_treated3,alcohol,H2O2_treated3,effect3,total_oxidative_stress_treatment3]=treatment_cortocosteriods(treatment_time3,tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
 
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_treated3)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_treated3)-0

%% plot results

figure
t=tiledlayout(2,2);
title(t,"Concentration loads")
ax1=nexttile;
hold on
plot(ax1, time, concentration_normal(:,1)/concentration_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax1,time, concentration_normal(:,2)/concentration_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax1,time, concentration_normal(:,3),'LineWidth',3,'Color',[0,0,1])
plot(ax1,time, concentration_normal(:,4),'LineWidth',3,'Color',[1,1,0])
plot(ax1,time, concentration_normal(:,5),'LineWidth',3,'Color',[0,1,1])
grid on;
hold off
xlabel(ax1,'time (min)');
ylabel(ax1,'concentration load');
title(ax1,'No treatment');
ax2=nexttile;
hold on
plot(ax2, time, concentration_treated1(:,1)/concentration_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax2,time, concentration_treated1(:,2)/concentration_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax2,time, concentration_treated1(:,3),'LineWidth',3,'Color',[0,0,1])
plot(ax2,time, concentration_treated1(:,4),'LineWidth',3,'Color',[1,1,0])
plot(ax2,time, concentration_treated1(:,5),'LineWidth',3,'Color',[0,1,1])
grid on;
xlabel(ax2,'time (min)');
ylabel(ax2,'concentration load');
title(ax2,'Treatment 1');
ax3=nexttile;
hold on
plot(ax3, time, concentration_treated2(:,1)/concentration_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax3,time, concentration_treated2(:,2)/concentration_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax3,time, concentration_treated2(:,3),'LineWidth',3,'Color',[0,0,1])
plot(ax3,time, concentration_treated2(:,4),'LineWidth',3,'Color',[1,1,0])
plot(ax3,time, concentration_treated2(:,5),'LineWidth',3,'Color',[0,1,1])
grid on;
xlabel(ax3,'time (min)');
ylabel(ax3,'concentration load');
title(ax3,'Treatment 2');
ax4=nexttile;
hold on
plot(ax4, time, concentration_treated3(:,1)/concentration_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax4,time, concentration_treated3(:,2)/concentration_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax4,time, concentration_treated3(:,3),'LineWidth',3,'Color',[0,0,1])
plot(ax4,time, concentration_treated3(:,4),'LineWidth',3,'Color',[1,1,0])
plot(ax4,time, concentration_treated3(:,5),'LineWidth',3,'Color',[0,1,1])
grid on;
xlabel(ax4,'time (min)');
ylabel(ax4,'concentration load');
title(ax4,'Treatment 3');
fontsize(18,"points") 
legend(["[GSH]_m","[MET]_c","[ammonia]_m","[a-CoA]_m","[malonyl-CoA]_m"])


figure
t=tiledlayout(2,2);
title(t,"Scaled velocities")
ax1=nexttile;
hold on
plot(ax1,time, velocities_normal(:,1)./velocities_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax1,time, (velocities_normal(:,2))./velocities_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax1,time, (velocities_normal(:,3))./velocities_IC(3),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
plot(ax1,time, (velocities_normal(:,4))./velocities_IC(4),'LineWidth',3,'Color',[1,1,0])
ylim([.8,1.2])
grid on;
hold off
xlabel(ax1,'time (min)');
ylabel(ax1,'scaled velocities');
title(ax1,'No treatment');
ax2=nexttile;
hold on
ylim([.8,1.2])
plot(ax2,time, velocities_treated1(:,1)./velocities_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax2,time, (velocities_treated1(:,2))./velocities_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax2,time, (velocities_treated1(:,3))./velocities_IC(3),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
plot(ax2,time, (velocities_treated1(:,4))./velocities_IC(4),'LineWidth',3,'Color',[1,1,0])
grid on;
xlabel(ax2,'time (min)');
ylabel(ax2,'scaled velocities');
title(ax2,'treatment 1');
ax3=nexttile;
hold on
plot(ax3,time, velocities_treated2(:,1)./velocities_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax3,time, (velocities_treated2(:,2))./velocities_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax3,time, (velocities_treated2(:,3))./velocities_IC(3),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
plot(ax3,time, (velocities_treated2(:,4))./velocities_IC(4),'LineWidth',3,'Color',[1,1,0])
ylim([.8,1.2])
grid on;
xlabel(ax3,'time (min)');
ylabel(ax3,'scaled velocities');
title(ax3,'treatment 2');
ax4=nexttile;
hold on
plot(ax4,time, velocities_treated3(:,1)./velocities_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax4,time, (velocities_treated3(:,2))./velocities_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax4,time, (velocities_treated3(:,3))./velocities_IC(3),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
plot(ax4,time, (velocities_treated3(:,4))./velocities_IC(4),'LineWidth',3,'Color',[1,1,0])
ylim([.8,1.2])
grid on;
xlabel(ax4,'time (min)');
ylabel(ax4,'scaled velocities');
title(ax4,'treatment 3');
fontsize(18,"points") 
legend(["GSH","MET","Urea","FABs"])

figure
t=tiledlayout(2,2);
title(t,"[H_2O_2]_m")
ax1=nexttile;
hold on
plot(ax1,time,total_oxidative_stress_normal(:),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax1,'time (min)');
ylabel(ax1,'concentration \mu M');
title(ax1,'no treatment');
ax2=nexttile;
hold on
plot(ax2,time,total_oxidative_stress_treatment1(:),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax2,'time (min)');
ylabel(ax2,'concentration \mu M');
title(ax2,'treatment1');
ax3=nexttile;
hold on
plot(ax3,time,total_oxidative_stress_treatment2(:),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax3,'time (min)');
ylabel(ax3,'concentration \mu M');
title(ax3,'treatment2');
ax4=nexttile;
hold on
plot(ax4,time,total_oxidative_stress_treatment3(:),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax4,'time (min)');
ylabel(ax4,'concentration \mu M');
title(ax4,'treatment3');
fontsize(18,"points") 


figure
t=tiledlayout(2,2);
title(t,"[H_2O_2] directly generated from ethanol metabolism")
ax1=nexttile;
hold on
plot(ax1,time,H2O2_normal(:)*1000,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax1,'time (min)');
ylabel(ax1,'[H_2O_2] mM');
title(ax1,'no treatment');
ax2=nexttile;
hold on
plot(ax2,time,H2O2_treated1(:)*1000,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax2,'time (min)');
ylabel(ax2,'[H_2O_2] mM');
title(ax2,'treatment1');
ax3=nexttile;
hold on
plot(ax3,time,H2O2_treated2(:)*1000,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax3,'time (min)');
ylabel(ax3,'[H_2O_2] mM');
title(ax3,'treatment2');
ax4=nexttile;
hold on
plot(ax4,time,H2O2_treated3(:)*1000,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax4,'time (min)');
ylabel(ax4,'[H_2O_2] mM');
title(ax4,'treatment3');
fontsize(18,"points") 


grid on;
hold off
figure
hold on
plot(time, effect1,'LineWidth',3,'Color',[1, 0, 0])
plot(time, effect2,'LineWidth',3,'Color',[1, 1, 0])
plot(time, effect3,'LineWidth',3,'Color',[0, 1, 0])
xlabel('Time (min)');
ylabel('k_c effect on GSH cycle');
title('k_c effect on GSH cycle');
grid on;
hold off





