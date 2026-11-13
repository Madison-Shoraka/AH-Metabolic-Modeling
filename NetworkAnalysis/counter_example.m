%% standard parameters
y3=1; %Time delay (min) between alcohol metabolism and acetyl-CoA production
tau=10000;
dt=.005;
T=2000;
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




%% Normal Run

mM_drink=17.4; %mM alcohol per drink
perturb=[200]; %minutes from start, can include multiple

% normal/healthy
decay=.01;
y1=.2;
y2=.25;
[concentration_healthy,velocities_healthy,alcohol,H2O2_healthy]=Patient_specific_function(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
disp('---original structure (healthy)---')
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_healthy)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_healthy)-0
disp("[GSH]_m")
dt*sum(concentration_healthy(:,1))-T*concentration_IC(1)



% diseased
decay=.005;
y1=.5;
y2=.6;

[concentration_diseased,velocities_diseased,alcohol,H2O2_diseased]=Patient_specific_function(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
disp('---original structure (diseased)---')
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_diseased)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_diseased)-0
disp("[GSH]_m")
dt*sum(concentration_diseased(:,1))-T*concentration_IC(1)


% failed compensated
decay=.003;
y1=.75;
y2=.5;

[concentration_failedC,velocities_failedC,alcohol,H2O2_failedC]=Patient_specific_function(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
disp('---original structure (failed, compensated)---')
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_failedC)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_failedC)-0
disp("[GSH]_m")
dt*sum(concentration_failedC(:,1))-T*concentration_IC(1)

% failed decompensated
decay=.001;
y1=.75;
y2=.6;

[concentration_failedD,velocities_failedD,alcohol,H2O2_failedD]=Patient_specific_function(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
disp('---original structure (failed, decompensated)---')
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_failedD)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_failedD)-0
disp("[GSH]_m")
dt*sum(concentration_failedD(:,1))-T*concentration_IC(1)

figure
t=tiledlayout(2,2);
title(t,"Load concentration")
ax1=nexttile;
hold on
plot(ax1,time, concentration_healthy(:,1)./concentration_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax1,time, (concentration_healthy(:,2))./concentration_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax1,time, (concentration_healthy(:,3)),'LineWidth',3,'Color',[0,0,1])
plot(ax1,time, (concentration_healthy(:,4)),'LineWidth',3,'Color',[1,1,0])
plot(ax1,time, (concentration_healthy(:,5)),'LineWidth',3,'Color',[0,1,1])
grid on;
hold off
xlabel(ax1,'time (min)');
ylabel(ax1,'shifted and scaled loads');
title(ax1,'Healthy');
ax2=nexttile;
hold on
plot(ax2,time, concentration_diseased(:,1)./concentration_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax2,time, (concentration_diseased(:,2))./concentration_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax2,time, (concentration_diseased(:,3)),'LineWidth',3,'Color',[0,0,1])
plot(ax2,time, (concentration_diseased(:,4)),'LineWidth',3,'Color',[1,1,0])
plot(ax2,time, (concentration_diseased(:,5)),'LineWidth',3,'Color',[0,1,1])
grid on;
hold off
xlabel(ax2,'time (min)');
ylabel(ax2,'shifted and scaled loads');
title(ax2,'Diseased');
ax3=nexttile;
hold on
plot(ax3,time, concentration_failedC(:,1)./concentration_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax3,time, (concentration_failedC(:,2))./concentration_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax3,time, (concentration_failedC(:,3)),'LineWidth',3,'Color',[0,0,1])
plot(ax3,time, (concentration_failedC(:,4)),'LineWidth',3,'Color',[1,1,0])
plot(ax3,time, (concentration_failedC(:,5)),'LineWidth',3,'Color',[0,1,1])
grid on;
hold off
xlabel(ax3,'time (min)');
ylabel(ax3,'shifted and scaled loads');
title(ax3,'Failed, compensated');
ax4=nexttile;
hold on
plot(ax4,time, concentration_failedD(:,1)./concentration_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax4,time, (concentration_failedD(:,2))./concentration_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax4,time, (concentration_failedD(:,3)),'LineWidth',3,'Color',[0,0,1])
plot(ax4,time, (concentration_failedD(:,4)),'LineWidth',3,'Color',[1,1,0])
plot(ax4,time, (concentration_failedD(:,5)),'LineWidth',3,'Color',[0,1,1])
grid on;
hold off
xlabel(ax4,'time (min)');
ylabel(ax4,'shifted and scaled loads');
title(ax4,'Failed, decompensated');
fontsize(18,"points") 
legend(["[GSH]_m","[SAMe]_c","[ammonia]_m","[a-CoA]_m","[malonyl-CoA]_m"])




figure
t=tiledlayout(2,2);
title(t,"Scaled velocities")
ax1=nexttile;
hold on
plot(ax1,time, velocities_healthy(:,1)./velocities_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax1,time, (velocities_healthy(:,2))./velocities_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax1,time, (velocities_healthy(:,3))./velocities_IC(3),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
plot(ax1,time, (velocities_healthy(:,4))./velocities_IC(4),'LineWidth',3,'Color',[1,1,0])
grid on;
hold off
xlabel(ax1,'time (min)');
ylabel(ax1,'scaled velocities');
title(ax1,'Healthy');
ax2=nexttile;
hold on
plot(ax2,time, velocities_diseased(:,1)./velocities_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax2,time, (velocities_diseased(:,2))./velocities_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax2,time, (velocities_diseased(:,3))./velocities_IC(3),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
plot(ax2,time, (velocities_diseased(:,4))./velocities_IC(4),'LineWidth',3,'Color',[1,1,0])
grid on;
hold off
xlabel(ax2,'time (min)');
ylabel(ax2,'scaled velocities');
title(ax2,'Diseased');
ax3=nexttile;
hold on
plot(ax3,time, velocities_failedC(:,1)./velocities_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax3,time, (velocities_failedC(:,2))./velocities_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax3,time, (velocities_failedC(:,3))./velocities_IC(3),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
plot(ax3,time, (velocities_failedC(:,4))./velocities_IC(4),'LineWidth',3,'Color',[1,1,0])
grid on;
hold off
xlabel(ax3,'time (min)');
ylabel(ax3,'scaled velocities');
title(ax3,'Failed, compensated');
ax4=nexttile;
hold on
plot(ax4,time, velocities_failedD(:,1)./velocities_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax4,time, (velocities_failedD(:,2))./velocities_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax4,time, (velocities_failedD(:,3))./velocities_IC(3),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
plot(ax4,time, (velocities_failedD(:,4))./velocities_IC(4),'LineWidth',3,'Color',[1,1,0])
grid on;
hold off
xlabel(ax4,'time (min)');
ylabel(ax4,'scaled velocities');
title(ax4,'Failed, decompensated');
fontsize(18,"points") 
legend(["GSH","MET","Urea","FABs"])



figure
t=tiledlayout(2,2);
title(t,"[H_2O_2]_m")
ax1=nexttile;
hold on
plot(ax1,time,H2O2_healthy(:).*1000,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax1,'time (min)');
ylabel(ax1,'shifted concentration \mu M');
title(ax1,'Healthy');
ax2=nexttile;
hold on
plot(ax2,time,H2O2_diseased(:).*1000,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax2,'time (min)');
ylabel(ax2,'shifted concentration \mu M');
title(ax2,'Diseased');
ax3=nexttile;
hold on
plot(ax3,time,H2O2_failedC(:).*1000,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax3,'time (min)');
ylabel(ax3,'shifted concentration \mu M');
title(ax3,'Failed, compensated');
ax4=nexttile;
hold on
plot(ax4,time,H2O2_failedD(:).*1000,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax4,'time (min)');
ylabel(ax4,'shifted concentration \mu M');
fontsize(18,"points") 
title(ax4,'Failed, decompensated');



%% Not Normal Run: partial decoupling
%{
ORIGINAL
Adj_mat=[0     7     4     2
         5     0     2     1
         6     5     0     4
         2     1     3     0];
%}

Adj_mat=[0     7     4     0
         5     0     0     0
         6     0     0     4
         0     0     3     0];


mM_drink=17.4; %mM alcohol per drink
perturb=[200]; %minutes from start, can include multiple

% normal/healthy
decay=.01;
y1=.2;
y2=.25;
[concentration_healthy,velocities_healthy,alcohol,H2O2_healthy]=Patient_specific_function(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);

disp('---partial decoupling (healthy)---')
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_healthy)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_healthy)-0
disp("[GSH]_m")
dt*sum(concentration_healthy(:,1))-T*concentration_IC(1)



% diseased
decay=.005;
y1=.5;
y2=.6;

[concentration_diseased,velocities_diseased,alcohol,H2O2_diseased]=Patient_specific_function(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
disp('---partial decoupling (diseased)---')
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_diseased)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_diseased)-0
disp("[GSH]_m")
dt*sum(concentration_diseased(:,1))-T*concentration_IC(1)


% failed compensated
decay=.003;
y1=.75;
y2=.5;

[concentration_failedC,velocities_failedC,alcohol,H2O2_failedC]=Patient_specific_function(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
disp('---partial decoupling (failed, compensated)---')
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_failedC)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_failedC)-0
disp("[GSH]_m")
dt*sum(concentration_failedC(:,1))-T*concentration_IC(1)

% failed decompensated
decay=.001;
y1=.75;
y2=.6;

[concentration_failedD,velocities_failedD,alcohol,H2O2_failedD]=Patient_specific_function(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
disp('---partial decoupling (failed, decompensated)---')
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_failedD)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_failedD)-0
disp("[GSH]_m")
dt*sum(concentration_failedD(:,1))-T*concentration_IC(1)

figure
t=tiledlayout(2,2);
title(t,"Load concentration")
ax1=nexttile;
hold on
plot(ax1,time, concentration_healthy(:,1)./concentration_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax1,time, (concentration_healthy(:,2))./concentration_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax1,time, (concentration_healthy(:,3)),'LineWidth',3,'Color',[0,0,1])
plot(ax1,time, (concentration_healthy(:,4)),'LineWidth',3,'Color',[1,1,0])
plot(ax1,time, (concentration_healthy(:,5)),'LineWidth',3,'Color',[0,1,1])
grid on;
hold off
xlabel(ax1,'time (min)');
ylabel(ax1,'shifted and scaled loads');
title(ax1,'Healthy');
ax2=nexttile;
hold on
plot(ax2,time, concentration_diseased(:,1)./concentration_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax2,time, (concentration_diseased(:,2))./concentration_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax2,time, (concentration_diseased(:,3)),'LineWidth',3,'Color',[0,0,1])
plot(ax2,time, (concentration_diseased(:,4)),'LineWidth',3,'Color',[1,1,0])
plot(ax2,time, (concentration_diseased(:,5)),'LineWidth',3,'Color',[0,1,1])
grid on;
hold off
xlabel(ax2,'time (min)');
ylabel(ax2,'shifted and scaled loads');
title(ax2,'Diseased');
ax3=nexttile;
hold on
plot(ax3,time, concentration_failedC(:,1)./concentration_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax3,time, (concentration_failedC(:,2))./concentration_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax3,time, (concentration_failedC(:,3)),'LineWidth',3,'Color',[0,0,1])
plot(ax3,time, (concentration_failedC(:,4)),'LineWidth',3,'Color',[1,1,0])
plot(ax3,time, (concentration_failedC(:,5)),'LineWidth',3,'Color',[0,1,1])
grid on;
hold off
xlabel(ax3,'time (min)');
ylabel(ax3,'shifted and scaled loads');
title(ax3,'Failed, compensated');
ax4=nexttile;
hold on
plot(ax4,time, concentration_failedD(:,1)./concentration_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax4,time, (concentration_failedD(:,2))./concentration_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax4,time, (concentration_failedD(:,3)),'LineWidth',3,'Color',[0,0,1])
plot(ax4,time, (concentration_failedD(:,4)),'LineWidth',3,'Color',[1,1,0])
plot(ax4,time, (concentration_failedD(:,5)),'LineWidth',3,'Color',[0,1,1])
grid on;
hold off
xlabel(ax4,'time (min)');
ylabel(ax4,'shifted and scaled loads');
title(ax4,'Failed, decompensated');
fontsize(18,"points") 
legend(["[GSH]_m","[SAMe]_c","[ammonia]_m","[a-CoA]_m","[malonyl-CoA]_m"])




figure
t=tiledlayout(2,2);
title(t,"Scaled velocities")
ax1=nexttile;
hold on
plot(ax1,time, velocities_healthy(:,1)./velocities_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax1,time, (velocities_healthy(:,2))./velocities_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax1,time, (velocities_healthy(:,3))./velocities_IC(3),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
plot(ax1,time, (velocities_healthy(:,4))./velocities_IC(4),'LineWidth',3,'Color',[1,1,0])
grid on;
hold off
xlabel(ax1,'time (min)');
ylabel(ax1,'scaled velocities');
title(ax1,'Healthy');
ax2=nexttile;
hold on
plot(ax2,time, velocities_diseased(:,1)./velocities_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax2,time, (velocities_diseased(:,2))./velocities_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax2,time, (velocities_diseased(:,3))./velocities_IC(3),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
plot(ax2,time, (velocities_diseased(:,4))./velocities_IC(4),'LineWidth',3,'Color',[1,1,0])
grid on;
hold off
xlabel(ax2,'time (min)');
ylabel(ax2,'scaled velocities');
title(ax2,'Diseased');
ax3=nexttile;
hold on
plot(ax3,time, velocities_failedC(:,1)./velocities_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax3,time, (velocities_failedC(:,2))./velocities_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax3,time, (velocities_failedC(:,3))./velocities_IC(3),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
plot(ax3,time, (velocities_failedC(:,4))./velocities_IC(4),'LineWidth',3,'Color',[1,1,0])
grid on;
hold off
xlabel(ax3,'time (min)');
ylabel(ax3,'scaled velocities');
title(ax3,'Failed, compensated');
ax4=nexttile;
hold on
plot(ax4,time, velocities_failedD(:,1)./velocities_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax4,time, (velocities_failedD(:,2))./velocities_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax4,time, (velocities_failedD(:,3))./velocities_IC(3),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
plot(ax4,time, (velocities_failedD(:,4))./velocities_IC(4),'LineWidth',3,'Color',[1,1,0])
grid on;
hold off
xlabel(ax4,'time (min)');
ylabel(ax4,'scaled velocities');
title(ax4,'Failed, decompensated');
fontsize(18,"points") 
legend(["GSH","MET","Urea","FABs"])



figure
t=tiledlayout(2,2);
title(t,"[H_2O_2]_m")
ax1=nexttile;
hold on
plot(ax1,time,H2O2_healthy(:).*1000,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax1,'time (min)');
ylabel(ax1,'shifted concentration \mu M');
title(ax1,'Healthy');
ax2=nexttile;
hold on
plot(ax2,time,H2O2_diseased(:).*1000,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax2,'time (min)');
ylabel(ax2,'shifted concentration \mu M');
title(ax2,'Diseased');
ax3=nexttile;
hold on
plot(ax3,time,H2O2_failedC(:).*1000,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax3,'time (min)');
ylabel(ax3,'shifted concentration \mu M');
title(ax3,'Failed, compensated');
ax4=nexttile;
hold on
plot(ax4,time,H2O2_failedD(:).*1000,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax4,'time (min)');
ylabel(ax4,'shifted concentration \mu M');
fontsize(18,"points") 
title(ax4,'Failed, decompensated');


%% Not Normal Run: Full decoupling


Adj_mat=[0     1     0     0
         1     0     0     0
         0     0     0     1
         0     0     1     0];

mM_drink=17.4; %mM alcohol per drink
perturb=[200]; %minutes from start, can include multiple

% normal/healthy
decay=.01;
y1=.2;
y2=.25;
[concentration_healthy,velocities_healthy,alcohol,H2O2_healthy]=Patient_specific_function(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
disp('---full decoupling (healthy)---')
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_healthy)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_healthy)-0
disp("[GSH]_m")
dt*sum(concentration_healthy(:,1))-T*concentration_IC(1)



% diseased
decay=.005;
y1=.5;
y2=.6;

[concentration_diseased,velocities_diseased,alcohol,H2O2_diseased]=Patient_specific_function(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
disp('---full decoupling (diseased)---')
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_diseased)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_diseased)-0
disp("[GSH]_m")
dt*sum(concentration_diseased(:,1))-T*concentration_IC(1)


% failed compensated
decay=.003;
y1=.75;
y2=.5;

[concentration_failedC,velocities_failedC,alcohol,H2O2_failedC]=Patient_specific_function(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
disp('---full decoupling (failed, compensated)---')
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_failedC)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_failedC)-0
disp("[GSH]_m")
dt*sum(concentration_failedC(:,1))-T*concentration_IC(1)

% failed decompensated
decay=.001;
y1=.75;
y2=.6;

[concentration_failedD,velocities_failedD,alcohol,H2O2_failedD]=Patient_specific_function(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
disp('---full decoupling (failed, decompensated)---')
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_failedD)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_failedD)-0
disp("[GSH]_m")
dt*sum(concentration_failedD(:,1))-T*concentration_IC(1)

figure
t=tiledlayout(2,2);
title(t,"Load concentration")
ax1=nexttile;
hold on
plot(ax1,time, concentration_healthy(:,1)./concentration_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax1,time, (concentration_healthy(:,2))./concentration_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax1,time, (concentration_healthy(:,3)),'LineWidth',3,'Color',[0,0,1])
plot(ax1,time, (concentration_healthy(:,4)),'LineWidth',3,'Color',[1,1,0])
plot(ax1,time, (concentration_healthy(:,5)),'LineWidth',3,'Color',[0,1,1])
grid on;
hold off
xlabel(ax1,'time (min)');
ylabel(ax1,'shifted and scaled loads');
title(ax1,'Healthy');
ax2=nexttile;
hold on
plot(ax2,time, concentration_diseased(:,1)./concentration_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax2,time, (concentration_diseased(:,2))./concentration_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax2,time, (concentration_diseased(:,3)),'LineWidth',3,'Color',[0,0,1])
plot(ax2,time, (concentration_diseased(:,4)),'LineWidth',3,'Color',[1,1,0])
plot(ax2,time, (concentration_diseased(:,5)),'LineWidth',3,'Color',[0,1,1])
grid on;
hold off
xlabel(ax2,'time (min)');
ylabel(ax2,'shifted and scaled loads');
title(ax2,'Diseased');
ax3=nexttile;
hold on
plot(ax3,time, concentration_failedC(:,1)./concentration_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax3,time, (concentration_failedC(:,2))./concentration_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax3,time, (concentration_failedC(:,3)),'LineWidth',3,'Color',[0,0,1])
plot(ax3,time, (concentration_failedC(:,4)),'LineWidth',3,'Color',[1,1,0])
plot(ax3,time, (concentration_failedC(:,5)),'LineWidth',3,'Color',[0,1,1])
grid on;
hold off
xlabel(ax3,'time (min)');
ylabel(ax3,'shifted and scaled loads');
title(ax3,'Failed, compensated');
ax4=nexttile;
hold on
plot(ax4,time, concentration_failedD(:,1)./concentration_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax4,time, (concentration_failedD(:,2))./concentration_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax4,time, (concentration_failedD(:,3)),'LineWidth',3,'Color',[0,0,1])
plot(ax4,time, (concentration_failedD(:,4)),'LineWidth',3,'Color',[1,1,0])
plot(ax4,time, (concentration_failedD(:,5)),'LineWidth',3,'Color',[0,1,1])
grid on;
hold off
xlabel(ax4,'time (min)');
ylabel(ax4,'shifted and scaled loads');
title(ax4,'Failed, decompensated');
fontsize(18,"points") 
legend(["[GSH]_m","[SAMe]_c","[ammonia]_m","[a-CoA]_m","[malonyl-CoA]_m"])




figure
t=tiledlayout(2,2);
title(t,"Scaled velocities")
ax1=nexttile;
hold on
plot(ax1,time, velocities_healthy(:,1)./velocities_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax1,time, (velocities_healthy(:,2))./velocities_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax1,time, (velocities_healthy(:,3))./velocities_IC(3),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
plot(ax1,time, (velocities_healthy(:,4))./velocities_IC(4),'LineWidth',3,'Color',[1,1,0])
grid on;
hold off
xlabel(ax1,'time (min)');
ylabel(ax1,'scaled velocities');
title(ax1,'Healthy');
ax2=nexttile;
hold on
plot(ax2,time, velocities_diseased(:,1)./velocities_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax2,time, (velocities_diseased(:,2))./velocities_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax2,time, (velocities_diseased(:,3))./velocities_IC(3),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
plot(ax2,time, (velocities_diseased(:,4))./velocities_IC(4),'LineWidth',3,'Color',[1,1,0])
grid on;
hold off
xlabel(ax2,'time (min)');
ylabel(ax2,'scaled velocities');
title(ax2,'Diseased');
ax3=nexttile;
hold on
plot(ax3,time, velocities_failedC(:,1)./velocities_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax3,time, (velocities_failedC(:,2))./velocities_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax3,time, (velocities_failedC(:,3))./velocities_IC(3),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
plot(ax3,time, (velocities_failedC(:,4))./velocities_IC(4),'LineWidth',3,'Color',[1,1,0])
grid on;
hold off
xlabel(ax3,'time (min)');
ylabel(ax3,'scaled velocities');
title(ax3,'Failed, compensated');
ax4=nexttile;
hold on
plot(ax4,time, velocities_failedD(:,1)./velocities_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax4,time, (velocities_failedD(:,2))./velocities_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax4,time, (velocities_failedD(:,3))./velocities_IC(3),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
plot(ax4,time, (velocities_failedD(:,4))./velocities_IC(4),'LineWidth',3,'Color',[1,1,0])
grid on;
hold off
xlabel(ax4,'time (min)');
ylabel(ax4,'scaled velocities');
title(ax4,'Failed, decompensated');
fontsize(18,"points") 
legend(["GSH","MET","Urea","FABs"])



figure
t=tiledlayout(2,2);
title(t,"[H_2O_2]_m")
ax1=nexttile;
hold on
plot(ax1,time,H2O2_healthy(:).*1000,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax1,'time (min)');
ylabel(ax1,'shifted concentration \mu M');
title(ax1,'Healthy');
ax2=nexttile;
hold on
plot(ax2,time,H2O2_diseased(:).*1000,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax2,'time (min)');
ylabel(ax2,'shifted concentration \mu M');
title(ax2,'Diseased');
ax3=nexttile;
hold on
plot(ax3,time,H2O2_failedC(:).*1000,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax3,'time (min)');
ylabel(ax3,'shifted concentration \mu M');
title(ax3,'Failed, compensated');
ax4=nexttile;
hold on
plot(ax4,time,H2O2_failedD(:).*1000,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax4,'time (min)');
ylabel(ax4,'shifted concentration \mu M');
fontsize(18,"points") 
title(ax4,'Failed, decompensated');


grid on;
hold off
figure
hold on
plot(time, alcohol,'LineWidth',3,'Color',[1, 0, 0])
xlabel('Time (min)');
ylabel('ethanol(mM)');
title('ethanol concentration (all diseased states)');
grid on;



%% Not Normal Run: one way chain

%{
Adj_mat=[0     1     1     0
         1     0     0     0
         1     0     0     1
         0     0     1     0];
%}

Adj_mat=[0     0     0     1
         0     0     1     0
         1     0     0     1
         0     1     0     0];

mM_drink=17.4; %mM alcohol per drink
perturb=[200]; %minutes from start, can include multiple

% normal/healthy
decay=.01;
y1=.2;
y2=.25;
[concentration_healthy,velocities_healthy,alcohol,H2O2_healthy]=Patient_specific_function(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
disp('---full decoupling (healthy)---')
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_healthy)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_healthy)-0
disp("[GSH]_m")
dt*sum(concentration_healthy(:,1))-T*concentration_IC(1)



% diseased
decay=.005;
y1=.5;
y2=.6;

[concentration_diseased,velocities_diseased,alcohol,H2O2_diseased]=Patient_specific_function(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
disp('---full decoupling (diseased)---')
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_diseased)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_diseased)-0
disp("[GSH]_m")
dt*sum(concentration_diseased(:,1))-T*concentration_IC(1)


% failed compensated
decay=.003;
y1=.75;
y2=.5;

[concentration_failedC,velocities_failedC,alcohol,H2O2_failedC]=Patient_specific_function(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
disp('---full decoupling (failed, compensated)---')
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_failedC)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_failedC)-0
disp("[GSH]_m")
dt*sum(concentration_failedC(:,1))-T*concentration_IC(1)

% failed decompensated
decay=.001;
y1=.75;
y2=.6;

[concentration_failedD,velocities_failedD,alcohol,H2O2_failedD]=Patient_specific_function(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
disp('---full decoupling (failed, decompensated)---')
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_failedD)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_failedD)-0
disp("[GSH]_m")
dt*sum(concentration_failedD(:,1))-T*concentration_IC(1)

figure
t=tiledlayout(2,2);
title(t,"Load concentration")
ax1=nexttile;
hold on
plot(ax1,time, concentration_healthy(:,1)./concentration_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax1,time, (concentration_healthy(:,2))./concentration_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax1,time, (concentration_healthy(:,3)),'LineWidth',3,'Color',[0,0,1])
plot(ax1,time, (concentration_healthy(:,4)),'LineWidth',3,'Color',[1,1,0])
plot(ax1,time, (concentration_healthy(:,5)),'LineWidth',3,'Color',[0,1,1])
grid on;
hold off
xlabel(ax1,'time (min)');
ylabel(ax1,'shifted and scaled loads');
title(ax1,'Healthy');
ax2=nexttile;
hold on
plot(ax2,time, concentration_diseased(:,1)./concentration_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax2,time, (concentration_diseased(:,2))./concentration_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax2,time, (concentration_diseased(:,3)),'LineWidth',3,'Color',[0,0,1])
plot(ax2,time, (concentration_diseased(:,4)),'LineWidth',3,'Color',[1,1,0])
plot(ax2,time, (concentration_diseased(:,5)),'LineWidth',3,'Color',[0,1,1])
grid on;
hold off
xlabel(ax2,'time (min)');
ylabel(ax2,'shifted and scaled loads');
title(ax2,'Diseased');
ax3=nexttile;
hold on
plot(ax3,time, concentration_failedC(:,1)./concentration_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax3,time, (concentration_failedC(:,2))./concentration_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax3,time, (concentration_failedC(:,3)),'LineWidth',3,'Color',[0,0,1])
plot(ax3,time, (concentration_failedC(:,4)),'LineWidth',3,'Color',[1,1,0])
plot(ax3,time, (concentration_failedC(:,5)),'LineWidth',3,'Color',[0,1,1])
grid on;
hold off
xlabel(ax3,'time (min)');
ylabel(ax3,'shifted and scaled loads');
title(ax3,'Failed, compensated');
ax4=nexttile;
hold on
plot(ax4,time, concentration_failedD(:,1)./concentration_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax4,time, (concentration_failedD(:,2))./concentration_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax4,time, (concentration_failedD(:,3)),'LineWidth',3,'Color',[0,0,1])
plot(ax4,time, (concentration_failedD(:,4)),'LineWidth',3,'Color',[1,1,0])
plot(ax4,time, (concentration_failedD(:,5)),'LineWidth',3,'Color',[0,1,1])
grid on;
hold off
xlabel(ax4,'time (min)');
ylabel(ax4,'shifted and scaled loads');
title(ax4,'Failed, decompensated');
fontsize(18,"points") 
legend(["[GSH]_m","[SAMe]_c","[ammonia]_m","[a-CoA]_m","[malonyl-CoA]_m"])




figure
t=tiledlayout(2,2);
title(t,"Scaled velocities")
ax1=nexttile;
hold on
plot(ax1,time, velocities_healthy(:,1)./velocities_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax1,time, (velocities_healthy(:,2))./velocities_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax1,time, (velocities_healthy(:,3))./velocities_IC(3),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
plot(ax1,time, (velocities_healthy(:,4))./velocities_IC(4),'LineWidth',3,'Color',[1,1,0])
grid on;
hold off
xlabel(ax1,'time (min)');
ylabel(ax1,'scaled velocities');
title(ax1,'Healthy');
ax2=nexttile;
hold on
plot(ax2,time, velocities_diseased(:,1)./velocities_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax2,time, (velocities_diseased(:,2))./velocities_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax2,time, (velocities_diseased(:,3))./velocities_IC(3),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
plot(ax2,time, (velocities_diseased(:,4))./velocities_IC(4),'LineWidth',3,'Color',[1,1,0])
grid on;
hold off
xlabel(ax2,'time (min)');
ylabel(ax2,'scaled velocities');
title(ax2,'Diseased');
ax3=nexttile;
hold on
plot(ax3,time, velocities_failedC(:,1)./velocities_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax3,time, (velocities_failedC(:,2))./velocities_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax3,time, (velocities_failedC(:,3))./velocities_IC(3),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
plot(ax3,time, (velocities_failedC(:,4))./velocities_IC(4),'LineWidth',3,'Color',[1,1,0])
grid on;
hold off
xlabel(ax3,'time (min)');
ylabel(ax3,'scaled velocities');
title(ax3,'Failed, compensated');
ax4=nexttile;
hold on
plot(ax4,time, velocities_failedD(:,1)./velocities_IC(1),'LineWidth',3,'Color',[1, 0, 0])
plot(ax4,time, (velocities_failedD(:,2))./velocities_IC(2),'LineWidth',3,'Color',[0, 0.4470, 0.7410])
plot(ax4,time, (velocities_failedD(:,3))./velocities_IC(3),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
plot(ax4,time, (velocities_failedD(:,4))./velocities_IC(4),'LineWidth',3,'Color',[1,1,0])
grid on;
hold off
xlabel(ax4,'time (min)');
ylabel(ax4,'scaled velocities');
title(ax4,'Failed, decompensated');
fontsize(18,"points") 
legend(["GSH","MET","Urea","FABs"])



figure
t=tiledlayout(2,2);
title(t,"[H_2O_2]_m")
ax1=nexttile;
hold on
plot(ax1,time,H2O2_healthy(:).*1000,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax1,'time (min)');
ylabel(ax1,'shifted concentration \mu M');
title(ax1,'Healthy');
ax2=nexttile;
hold on
plot(ax2,time,H2O2_diseased(:).*1000,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax2,'time (min)');
ylabel(ax2,'shifted concentration \mu M');
title(ax2,'Diseased');
ax3=nexttile;
hold on
plot(ax3,time,H2O2_failedC(:).*1000,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax3,'time (min)');
ylabel(ax3,'shifted concentration \mu M');
title(ax3,'Failed, compensated');
ax4=nexttile;
hold on
plot(ax4,time,H2O2_failedD(:).*1000,'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax4,'time (min)');
ylabel(ax4,'shifted concentration \mu M');
fontsize(18,"points") 
title(ax4,'Failed, decompensated');


grid on;
hold off
figure
hold on
plot(time, alcohol,'LineWidth',3,'Color',[1, 0, 0])
xlabel('Time (min)');
ylabel('ethanol(mM)');
title('ethanol concentration (all diseased states)');
grid on;







