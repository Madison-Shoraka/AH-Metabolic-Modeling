
%% standard parameters
y3=1; %Time delay (min) between alcohol metabolism and acetyl-CoA production
tau=10000;
dt=.005;
T=1640+2*1440;
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

%% Set up grid
decay_min=.001;
decay_max=.02;
y1_min=.01;
y1_max=.75;
y2_min=.1;
y2_max=1;

dimension=3;
m=10;
resolution=(m+1)^dimension;


decay_vec=decay_min:(decay_max-decay_min)/m:decay_max;
y1_vec=y1_min:(y1_max-y1_min)/m:y1_max;
y2_vec=y2_min:(y2_max-y2_min)/m:y2_max;

[decay_list,y1_list,y2_list]=meshgrid(decay_vec,y1_vec,y2_vec);
decay_list=reshape(decay_list,resolution,1);
y1_list=reshape(y1_list,resolution,1);
y2_list=reshape(y2_list,resolution,1);


%% Run cases

%case 1
total_concentration_case1=zeros(resolution,5);
total_velocities_case1=zeros(resolution,4);
total_H2O2_case1=zeros(resolution,1);
total_GSH_case1=zeros(resolution,1);


mM_drink=17.4; %mM alcohol per drink
perturb=[200]; %minutes from start, can include multiple
parfor j=1:resolution
        decay=decay_list(j)
        y1=y1_list(j)
        y2=y2_list(j)
        
       [concentration,velocities,alcohol,H2O2,total_GSH]=full_model(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
        total_concentration_case1(j,:)=dt*sum(concentration)-T*concentration_IC;
        total_velocities_case1(j,:)=dt*sum(velocities)-T*velocities_IC;
        total_H2O2_case1(j)=dt*sum(H2O2)-0;
        
        total_GSH_case1(j)=dt*sum(total_GSH)-T*concentration_IC(1);
end 
   

%case 2
total_concentration_case2=zeros(resolution,5);
total_velocities_case2=zeros(resolution,4);
total_H2O2_case2=zeros(resolution,1);
total_GSH_case2=zeros(resolution,1);


mM_drink=17.4; %mM alcohol per drink
perturb=[200,1640]; %minutes from start, can include multiple (24 hours later)
parfor j=1:resolution
        decay=decay_list(j)
        y1=y1_list(j)
        y2=y2_list(j)
       [concentration,velocities,alcohol,H2O2,total_GSH]=full_model(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
        total_concentration_case2(j,:)=dt*sum(concentration)-T*concentration_IC;
        total_velocities_case2(j,:)=dt*sum(velocities)-T*velocities_IC;
        total_H2O2_case2(j)=dt*sum(H2O2)-0;
        total_GSH_case2(j)=dt*sum(total_GSH)-T*concentration_IC(1);
end 



%% Plot results
     

%case 1
figure()
t=tiledlayout(2,2);
title(t,"Net velocity change one-hit")
ax1 = nexttile;
scatter3(ax1,decay_list,y2_list,y1_list,40,total_velocities_case1(:,1),'filled')
title(ax1, "GSH cycle")
xlabel(ax1,"k_d")
ylabel(ax1,"y_2")
zlabel(ax1,"y_1")
colorbar()
clim([-max(abs(total_velocities_case1(:,1))),max(abs(total_velocities_case1(:,1)))])
ax2 = nexttile;
scatter3(ax2,decay_list,y2_list,y1_list,40,total_velocities_case1(:,2),'filled')
title(ax2, "MET cycle")
xlabel(ax2,"k_d")
ylabel(ax2,"y_2")
zlabel(ax2,"y_1")
colorbar()
clim([-max(abs(total_velocities_case1(:,2))),max(abs(total_velocities_case1(:,2)))])
ax3 = nexttile;
scatter3(ax3,decay_list,y2_list,y1_list,40,total_velocities_case1(:,3),'filled')
title(ax3, "Urea cycle")
xlabel(ax3,"k_d")
ylabel(ax3,"y_2")
zlabel(ax3,"y_1")
colorbar()
clim([-max(abs(total_velocities_case1(:,3))),max(abs(total_velocities_case1(:,3)))])
ax4 = nexttile;
scatter3(ax4,decay_list,y2_list,y1_list,40,total_velocities_case1(:,4),'filled')
title(ax4, "FABs cycle")
xlabel(ax4,"k_d")
ylabel(ax4,"y_2")
zlabel(ax4,"y_1")
colorbar()
clim([-max(abs(total_velocities_case1(:,4))),max(abs(total_velocities_case1(:,4)))])
linkaxes([ax1,ax3],'x');
linkaxes([ax1,ax2],'y');
linkaxes([ax2,ax4],'x');
linkaxes([ax3,ax4],'y');

%case 2
figure()
t=tiledlayout(2,2);
title(t,"Net velocity changes two-hit")
ax1 = nexttile;
scatter3(ax1,decay_list,y2_list,y1_list,40,total_velocities_case2(:,1),'filled')
title(ax1, "GSH cycle")
xlabel(ax1,"k_d")
ylabel(ax1,"y_2")
zlabel(ax1,"y_1")
colorbar()
clim([-max(abs(total_velocities_case2(:,1))),max(abs(total_velocities_case2(:,1)))])
ax2 = nexttile;
scatter3(ax2,decay_list,y2_list,y1_list,40,total_velocities_case2(:,2),'filled')
title(ax2, "MET cycle")
xlabel(ax2,"k_d")
ylabel(ax2,"y_2")
zlabel(ax2,"y_1")
colorbar()
clim([-max(abs(total_velocities_case2(:,2))),max(abs(total_velocities_case2(:,2)))])
ax3 = nexttile;
scatter3(ax3,decay_list,y2_list,y1_list,40,total_velocities_case2(:,3),'filled')
title(ax3, "Urea cycle")
xlabel(ax3,"k_d")
ylabel(ax3,"y_2")
zlabel(ax3,"y_1")
colorbar()
clim([-max(abs(total_velocities_case2(:,3))),max(abs(total_velocities_case2(:,3)))])
ax4 = nexttile;
scatter3(ax4,decay_list,y2_list,y1_list,40,total_velocities_case2(:,4),'filled')
title(ax4, "FABs cycle")
xlabel(ax4,"k_d")
ylabel(ax4,"y_2")
zlabel(ax4,"y_1")
colorbar()
clim([-max(abs(total_velocities_case2(:,4))),max(abs(total_velocities_case2(:,4)))])
linkaxes([ax1,ax3],'x');
linkaxes([ax1,ax2],'y');
linkaxes([ax2,ax4],'x');
linkaxes([ax3,ax4],'y');


%H2O2
figure()
t=tiledlayout(2,1);
title(t,"Net H2O2 changes")
ax1 = nexttile;
scatter3(ax1,decay_list,y2_list,y1_list,40,total_H2O2_case1.*1000,'filled')
title(ax1, "one-hit")
xlabel(ax1,"k_d")
ylabel(ax1,"y_2")
zlabel(ax1,"y_1")
colorbar()
clim([-max(abs(total_H2O2_case1)),max(abs(total_H2O2_case1))])
ax2 = nexttile;
scatter3(ax2,decay_list,y2_list,y1_list,40,total_H2O2_case2.*1000,'filled')
title(ax2, "two-hit")
xlabel(ax2,"k_d")
ylabel(ax2,"y_2")
zlabel(ax2,"y_1")
colorbar()
clim([-max(abs(total_H2O2_case2)),max(abs(total_H2O2_case2))])


%GSH
figure()
t=tiledlayout(2,1);
title(t,"Net GSH changes")
ax1 = nexttile;
scatter3(ax1,decay_list,y2_list,y1_list,40,total_GSH_case1,'filled')
title(ax1, "one-hit")
xlabel(ax1,"k_d")
ylabel(ax1,"y_2")
zlabel(ax1,"y_1")
colorbar()
clim([-max(abs(total_GSH_case1)),max(abs(total_GSH_case1))])
ax2 = nexttile;
scatter3(ax2,decay_list,y2_list,y1_list,40,total_GSH_case2,'filled')
title(ax2, "two-hit")
xlabel(ax2,"k_d")
ylabel(ax2,"y_2")
zlabel(ax2,"y_1")
colorbar()
clim([-max(abs(total_GSH_case2)),max(abs(total_GSH_case2))])


%% Data processing


%Cycles emergent bahivor 
cycleRegulation = zeros(2, 6); % 4 cases, 4 metrics
for i = 1:2
    cycleRegulation(i, 1) = length(find(eval(['total_velocities_case' num2str(i) '(:,3)<0'])))/resolution; % Urea cycle downregulated
    cycleRegulation(i, 3) = length(find(eval(['total_velocities_case' num2str(i) '(:,2)<0'])))/resolution; % MET cycle downregulated
    cycleRegulation(i, 4) = length(find(eval(['total_velocities_case' num2str(i) '(:,4)>0'])))/resolution; % FABs cycle upregulated
    cycleRegulation(i, 5) = length(find(eval(['total_concentration_case' num2str(i) '(:,1)<0'])))/resolution; % total GSH_m depeletion 
    cycleRegulation(i, 6) = length(find(eval(['total_H2O2_case' num2str(i) '>0'])))/resolution; % H2O2 outcompetes GSH
end
disp(cycleRegulation);








%% assess instability
%unstable cases
unstable_param_case1=[decay_list(find(total_H2O2_case1>10)),y1_list(find(total_H2O2_case1>10)),y2_list(find(total_H2O2_case1>10))];
unstable_param_case2=[decay_list(find(total_H2O2_case2>10)),y1_list(find(total_H2O2_case2>10)),y2_list(find(total_H2O2_case2>10))];

% Assessing the stability of the system based on unstable parameters
unstableCount_case1 = size(unstable_param_case1, 1);
unstableCount_case2 = size(unstable_param_case2, 1);
disp(['Unstable cases for one-hit: ', num2str(unstableCount_case1)]);
disp(['Unstable cases for two-hit: ', num2str(unstableCount_case2)]);

%assess systems that are unstable
indicies_unstable=find(total_H2O2_case1>10);
res_refined=length(indicies_unstable);
%Cycles emergent bahivor 
cycleRegulation_unstable = zeros(1, 6); % 4 cases, 4 metrics
for i = 1:1
    cycleRegulation_unstable(i, 1) = length(find(eval(['total_velocities_case' num2str(i) '(indicies_unstable,3)<0'])))/res_refined; % Urea cycle downregulated
    cycleRegulation_unstable(i, 3) = length(find(eval(['total_velocities_case' num2str(i) '(indicies_unstable,2)<0'])))/res_refined; % MET cycle downregulated
    cycleRegulation_unstable(i, 4) = length(find(eval(['total_velocities_case' num2str(i) '(indicies_unstable,4)>0'])))/res_refined; % FABs cycle upregulated
    cycleRegulation_unstable(i, 5) = length(find(eval(['total_concentration_case' num2str(i) '(indicies_unstable,1)<0'])))/res_refined; % total GSH_m depeletion 
    cycleRegulation_unstable(i, 6) = length(find(eval(['total_H2O2_case' num2str(i) '(indicies_unstable)>0'])))/res_refined; % H2O2 outcompetes GSH
end
disp(cycleRegulation_unstable);


%% assess conditional stability

indicies_stable=find(total_H2O2_case1<10);
res_refined=length(indicies_stable);

%Cycles emergent bahivor 
cycleRegulation_stable = zeros(1, 5); % 1 case, 4 metrics
for i = 1:1
    cycleRegulation_stable(i, 1) = length(find(eval(['total_velocities_case' num2str(i) '(indicies_stable,3)<0'])))/res_refined; % Urea cycle downregulated
    cycleRegulation_stable(i, 2) = length(find(eval(['total_velocities_case' num2str(i) '(indicies_stable,1)>0'])))/res_refined; % GSH cycle upregulated
    cycleRegulation_stable(i, 3) = length(find(eval(['total_velocities_case' num2str(i) '(indicies_stable,2)<0'])))/res_refined; % MET cycle downregulated
    cycleRegulation_stable(i, 4) = length(find(eval(['total_velocities_case' num2str(i) '(indicies_stable,4)>0'])))/res_refined; % FABs cycle upregulated
    cycleRegulation_stable(i, 5) = length(find(eval(['total_concentration_case' num2str(i) '(indicies_stable,1)<0'])))/res_refined; % total GSH_m depeletion 
end
disp(cycleRegulation_stable);


%% assess always stability
%test if there are y1,y2,decay parameters that are always stable
setxor(y1_list,y1_list(find(total_H2O2_case1>0)));
decay_stable=setxor(decay_list,decay_list(find(total_H2O2_case1>0)));
setxor(y2_list,y2_list(find(total_H2O2_case1>0)));

%assess systems that are always stable
indicies_always_stable=find(ismember(decay_list,decay_stable));
res_refined=length(indicies_always_stable);
%Cycles emergent bahivor 
cycleRegulation_always_stable = zeros(1, 6); % 4 cases, 4 metrics
for i = 1:1
    cycleRegulation_always_stable(i, 1) = length(find(eval(['total_velocities_case' num2str(i) '(indicies_always_stable,3)<0'])))/res_refined; % Urea cycle downregulated
    cycleRegulation_always_stable(i, 2) = length(find(eval(['total_velocities_case' num2str(i) '(indicies_always_stable,1)>0'])))/res_refined; % GSH cycle upregulated
    cycleRegulation_always_stable(i, 3) = length(find(eval(['total_velocities_case' num2str(i) '(indicies_always_stable,2)<0'])))/res_refined; % MET cycle downregulated
    cycleRegulation_always_stable(i, 4) = length(find(eval(['total_velocities_case' num2str(i) '(indicies_always_stable,4)>0'])))/res_refined; % FABs cycle upregulated
    cycleRegulation_always_stable(i, 5) = length(find(eval(['total_concentration_case' num2str(i) '(indicies_always_stable,1)<0'])))/res_refined; % total GSH_m depeletion 
    cycleRegulation_always_stable(i, 6) = length(find(eval(['total_H2O2_case' num2str(i) '(indicies_always_stable)>0'])))/res_refined; % H2O2 outcompetes GSH
end
disp(cycleRegulation_always_stable);




%% Examples
mM_drink=17.4; %mM alcohol per drink
perturb=[200]; %minutes from start, can include multiple

% normal/healthy
decay=.01;
y1=.2;
y2=.25;
[concentration_healthy,velocities_healthy,alcohol,H2O2_healthy,total_GSH_healthy]=full_model(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);

disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_healthy)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_healthy)-0
disp("[GSH]_m^T")
dt*sum(total_GSH_healthy)-T*concentration_IC(1)



% diseased
decay=.005;
y1=.5;
y2=.6;

[concentration_diseased,velocities_diseased,alcohol,H2O2_diseased,total_GSH_diseased]=full_model(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_diseased)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_diseased)-0
disp("[GSH]_m^T")
dt*sum(total_GSH_diseased)-T*concentration_IC(1)


% failed compensated
decay=.003;
y1=.75;
y2=.5;

[concentration_failedC,velocities_failedC,alcohol,H2O2_failedC,total_GSH_failedC]=full_model(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_failedC)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_failedC)-0
disp("[GSH]_m^T")
dt*sum(total_GSH_failedC)-T*concentration_IC(1)

% failed decompensated
decay=.001;
y1=.75;
y2=.6;

[concentration_failedD,velocities_failedD,alcohol,H2O2_failedD,total_GSH_failedD]=full_model(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_failedD)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_failedD)-0
disp("[GSH]_m^T")
dt*sum(total_GSH_failedD)-T*concentration_IC(1)



figure
t=tiledlayout(2,2);
title(t,"Concentration Burden per Metabolite")
ax1=nexttile;
hold on
plot(ax1,time, (concentration_IC(1)-concentration_healthy(:,1))./(critical_concentration_vec(1)),'LineWidth',3)
plot(ax1,time, (concentration_IC(2)-concentration_healthy(:,2))./(critical_concentration_vec(2)),'LineWidth',3)
plot(ax1,time, (concentration_healthy(:,3))./critical_concentration_vec(3),'LineWidth',3)
plot(ax1,time, (concentration_healthy(:,4))./critical_concentration_vec(4),'LineWidth',3)
plot(ax1,time, (concentration_healthy(:,5))./critical_concentration_vec(5),'LineWidth',3)
grid on;
hold off
xlabel(ax1,'time (min)');
ylabel(ax1,'');
title(ax1,'Healthy');
ax2=nexttile;
hold on
plot(ax2,time, (concentration_IC(1)-concentration_diseased(:,1))./(critical_concentration_vec(1)),'LineWidth',3)
plot(ax2,time, (concentration_IC(2)-concentration_diseased(:,2))./(critical_concentration_vec(2)),'LineWidth',3)
plot(ax2,time, (concentration_diseased(:,3))./critical_concentration_vec(3),'LineWidth',3)
plot(ax2,time, (concentration_diseased(:,4))./critical_concentration_vec(4),'LineWidth',3)
plot(ax2,time, (concentration_diseased(:,5))./critical_concentration_vec(5),'LineWidth',3)
grid on;
hold off
xlabel(ax2,'time (min)');
ylabel(ax2,'shifted and scaled loads');
title(ax2,'Diseased');
ax3=nexttile;
hold on
plot(ax3,time, (concentration_IC(1)-concentration_failedC(:,1))./(critical_concentration_vec(1)),'LineWidth',3)
plot(ax3,time, (concentration_IC(2)-concentration_failedC(:,2))./(critical_concentration_vec(2)),'LineWidth',3)
plot(ax3,time, (concentration_failedC(:,3))./critical_concentration_vec(3),'LineWidth',3)
plot(ax3,time, (concentration_failedC(:,4))./critical_concentration_vec(4),'LineWidth',3)
plot(ax3,time, (concentration_failedC(:,5))./critical_concentration_vec(5),'LineWidth',3)
grid on;
hold off
xlabel(ax3,'time (min)');
ylabel(ax3,'shifted and scaled loads');
title(ax3,'Failed, compensated');
ax4=nexttile;
hold on
plot(ax4,time, (concentration_IC(1)-concentration_failedD(:,1))./(critical_concentration_vec(1)),'LineWidth',3)
plot(ax4,time, (concentration_IC(2)-concentration_failedD(:,2))./(critical_concentration_vec(2)),'LineWidth',3)
plot(ax4,time, (concentration_failedD(:,3))./critical_concentration_vec(3),'LineWidth',3)
plot(ax4,time, (concentration_failedD(:,4))./critical_concentration_vec(4),'LineWidth',3)
plot(ax4,time, (concentration_failedD(:,5))./critical_concentration_vec(5),'LineWidth',3)
grid on;
hold off
xlabel(ax4,'time (min)');
ylabel(ax4,'shifted and scaled loads');
title(ax4,'Failed, decompensated');
fontsize(18,"points") 
legend(["[GSH]_{mito}","[SAMe]_{cyto}","[ammonia]_{mito}","[a-CoA]_{mito}","[malonyl-CoA]_{mito}"])




figure
t=tiledlayout(2,2);
title(t,"Scaled velocities")
ax1=nexttile;
hold on
plot(ax1,time, velocities_healthy(:,1)./velocities_IC(1),'LineWidth',3)
plot(ax1,time, (velocities_healthy(:,2))./velocities_IC(2),'LineWidth',3)
plot(ax1,time, (velocities_healthy(:,3))./velocities_IC(3),'LineWidth',3)
plot(ax1,time, (velocities_healthy(:,4))./velocities_IC(4),'LineWidth',3)
grid on;
hold off
xlabel(ax1,'time (min)');
ylabel(ax1,'scaled velocities');
title(ax1,'Healthy');
ax2=nexttile;
hold on
plot(ax2,time, velocities_diseased(:,1)./velocities_IC(1),'LineWidth',3)
plot(ax2,time, (velocities_diseased(:,2))./velocities_IC(2),'LineWidth',3)
plot(ax2,time, (velocities_diseased(:,3))./velocities_IC(3),'LineWidth',3)
plot(ax2,time, (velocities_diseased(:,4))./velocities_IC(4),'LineWidth',3)
grid on;
hold off
xlabel(ax2,'time (min)');
ylabel(ax2,'scaled velocities');
title(ax2,'Diseased');
ax3=nexttile;
hold on
plot(ax3,time, velocities_failedC(:,1)./velocities_IC(1),'LineWidth',3)
plot(ax3,time, (velocities_failedC(:,2))./velocities_IC(2),'LineWidth',3)
plot(ax3,time, (velocities_failedC(:,3))./velocities_IC(3),'LineWidth',3)
plot(ax3,time, (velocities_failedC(:,4))./velocities_IC(4),'LineWidth',3)
grid on;
hold off
xlabel(ax3,'time (min)');
ylabel(ax3,'scaled velocities');
title(ax3,'Failed, compensated');
ax4=nexttile;
hold on
plot(ax4,time, velocities_failedD(:,1)./velocities_IC(1),'LineWidth',3)
plot(ax4,time, (velocities_failedD(:,2))./velocities_IC(2),'LineWidth',3)
plot(ax4,time, (velocities_failedD(:,3))./velocities_IC(3),'LineWidth',3)
plot(ax4,time, (velocities_failedD(:,4))./velocities_IC(4),'LineWidth',3)
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
ylabel(t,"[H_2O_2]_m (mM)")
xlabel(t,"time (min)")
ax1=nexttile;
hold on
plot(ax1,time,H2O2_healthy(:),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
title(ax1,'Healthy');
ax2=nexttile;
hold on
plot(ax2,time,H2O2_diseased(:),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
title(ax2,'Diseased');
ax3=nexttile;
hold on
plot(ax3,time,H2O2_failedC(:),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
xlabel(ax3,'time (min)');
ylabel(ax3,'mM');
title(ax3,'Failed, compensated');
ax4=nexttile;
hold on
plot(ax4,time,H2O2_failedD(:),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
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
hold off



%% Example: Metabolic decomp without oxidative failure
mM_drink=17.4; %mM alcohol per drink
perturb=[200]; %minutes from start, can include multiple


decay=.001;
y1=.6;
y2=.1;

[concentration_met,velocities_met,alcohol,H2O2_met,total_GSH_met]=full_model(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
disp("Net velocity changes (GSH,MET,Urea,FABs)")
dt*sum(velocities_met)-T*velocities_IC
disp("Net H_2O_2 accumulation")
dt*sum(H2O2_met)-0
disp("[GSH]_m^T")
dt*sum(total_GSH_met)-T*concentration_IC(1)


figure

title("Concentration Burden per Metabolite")
hold on
plot(time, (concentration_IC(1)-concentration_met(:,1))./(critical_concentration_vec(1)),'LineWidth',3)
plot(time, (concentration_IC(2)-concentration_met(:,2))./(critical_concentration_vec(2)),'LineWidth',3)
plot(time, (concentration_met(:,3))./critical_concentration_vec(3),'LineWidth',3)
plot(time, (concentration_met(:,4))./critical_concentration_vec(4),'LineWidth',3)
plot(time, (concentration_met(:,5))./critical_concentration_vec(5),'LineWidth',3)
grid on;
hold off

figure
title("[H_2O_2]_m")
ylabel("[H_2O_2]_m (mM)")
xlabel("time (min)")
hold on
plot(time,H2O2_met(:),'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250])
grid on;
hold off
