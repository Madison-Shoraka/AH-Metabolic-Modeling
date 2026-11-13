
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
metabolic_decomp_case1=zeros(resolution,4);
total_oxidative_stress_case1=zeros(resolution,1);
epsilon_t=100;
epsilon_t=T/dt-epsilon_t/dt;


mM_drink=17.4; %mM alcohol per drink
perturb=[200]; %minutes from start, can include multiple
parfor j=1:resolution
        decay=decay_list(j)
        y1=y1_list(j)
        y2=y2_list(j)        
       [concentration,velocities,alcohol,ROS,total_oxidative_stress]=full_model(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
        total_concentration_case1(j,:)=dt*sum(concentration)-T*concentration_IC;
        total_velocities_case1(j,:)=dt*sum(velocities)-T*velocities_IC;
        total_H2O2_case1(j)=dt*sum(ROS)-0;
        total_GSH_case1(j)=dt*sum(concentration(:,1))-T*concentration_IC(1);
        metabolic_decomp_case1(j,:)=dt*sum(abs(velocities(epsilon_t:end,:)-velocities_IC));
        total_oxidative_stress_case1(j)=dt*sum(total_oxidative_stress);
end 
   

%case 2
total_concentration_case2=zeros(resolution,5);
total_velocities_case2=zeros(resolution,4);
total_H2O2_case2=zeros(resolution,1);
total_GSH_case2=zeros(resolution,1);
metabolic_decomp_case2=zeros(resolution,4);
total_oxidative_stress_case2=zeros(resolution,1);


mM_drink=17.4; %mM alcohol per drink
perturb=[200,1640]; %minutes from start, can include multiple (24 hours later)
parfor j=1:resolution
        decay=decay_list(j)
        y1=y1_list(j)
        y2=y2_list(j)
       [concentration,velocities,alcohol,ROS,total_oxidative_stress]=full_model(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
        total_concentration_case2(j,:)=dt*sum(concentration)-T*concentration_IC;
        total_velocities_case2(j,:)=dt*sum(velocities)-T*velocities_IC;
        total_H2O2_case2(j)=dt*sum(ROS)-0;
        total_GSH_case2(j)=dt*sum(concentration(:,1))-T*concentration_IC(1);
        metabolic_decomp_case2(j,:)=dt*sum(abs(velocities(epsilon_t:end,:)-velocities_IC))
        total_oxidative_stress_case2(j)=dt*sum(total_oxidative_stress);


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
cycleRegulation = zeros(2, 6); % 2 cases, 6 metrics
for i = 1:2
    cycleRegulation(i, 1) = length(find(eval(['total_velocities_case' num2str(i) '(:,3)<0'])))/resolution; % Urea cycle downregulated
    cycleRegulation(i, 2) = length(find(eval(['total_velocities_case' num2str(i) '(:,2)<0'])))/resolution; % MET cycle downregulated
    cycleRegulation(i, 3) = length(find(eval(['total_velocities_case' num2str(i) '(:,4)>0'])))/resolution; % FABs cycle upregulated
    cycleRegulation(i, 4) = length(find(eval(['total_concentration_case' num2str(i) '(:,1)<0'])))/resolution; % total GSH_m depeletion 
    cycleRegulation(i, 5) = length(find(eval(['total_H2O2_case' num2str(i) '>0'])))/resolution; % H2O2 outcompetes GSH
    cycleRegulation(i, 6) = length(find(any(eval(['metabolic_decomp_case' num2str(i) '>.1'])')'))/resolution; % H2O2 outcompetes GSH
    
end
disp(cycleRegulation);


    oxidative_failure_case1=(total_H2O2_case1 >0); 
    met_decomp_case1=(any((metabolic_decomp_case1>.1)')');


    %comparing
    oxidative_failure_AND_met_failure1=oxidative_failure_case1 & met_decomp_case1;
    oxidative_failure_NO_met_failure1=oxidative_failure_case1 & ~met_decomp_case1;
    met_failure_NO_oxidative_failure1=~oxidative_failure_case1 & met_decomp_case1;
    
    oxidative_failure_case2=(total_H2O2_case2 >0); 
    met_decomp_case2=(any((metabolic_decomp_case2>1)')');
    
    index_failure_case1=unique([find(oxidative_failure_AND_met_failure1);find(met_failure_NO_oxidative_failure1)]);


    %comparing
    oxidative_failure_AND_met_failure2=oxidative_failure_case2 & met_decomp_case2;
    
    oxidative_failure_NO_met_failure2=oxidative_failure_case2 & ~met_decomp_case2;
    
    met_failure_NO_oxidative_failure2=~oxidative_failure_case2 & met_decomp_case2;
    
    index_failure_case2=unique([find(oxidative_failure_AND_met_failure2);find(met_failure_NO_oxidative_failure2)]);


Failure_mat=[sum(oxidative_failure_AND_met_failure1)/resolution, sum(oxidative_failure_NO_met_failure1)/resolution, sum(met_failure_NO_oxidative_failure1)/resolution;
    sum(oxidative_failure_AND_met_failure2)/resolution, sum(oxidative_failure_NO_met_failure2)/resolution, sum(met_failure_NO_oxidative_failure2)/resolution]


%% Data plotting for failure types


figure()
hold on
t = tiledlayout(3,4);
title(t,"BOTH metabolic decomp and oxidative failure by decay value (k_d) one-hit")

for i = 1:length(decay_vec)
    ax = nexttile;
    
    idx = find(abs(decay_list - decay_vec(i)) < 1e-12);
    
    slice = reshape(oxidative_failure_AND_met_failure1(idx), 11, 11);
    corresponding_y1 = reshape(y1_list(idx), 11, 11);
    corresponding_y2 = reshape(y2_list(idx), 11, 11);
    
    
    imagesc(corresponding_y2(1,:),corresponding_y1(:,1), (slice))
    set(ax, 'YDir', 'normal'); % Set to normal (bottom-to-top)
    title(ax, sprintf("k_d = %g", decay_vec(i)))
    xlabel(ax, "y_2")
    ylabel(ax, "y_1")
    colorbar()
    clim([0, 1])
end
fontsize(gcf, "scale", 1.32);
hold off

figure()
hold on
t = tiledlayout(3,4);
title(t,"ONLY metabolic decomp by decay value (k_d) one-hit")

for i = 1:length(decay_vec)
    ax = nexttile;
    
    idx = find(abs(decay_list - decay_vec(i)) < 1e-12);
    
    slice = reshape(met_failure_NO_oxidative_failure1(idx), 11, 11);
    corresponding_y1 = reshape(y1_list(idx), 11, 11);
    corresponding_y2 = reshape(y2_list(idx), 11, 11);
    
    
    imagesc(corresponding_y2(1,:),corresponding_y1(:,1), (slice))
    set(ax, 'YDir', 'normal'); % Set to normal (bottom-to-top)
    title(ax, sprintf("k_d = %g", decay_vec(i)))
    xlabel(ax, "y_2")
    ylabel(ax, "y_1")
    colorbar()
    clim([0, 1])
end
fontsize(gcf, "scale", 1.32);
hold off



figure()
hold on
t = tiledlayout(3,4);
title(t,"ONLY oxidative failure by decay value (k_d) one-hit")

for i = 1:length(decay_vec)
    ax = nexttile;
    
    idx = find(abs(decay_list - decay_vec(i)) < 1e-12);
    
    slice = reshape(oxidative_failure_NO_met_failure1(idx), 11, 11);
    corresponding_y1 = reshape(y1_list(idx), 11, 11);
    corresponding_y2 = reshape(y2_list(idx), 11, 11);
    
    
    imagesc(corresponding_y2(1,:),corresponding_y1(:,1), (slice))
    set(ax, 'YDir', 'normal'); % Set to normal (bottom-to-top)
    title(ax, sprintf("k_d = %g", decay_vec(i)))
    xlabel(ax, "y_2")
    ylabel(ax, "y_1")
    colorbar()
    clim([0, 1])
end
fontsize(gcf, "scale", 1.32);
hold off

figure()
hold on
t = tiledlayout(3,4);
title(t,"BOTH metabolic decomp and oxidative failure by decay value (k_d) two-hit")

for i = 1:length(decay_vec)
    ax = nexttile;
    
    idx = find(abs(decay_list - decay_vec(i)) < 1e-12);
    
    slice = reshape(oxidative_failure_AND_met_failure2(idx), 11, 11);
    corresponding_y1 = reshape(y1_list(idx), 11, 11);
    corresponding_y2 = reshape(y2_list(idx), 11, 11);
    
    
    imagesc(corresponding_y2(1,:),corresponding_y1(:,1), (slice))
    set(ax, 'YDir', 'normal'); % Set to normal (bottom-to-top)
    title(ax, sprintf("k_d = %g", decay_vec(i)))
    xlabel(ax, "y_2")
    ylabel(ax, "y_1")
    colorbar()
    clim([0, 1])
end
fontsize(gcf, "scale", 1.32);
hold off

figure()
hold on
t = tiledlayout(3,4);
title(t,"ONLY metabolic decomp by decay value (k_d) two-hit")

for i = 1:length(decay_vec)
    ax = nexttile;
    
    idx = find(abs(decay_list - decay_vec(i)) < 1e-12);
    
    slice = reshape(met_failure_NO_oxidative_failure2(idx), 11, 11);
    corresponding_y1 = reshape(y1_list(idx), 11, 11);
    corresponding_y2 = reshape(y2_list(idx), 11, 11);
    
    
    imagesc(corresponding_y2(1,:),corresponding_y1(:,1), (slice))
    set(ax, 'YDir', 'normal'); % Set to normal (bottom-to-top)
    title(ax, sprintf("k_d = %g", decay_vec(i)))
    xlabel(ax, "y_2")
    ylabel(ax, "y_1")
    colorbar()
    clim([0, 1])
end
fontsize(gcf, "scale", 1.32);
hold off



figure()
hold on
t = tiledlayout(3,4);
title(t,"ONLY oxidative failure by decay value (k_d) two-hit")

for i = 1:length(decay_vec)
    ax = nexttile;
    
    idx = find(abs(decay_list - decay_vec(i)) < 1e-12);
    
    slice = reshape(oxidative_failure_NO_met_failure2(idx), 11, 11);
    corresponding_y1 = reshape(y1_list(idx), 11, 11);
    corresponding_y2 = reshape(y2_list(idx), 11, 11);
    
    
    imagesc(corresponding_y2(1,:),corresponding_y1(:,1), (slice))
    set(ax, 'YDir', 'normal'); % Set to normal (bottom-to-top)
    title(ax, sprintf("k_d = %g", decay_vec(i)))
    xlabel(ax, "y_2")
    ylabel(ax, "y_1")
    colorbar()
    clim([0, 1])
end
fontsize(gcf, "scale", 1.32);
hold off
