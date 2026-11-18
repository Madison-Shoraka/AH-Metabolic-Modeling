
%% standard parameters
y3=1; %Time delay (min) between alcohol metabolism and acetyl-CoA production
tau=10000;
dt=.005;
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


%velocities_IC_list=[general_cubic_poly(health_state_list,-0.0002,0.0012,0.025,1.9659) general_cubic_poly(health_state_list, -0.0002,0.003,-0.0126,0.4253) general_cubic_poly(health_state_list,0.0006,-0.0061,-0.0113,0.7215) general_cubic_poly(health_state_list, -0.0002,0.002,-0.0011 ,0.8874)];

%% Run cases

%case 1
total_concentration_case1=zeros(resolution,5);
total_velocities_case1=zeros(resolution,4);
total_H2O2_case1=zeros(resolution,1);
metabolic_decomp_case1=zeros(resolution,4);
mM_drink=17.4; %mM alcohol per drink
perturb=[200]; %minutes from start, can include multiple
T=perturb(end)+3*14 40;
epsilon_t=100;
epsilon_t=T/dt-epsilon_t/dt;
time=0:dt:T-dt;
nt=length(time);
parfor j=1:resolution
        decay=decay_list(j)
        y1=y1_list(j)
        y2=y2_list(j)
       [concentration,velocities,alcohol,H2O2]=full_model(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
        total_concentration_case1(j,:)=dt*sum(concentration)-T*concentration_IC;
        total_velocities_case1(j,:)=dt*sum(velocities)-T*velocities_IC;
        total_H2O2_case1(j)=dt*sum(H2O2)-0;
        metabolic_decomp_case1(j,:)=dt*sum(abs(velocities(epsilon_t:end,:)-velocities_IC));

end 
   

%case 2
total_concentration_case2=zeros(resolution,5);
total_velocities_case2=zeros(resolution,4);
total_H2O2_case2=zeros(resolution,1);
metabolic_decomp_case2=zeros(resolution,4);
mM_drink=17.4; %mM alcohol per drink
perturb=[200,260]; %minutes from start, can include multiple
T=perturb(end)+3*1440;
epsilon_t=100;
epsilon_t=T/dt-epsilon_t/dt;
time=0:dt:T-dt;
nt=length(time);
parfor j=1:resolution
        decay=decay_list(j)
        y1=y1_list(j)
        y2=y2_list(j)
       
       [concentration,velocities,alcohol,H2O2]=full_model(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
        total_concentration_case2(j,:)=dt*sum(concentration)-T*concentration_IC;
        total_velocities_case2(j,:)=dt*sum(velocities)-T*velocities_IC;
        total_H2O2_case2(j)=dt*sum(H2O2)-0;
        metabolic_decomp_case2(j,:)=dt*sum(abs(velocities(epsilon_t:end,:)-velocities_IC));
end 


%case 3
total_concentration_case3=zeros(resolution,5);
total_velocities_case3=zeros(resolution,4);
total_H2O2_case3=zeros(resolution,1);
metabolic_decomp_case3=zeros(resolution,4);
mM_drink=17.4; %mM alcohol per drink
perturb=[200,260,1640]; %minutes from start, can include multiple
T=perturb(end)+3*1440;
epsilon_t=100;
epsilon_t=T/dt-epsilon_t/dt;
time=0:dt:T-dt;
nt=length(time);
parfor j=1:resolution
        decay=decay_list(j)
        y1=y1_list(j)
        y2=y2_list(j)
        
       [concentration,velocities,alcohol,H2O2]=full_model(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
        total_concentration_case3(j,:)=dt*sum(concentration)-T*concentration_IC;
        total_velocities_case3(j,:)=dt*sum(velocities)-T*velocities_IC;
        total_H2O2_case3(j)=dt*sum(H2O2)-0;
        metabolic_decomp_case3(j,:)=dt*sum(abs(velocities(epsilon_t:end,:)-velocities_IC));
end 

%case 4
total_concentration_case4=zeros(resolution,5);
total_velocities_case4=zeros(resolution,4);
total_H2O2_case4=zeros(resolution,1);
metabolic_decomp_case4=zeros(resolution,4);
mM_drink=17.4; %mM alcohol per drink
perturb=[200,260,1640,1700]; %minutes from start, can include multiple
T=perturb(end)+3*1440;
epsilon_t=100;
epsilon_t=T/dt-epsilon_t/dt;
time=0:dt:T-dt;
nt=length(time);
parfor j=1:resolution
        decay=decay_list(j)
        y1=y1_list(j)
        y2=y2_list(j)
       
       [concentration,velocities,alcohol,H2O2]=full_model(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
       total_concentration_case4(j,:)=dt*sum(concentration)-T*concentration_IC;
        total_velocities_case4(j,:)=dt*sum(velocities)-T*velocities_IC;
        total_H2O2_case4(j)=dt*sum(H2O2)-0;
        metabolic_decomp_case4(j,:)=dt*sum(abs(velocities(epsilon_t:end,:)-velocities_IC));
end

%% Plot results
     

%case 1

figure()
t=tiledlayout(2,2);
title(t,"Net velocity change case 1")
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
title(t,"Net velocity changes case 2")
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

%case 3
figure()
t=tiledlayout(2,2);
title(t,"Net velocity changes case 3")
ax1 = nexttile;
scatter3(ax1,decay_list,y2_list,y1_list,40,total_velocities_case3(:,1),'filled')
title(ax1, "GSH cycle")
xlabel(ax1,"k_d")
ylabel(ax1,"y_2")
zlabel(ax1,"y_1")
colorbar()
clim([-max(abs(total_velocities_case3(:,1))),max(abs(total_velocities_case3(:,1)))])
ax2 = nexttile;
scatter3(ax2,decay_list,y2_list,y1_list,40,total_velocities_case3(:,2),'filled')
title(ax2, "MET cycle")
xlabel(ax2,"k_d")
ylabel(ax2,"y_2")
zlabel(ax2,"y_1")
colorbar()
clim([-max(abs(total_velocities_case3(:,2))),max(abs(total_velocities_case3(:,2)))])
ax3 = nexttile;
scatter3(ax3,decay_list,y2_list,y1_list,40,total_velocities_case3(:,3),'filled')
title(ax3, "Urea cycle")
xlabel(ax3,"k_d")
ylabel(ax3,"y_2")
zlabel(ax3,"y_1")
colorbar()
clim([-max(abs(total_velocities_case3(:,3))),max(abs(total_velocities_case3(:,3)))])
ax4 = nexttile;
scatter3(ax4,decay_list,y2_list,y1_list,40,total_velocities_case3(:,4),'filled')
title(ax4, "FABs cycle")
xlabel(ax4,"k_d")
ylabel(ax4,"y_2")
zlabel(ax4,"y_1")
colorbar()
clim([-max(abs(total_velocities_case3(:,4))),max(abs(total_velocities_case3(:,4)))])
linkaxes([ax1,ax3],'x');
linkaxes([ax1,ax2],'y');
linkaxes([ax2,ax4],'x');
linkaxes([ax3,ax4],'y');


%case 4
figure()
t=tiledlayout(2,2);
title(t,"Net velocity changes case 4")
ax1 = nexttile;
scatter3(ax1,decay_list,y2_list,y1_list,40,total_velocities_case4(:,1),'filled')
title(ax1, "GSH cycle")
xlabel(ax1,"k_d")
ylabel(ax1,"y_2")
zlabel(ax1,"y_1")
colorbar()
clim([-max(abs(total_velocities_case4(:,1))),max(abs(total_velocities_case4(:,1)))])
ax2 = nexttile;
scatter3(ax2,decay_list,y2_list,y1_list,40,total_velocities_case4(:,2),'filled')
title(ax2, "MET cycle")
xlabel(ax2,"k_d")
ylabel(ax2,"y_2")
zlabel(ax2,"y_1")
colorbar()
clim([-max(abs(total_velocities_case4(:,2))),max(abs(total_velocities_case4(:,2)))])
ax3 = nexttile;
scatter3(ax3,decay_list,y2_list,y1_list,40,total_velocities_case4(:,3),'filled')
title(ax3, "Urea cycle")
xlabel(ax3,"k_d")
ylabel(ax3,"y_2")
zlabel(ax3,"y_1")
colorbar()
clim([-max(abs(total_velocities_case4(:,3))),max(abs(total_velocities_case4(:,3)))])
ax4 = nexttile;
scatter3(ax4,decay_list,y2_list,y1_list,40,total_velocities_case4(:,4),'filled')
title(ax4, "FABs cycle")
xlabel(ax4,"k_d")
ylabel(ax4,"y_2")
zlabel(ax4,"y_1")
colorbar()
clim([-max(abs(total_velocities_case4(:,4))),max(abs(total_velocities_case4(:,4)))])
linkaxes([ax1,ax3],'x');
linkaxes([ax1,ax2],'y');
linkaxes([ax2,ax4],'x');
linkaxes([ax3,ax4],'y');


%case 1
figure()
t=tiledlayout(2,3);
title(t,"Net concentration changes case 1")
ax1 = nexttile;
scatter3(ax1,decay_list,y2_list,y1_list,40,total_concentration_case1(:,1),'filled')
title(ax1, "[GSH]_m")
xlabel(ax1,"k_d")
ylabel(ax1,"y_2")
zlabel(ax1,"y_1")
colorbar()
clim([-max(abs(total_concentration_case1(:,1))),max(abs(total_concentration_case1(:,1)))])
ax2 = nexttile;
scatter3(ax2,decay_list,y2_list,y1_list,40,total_concentration_case1(:,2),'filled')
title(ax2, "[SAM]_c")
xlabel(ax2,"k_d")
ylabel(ax2,"y_2")
zlabel(ax2,"y_1")
colorbar()
clim([-max(abs(total_concentration_case1(:,2))),max(abs(total_concentration_case1(:,2)))])
ax3 = nexttile;
scatter3(ax3,decay_list,y2_list,y1_list,40,total_concentration_case1(:,3),'filled')
title(ax3, "[ammonia]_m")
xlabel(ax3,"k_d")
ylabel(ax3,"y_2")
zlabel(ax3,"y_1")
colorbar()
clim([-max(abs(total_concentration_case1(:,3))),max(abs(total_concentration_case1(:,3)))])
ax4 = nexttile;
scatter3(ax4,decay_list,y2_list,y1_list,40,total_concentration_case1(:,4),'filled')
title(ax4, "[aCoA]_m")
xlabel(ax4,"k_d")
ylabel(ax4,"y_2")
zlabel(ax4,"y_1")
colorbar()
clim([-max(abs(total_concentration_case1(:,4))),max(abs(total_concentration_case1(:,4)))])
ax5 = nexttile;
scatter3(ax5,decay_list,y2_list,y1_list,40,total_concentration_case1(:,5),'filled')
title(ax5, "[malonyl-CoA]_m")
xlabel(ax5,"k_d")
ylabel(ax5,"y_2")
zlabel(ax5,"y_1")
colorbar()
clim([-max(abs(total_concentration_case1(:,5))),max(abs(total_concentration_case1(:,5)))])
linkaxes([ax1,ax3],'x');
linkaxes([ax1,ax2],'y');
linkaxes([ax2,ax4],'x');
linkaxes([ax3,ax4],'y');

%case 2
figure()
t=tiledlayout(2,3);
title(t,"Net concentration changes case 2")
ax1 = nexttile;
scatter3(ax1,decay_list,y2_list,y1_list,40,total_concentration_case2(:,1),'filled')
title(ax1, "[GSH]_m")
xlabel(ax1,"k_d")
ylabel(ax1,"y_2")
zlabel(ax1,"y_1")
colorbar()
clim([-max(abs(total_concentration_case2(:,1))),max(abs(total_concentration_case2(:,1)))])
ax2 = nexttile;
scatter3(ax2,decay_list,y2_list,y1_list,40,total_concentration_case2(:,2),'filled')
title(ax2, "[SAM]_c")
xlabel(ax2,"k_d")
ylabel(ax2,"y_2")
zlabel(ax2,"y_1")
colorbar()
clim([-max(abs(total_concentration_case2(:,2))),max(abs(total_concentration_case2(:,2)))])
ax3 = nexttile;
scatter3(ax3,decay_list,y2_list,y1_list,40,total_concentration_case2(:,3),'filled')
title(ax3, "[ammonia]_m")
xlabel(ax3,"k_d")
ylabel(ax3,"y_2")
zlabel(ax3,"y_1")
colorbar()
clim([-max(abs(total_concentration_case2(:,3))),max(abs(total_concentration_case2(:,3)))])
ax4 = nexttile;
scatter3(ax4,decay_list,y2_list,y1_list,40,total_concentration_case2(:,4),'filled')
title(ax4, "[aCoA]_m")
xlabel(ax4,"k_d")
ylabel(ax4,"y_2")
zlabel(ax4,"y_1")
colorbar()
clim([-max(abs(total_concentration_case2(:,4))),max(abs(total_concentration_case2(:,4)))])
ax5 = nexttile;
scatter3(ax5,decay_list,y2_list,y1_list,40,total_concentration_case2(:,5),'filled')
title(ax5, "[malonyl-CoA]_m")
xlabel(ax5,"k_d")
ylabel(ax5,"y_2")
zlabel(ax5,"y_1")
colorbar()
clim([-max(abs(total_concentration_case2(:,5))),max(abs(total_concentration_case2(:,5)))])
linkaxes([ax1,ax3],'x');
linkaxes([ax1,ax2],'y');
linkaxes([ax2,ax4],'x');
linkaxes([ax3,ax4],'y');

%case 3
figure()
t=tiledlayout(2,3);
title(t,"Net concentration changes case 3")
ax1 = nexttile;
scatter3(ax1,decay_list,y2_list,y1_list,40,total_concentration_case3(:,1),'filled')
title(ax1, "[GSH]_m")
xlabel(ax1,"k_d")
ylabel(ax1,"y_2")
zlabel(ax1,"y_1")
colorbar()
clim([-max(abs(total_concentration_case3(:,1))),max(abs(total_concentration_case3(:,1)))])
ax2 = nexttile;
scatter3(ax2,decay_list,y2_list,y1_list,40,total_concentration_case3(:,2),'filled')
title(ax2, "[SAM]_c")
xlabel(ax2,"k_d")
ylabel(ax2,"y_2")
zlabel(ax2,"y_1")
colorbar()
clim([-max(abs(total_concentration_case3(:,2))),max(abs(total_concentration_case3(:,2)))])
ax3 = nexttile;
scatter3(ax3,decay_list,y2_list,y1_list,40,total_concentration_case3(:,3),'filled')
title(ax3, "[ammonia]_m")
xlabel(ax3,"k_d")
ylabel(ax3,"y_2")
zlabel(ax3,"y_1")
colorbar()
clim([-max(abs(total_concentration_case3(:,3))),max(abs(total_concentration_case3(:,3)))])
ax4 = nexttile;
scatter3(ax4,decay_list,y2_list,y1_list,40,total_concentration_case3(:,4),'filled')
title(ax4, "[aCoA]_m")
xlabel(ax4,"k_d")
ylabel(ax4,"y_2")
zlabel(ax4,"y_1")
colorbar()
clim([-max(abs(total_concentration_case3(:,4))),max(abs(total_concentration_case3(:,4)))])
ax5 = nexttile;
scatter3(ax5,decay_list,y2_list,y1_list,40,total_concentration_case3(:,5),'filled')
title(ax5, "[malonyl-CoA]_m")
xlabel(ax5,"k_d")
ylabel(ax5,"y_2")
zlabel(ax5,"y_1")
colorbar()
clim([-max(abs(total_concentration_case3(:,5))),max(abs(total_concentration_case3(:,5)))])
linkaxes([ax1,ax3],'x');
linkaxes([ax1,ax2],'y');
linkaxes([ax2,ax4],'x');
linkaxes([ax3,ax4],'y');


%case 4
figure()
t=tiledlayout(2,3);
title(t,"Net concentration changes case 4")
ax1 = nexttile;
scatter3(ax1,decay_list,y2_list,y1_list,40,total_concentration_case4(:,1),'filled')
title(ax1, "[GSH]_m")
xlabel(ax1,"k_d")
ylabel(ax1,"y_2")
zlabel(ax1,"y_1")
colorbar()
clim([-max(abs(total_concentration_case4(:,1))),max(abs(total_concentration_case4(:,1)))])
ax2 = nexttile;
scatter3(ax2,decay_list,y2_list,y1_list,40,total_concentration_case4(:,2),'filled')
title(ax2, "[SAM]_c")
xlabel(ax2,"k_d")
ylabel(ax2,"y_2")
zlabel(ax2,"y_1")
colorbar()
clim([-max(abs(total_concentration_case4(:,2))),max(abs(total_concentration_case4(:,2)))])
ax3 = nexttile;
scatter3(ax3,decay_list,y2_list,y1_list,40,total_concentration_case4(:,3),'filled')
title(ax3, "[ammonia]_m")
xlabel(ax3,"k_d")
ylabel(ax3,"y_2")
zlabel(ax3,"y_1")
colorbar()
clim([-max(abs(total_concentration_case4(:,3))),max(abs(total_concentration_case4(:,3)))])
ax4 = nexttile;
scatter3(ax4,decay_list,y2_list,y1_list,40,total_concentration_case4(:,4),'filled')
title(ax4, "[aCoA]_m")
xlabel(ax4,"k_d")
ylabel(ax4,"y_2")
zlabel(ax4,"y_1")
colorbar()
clim([-max(abs(total_concentration_case4(:,4))),max(abs(total_concentration_case4(:,4)))])
ax5 = nexttile;
scatter3(ax5,decay_list,y2_list,y1_list,40,total_concentration_case3(:,5),'filled')
title(ax5, "[malonyl-CoA]_m")
xlabel(ax5,"k_d")
ylabel(ax5,"y_2")
zlabel(ax5,"y_1")
colorbar()
clim([-max(abs(total_concentration_case4(:,5))),max(abs(total_concentration_case4(:,5)))])
linkaxes([ax1,ax3],'x');
linkaxes([ax1,ax2],'y');
linkaxes([ax2,ax4],'x');
linkaxes([ax3,ax4],'y');


%H2O2
figure()
t=tiledlayout(2,2);
title(t,"Net H2O2 changes")
ax1 = nexttile;
scatter3(ax1,decay_list,y2_list,y1_list,40,total_H2O2_case1.*1000,'filled')
title(ax1, "case 1")
xlabel(ax1,"k_d")
ylabel(ax1,"y_2")
zlabel(ax1,"y_1")
colorbar()
clim([-max(abs(total_H2O2_case1)),max(abs(total_H2O2_case1))])
ax2 = nexttile;
scatter3(ax2,decay_list,y2_list,y1_list,40,total_H2O2_case2.*1000,'filled')
title(ax2, "case 2")
xlabel(ax2,"k_d")
ylabel(ax2,"y_2")
zlabel(ax2,"y_1")
colorbar()
clim([-max(abs(total_H2O2_case2)),max(abs(total_H2O2_case2))])
ax3 = nexttile;
scatter3(ax3,decay_list,y2_list,y1_list,40,total_H2O2_case3.*1000,'filled')
title(ax3, "case 3")
xlabel(ax3,"k_d")
ylabel(ax3,"y_2")
zlabel(ax3,"y_1")
colorbar()
clim([-max(abs(total_H2O2_case3)),max(abs(total_H2O2_case3))])
ax4 = nexttile;
scatter3(ax4,decay_list,y2_list,y1_list,40,total_H2O2_case4.*1000,'filled')
title(ax4, "case 4")
xlabel(ax4,"k_d")
ylabel(ax4,"y_2")
zlabel(ax4,"y_1")
colorbar()
clim([-max(abs(total_H2O2_case4)),max(abs(total_H2O2_case4))])
linkaxes([ax1,ax3],'x');
linkaxes([ax1,ax2],'y');
linkaxes([ax2,ax4],'x');
linkaxes([ax3,ax4],'y');


%% Data processing

% Process relevant cases for H2O2

x1=total_H2O2_case1.*1000>1;
breakCases_case1=[y1_list(x1),y2_list(x1),decay_list(x1)];
x2=total_H2O2_case2.*1000>1;
breakCases_case2=[y1_list(x2),y2_list(x2),decay_list(x2)];
x3=total_H2O2_case3.*1000>1;
breakCases_case3=[y1_list(x3),y2_list(x3),decay_list(x3)];
x4=total_H2O2_case4.*1000>1;
breakCases_case4=[y1_list(x4),y2_list(x4),decay_list(x4)];
first_break=4-max([4*x1,3*x2,2*x3,x4]')';
decay_vec=decay_min:(decay_max-decay_min)/m:decay_max;



figure()
hold on

for i = 1:length(decay_vec)
    ax = nexttile;
    
    idx = find(abs(decay_list - decay_vec(i)) < 1e-12);
    
    slice = reshape(first_break(idx), 11, 11);
    corresponding_y1 = reshape(y1_list(idx), 11, 11);
    corresponding_y2 = reshape(y2_list(idx), 11, 11);
    
    
    imagesc(corresponding_y2(1,:),corresponding_y1(:,1), (slice))
    set(gca, 'YDir', 'normal')
    axis off
    
    clim([0,4])
end
hold off

%% break analysis
y1_vec=y1_min:.2:y1_max;
y2_vec=y2_min:.2:y1_max;
decay_vec=decay_min:.001:decay_max;
[X,Y,Z]=meshgrid(y1_vec,y2_vec,decay_vec);
areas=zeros(length(y1_vec),length(y2_vec),length(decay_vec));
for i=1:length(y1_vec)-1
    for j=1:length(y2_vec)-1
        for k=1:length(decay_vec)-1
            y1_found=((y1_list>y1_vec(i)).*(y1_list<=y1_vec(i+1)));
            y2_found=((y2_list>y2_vec(j)).*(y2_list<=y2_vec(j+1)));
            decay_found=((decay_list>decay_vec(k)).*(decay_list<=decay_vec(k+1)));
            if sum(y1_found.*y2_found.*decay_found)>0
                areas(i,j,k)=mean(first_break(find(y1_found.*y2_found.*decay_found)));
            end
        end 
    end
end 
%Cycles emergent bahivor 
cycleRegulation = zeros(4, 5); % 4 cases, 4 metrics
for i = 1:4
    cycleRegulation(i, 1) = length(find(eval(['total_velocities_case' num2str(i) '(:,3)<0'])))/resolution; % Urea cycle downregulated
    cycleRegulation(i, 2) = length(find(eval(['total_velocities_case' num2str(i) '(:,1)>0'])))/resolution; % GSH cycle upregulated
    cycleRegulation(i, 3) = length(find(eval(['total_velocities_case' num2str(i) '(:,2)<0'])))/resolution; % MET cycle downregulated
    cycleRegulation(i, 4) = length(find(eval(['total_velocities_case' num2str(i) '(:,4)>0'])))/resolution; % FABs cycle upregulated
     cycleRegulation(i, 5) = length(find(eval(['total_H2O2_case' num2str(i) '>0'])))/resolution; % H2O2 outcompetes GSH
     cycleRegulation(i, 6) = length(find(max(eval(['metabolic_decomp_case' num2str(i) '>.1'])')'))/resolution; % metabolic decomp
     cycleRegulation(i, 7) = length(find(max(eval(['metabolic_decomp_case' num2str(i) '>.1'])')'.*eval(['total_H2O2_case' num2str(i) '>0'])))/resolution; % metabolic decomp AND H2O2

end
disp(cycleRegulation);


%% slice images
decay_vec=decay_min:(decay_max-decay_min)/m:decay_max;

figure()
t = tiledlayout(3,4);
title(t,"Net H2O2 changes per decay value (k_d)")

for i = 1:length(decay_vec)
    ax = nexttile;
    
    idx = find(abs(decay_list - decay_vec(i)) < 1e-12);
    
    slice = reshape(total_H2O2_case1(idx), 11, 11);
    corresponding_y1 = reshape(y1_list(idx), 11, 11);
    corresponding_y2 = reshape(y2_list(idx), 11, 11);
    
    imagesc( corresponding_y2(1,:),corresponding_y1(:,1), (slice))
    set(ax, 'YDir', 'normal'); % Set to normal (bottom-to-top)
    title(ax, sprintf("k_d = %g", decay_vec(i)))
    xlabel(ax, "y_2")
    ylabel(ax, "y_1")
    colorbar()
end

figure()
t = tiledlayout(3,4);
title(t,"Net H2O2 changes per decay value (k_d)")

for i = 1:length(decay_vec)
    ax = nexttile;
    
    idx = find(abs(decay_list - decay_vec(i)) < 1e-12);
    
    slice = reshape(total_H2O2_case2(idx), 11, 11);
    corresponding_y1 = reshape(y1_list(idx), 11, 11);
    corresponding_y2 = reshape(y2_list(idx), 11, 11);
    
    imagesc( corresponding_y2(1,:),corresponding_y1(:,1), (slice))
    set(ax, 'YDir', 'normal'); % Set to normal (bottom-to-top)
    title(ax, sprintf("k_d = %g", decay_vec(i)))
    xlabel(ax, "y_2")
    ylabel(ax, "y_1")
    colorbar()
end

figure()
t = tiledlayout(3,4);
title(t,"Net H2O2 changes per decay value (k_d)")

for i = 1:length(decay_vec)
    ax = nexttile;
    
    idx = find(abs(decay_list - decay_vec(i)) < 1e-12);
    
    slice = reshape(total_H2O2_case3(idx), 11, 11);
    corresponding_y1 = reshape(y1_list(idx), 11, 11);
    corresponding_y2 = reshape(y2_list(idx), 11, 11);
    
    imagesc( corresponding_y2(1,:),corresponding_y1(:,1), (slice))
    set(ax, 'YDir', 'normal'); % Set to normal (bottom-to-top)
    title(ax, sprintf("k_d = %g", decay_vec(i)))
    xlabel(ax, "y_2")
    ylabel(ax, "y_1")
    colorbar()
    
end

figure()
t = tiledlayout(3,4);
title(t,"Net H2O2 changes per decay value (k_d)")

for i = 1:length(decay_vec)
    ax = nexttile;
    
    idx = find(abs(decay_list - decay_vec(i)) < 1e-12);
    
    slice = reshape(total_H2O2_case4(idx), 11, 11);
    corresponding_y1 = reshape(y1_list(idx), 11, 11);
    corresponding_y2 = reshape(y2_list(idx), 11, 11);
    
    imagesc( corresponding_y2(1,:),corresponding_y1(:,1), (slice))
    set(ax, 'YDir', 'normal'); % Set to normal (bottom-to-top)
    title(ax, sprintf("k_d = %g", decay_vec(i)))
    xlabel(ax, "y_2")
    ylabel(ax, "y_1")
    colorbar()
    
end


%% Check
figure()
t=tiledlayout(2,2);
title(t,"Net velocity change case 1")
ax1 = nexttile;
scatter3(ax1,decay_list,y2_list,y1_list,40,total_velocities_case1(:,1)>0,'filled')
title(ax1, "GSH cycle")
xlabel(ax1,"k_d")
ylabel(ax1,"y_2")
zlabel(ax1,"y_1")
colorbar()
ax2 = nexttile;
scatter3(ax2,decay_list,y2_list,y1_list,40,total_velocities_case1(:,2)<0,'filled')
title(ax2, "MET cycle")
xlabel(ax2,"k_d")
ylabel(ax2,"y_2")
zlabel(ax2,"y_1")
colorbar()
ax3 = nexttile;
scatter3(ax3,decay_list,y2_list,y1_list,40,total_velocities_case1(:,3)<0,'filled')
title(ax3, "Urea cycle")
xlabel(ax3,"k_d")
ylabel(ax3,"y_2")
zlabel(ax3,"y_1")
colorbar()
ax4 = nexttile;
scatter3(ax4,decay_list,y2_list,y1_list,40,total_velocities_case1(:,4)>0,'filled')
title(ax4, "FABs cycle")
xlabel(ax4,"k_d")
ylabel(ax4,"y_2")
zlabel(ax4,"y_1")
colorbar()
linkaxes([ax1,ax3],'x');
linkaxes([ax1,ax2],'y');
linkaxes([ax2,ax4],'x');
linkaxes([ax3,ax4],'y');

figure()
t=tiledlayout(2,2);
title(t,"Net velocity change case 1")
ax1 = nexttile;
scatter3(ax1,decay_list,y2_list,y1_list,40,total_velocities_case2(:,1)>0,'filled')
title(ax1, "GSH cycle")
xlabel(ax1,"k_d")
ylabel(ax1,"y_2")
zlabel(ax1,"y_1")
colorbar()
ax2 = nexttile;
scatter3(ax2,decay_list,y2_list,y1_list,40,total_velocities_case2(:,2)<0,'filled')
title(ax2, "MET cycle")
xlabel(ax2,"k_d")
ylabel(ax2,"y_2")
zlabel(ax2,"y_1")
colorbar()
ax3 = nexttile;
scatter3(ax3,decay_list,y2_list,y1_list,40,total_velocities_case2(:,3)<0,'filled')
title(ax3, "Urea cycle")
xlabel(ax3,"k_d")
ylabel(ax3,"y_2")
zlabel(ax3,"y_1")
colorbar()
ax4 = nexttile;
scatter3(ax4,decay_list,y2_list,y1_list,40,total_velocities_case2(:,4)>0,'filled')
title(ax4, "FABs cycle")
xlabel(ax4,"k_d")
ylabel(ax4,"y_2")
zlabel(ax4,"y_1")
colorbar()
linkaxes([ax1,ax3],'x');
linkaxes([ax1,ax2],'y');
linkaxes([ax2,ax4],'x');
linkaxes([ax3,ax4],'y');

figure()
t=tiledlayout(2,2);
title(t,"Net velocity change case 1")
ax1 = nexttile;
scatter3(ax1,decay_list,y2_list,y1_list,40,total_velocities_case3(:,1)>0,'filled')
title(ax1, "GSH cycle")
xlabel(ax1,"k_d")
ylabel(ax1,"y_2")
zlabel(ax1,"y_1")
colorbar()
ax2 = nexttile;
scatter3(ax2,decay_list,y2_list,y1_list,40,total_velocities_case3(:,2)<0,'filled')
title(ax2, "MET cycle")
xlabel(ax2,"k_d")
ylabel(ax2,"y_2")
zlabel(ax2,"y_1")
colorbar()
ax3 = nexttile;
scatter3(ax3,decay_list,y2_list,y1_list,40,total_velocities_case3(:,3)<0,'filled')
title(ax3, "Urea cycle")
xlabel(ax3,"k_d")
ylabel(ax3,"y_2")
zlabel(ax3,"y_1")
colorbar()
ax4 = nexttile;
scatter3(ax4,decay_list,y2_list,y1_list,40,total_velocities_case3(:,4)>0,'filled')
title(ax4, "FABs cycle")
xlabel(ax4,"k_d")
ylabel(ax4,"y_2")
zlabel(ax4,"y_1")
colorbar()
linkaxes([ax1,ax3],'x');
linkaxes([ax1,ax2],'y');
linkaxes([ax2,ax4],'x');
linkaxes([ax3,ax4],'y');

figure()
t=tiledlayout(2,2);
title(t,"Net velocity change case 1")
ax1 = nexttile;
scatter3(ax1,decay_list,y2_list,y1_list,40,total_velocities_case4(:,1)>0,'filled')
title(ax1, "GSH cycle")
xlabel(ax1,"k_d")
ylabel(ax1,"y_2")
zlabel(ax1,"y_1")
colorbar()
ax2 = nexttile;
scatter3(ax2,decay_list,y2_list,y1_list,40,total_velocities_case4(:,2)<0,'filled')
title(ax2, "MET cycle")
xlabel(ax2,"k_d")
ylabel(ax2,"y_2")
zlabel(ax2,"y_1")
colorbar()
ax3 = nexttile;
scatter3(ax3,decay_list,y2_list,y1_list,40,total_velocities_case4(:,3)<0,'filled')
title(ax3, "Urea cycle")
xlabel(ax3,"k_d")
ylabel(ax3,"y_2")
zlabel(ax3,"y_1")
colorbar()
ax4 = nexttile;
scatter3(ax4,decay_list,y2_list,y1_list,40,total_velocities_case4(:,4)>0,'filled')
title(ax4, "FABs cycle")
xlabel(ax4,"k_d")
ylabel(ax4,"y_2")
zlabel(ax4,"y_1")
colorbar()
linkaxes([ax1,ax3],'x');
linkaxes([ax1,ax2],'y');
linkaxes([ax2,ax4],'x');
linkaxes([ax3,ax4],'y');




