
%% BASE PARAMETERS
%drinking perturbations
mM_drink=17.4; %mM alcohol per drink
perturb=[150,390]; %minutes from start, can include multiple

%standard parameters
tau=100000;
y3=1; %Time delay (min) between alcohol metabolism and acetyl-CoA production
T=4000; %min
dt=.01;
time=0:dt:T-dt;
sigma=3;
time_to_detox=60;
velocities_IC= [2.0737    0.4002    0.6193    0.9068]; %GSH,MET,Urea,FABs
l=[1 0 0 0 %metabolite = GSH_m
   0 1 0 0 %metbolite = Met_c
   0 0 1 0 %metabolite = ammonia_m
   0 0 1 1 %metabolite = a-CoA_m
   0 0 0 1]; %metabolite=malonyl-CoA_m
type=[1 1 0 0 0]; %1=produced, 0=used


concentration_IC=[10 .08 0 0 0];
critical_concentration_vec=[7.5,.04,1,2,.3];
critical_concentrations=[1 0 0 0 %metabolite = GSH_m
    0 1 0 0 %metbolite = SAM_c
    0 0 1 0 %metabolite = ammonia_m
    0 0 1/2 1/2 %metabolite = a-CoA_m
    0 0 0 1]; %metabolite = malonyl-CoA_m 
critical_concentrations=critical_concentration_vec'.*critical_concentrations;
max_concentration=[concentration_IC(1) concentration_IC(2) 100 100 100];

beta = 7; % sharpness parameter
alcohol_met=-.05;
Adj_mat=[0     7     4     2
     5     0     2     1
     6     5     0     4
     2     1     3     0];

response_rate=.01;
    


%% UQ Setup 
dimension=3;
pts=sobolset(2*dimension);
m=10;
resolution=2^m;
samplePoints = net(pts, resolution);
decay_min=.001;
decay_max=.02;
y1_min=.01;
y1_max=.75;
y2_min=.1;
y2_max=1;


matrix_A=samplePoints(:,1:dimension);
decay_listA=decay_min+(decay_max-decay_min).*matrix_A(:,1);
y1_listA=y1_min+(y1_max-y1_min).*matrix_A(:,2);
y2_listA=y2_min+(y2_max-y2_min).*matrix_A(:,3);
matrix_B=samplePoints(:,(dimension+1):(2*dimension));
decay_listB=decay_min+(decay_max-decay_min).*matrix_B(:,1);
y1_listB=y1_min+(y1_max-y1_min).*matrix_B(:,2);
y2_listB=y2_min+(y2_max-y2_min).*matrix_B(:,3);


%% Running UQ
sensitivity_indicies_concentrations=zeros(dimension,length(concentration_IC));
sensitivity_indicies_velocities=zeros(dimension,length(velocities_IC));
total_effect_indicies_concentrations=zeros(dimension,length(concentration_IC));
total_effect_indicies_velocities=zeros(dimension,length(velocities_IC));
sensitivity_indicies_H2O2=zeros(dimension,1);
total_effect_indicies_H2O2=zeros(dimension,1);
total_concentration_A=zeros(resolution,length(concentration_IC));
total_concentration_B=zeros(resolution,length(concentration_IC));
total_velocity_A=zeros(resolution,length(velocities_IC));
total_velocity_B=zeros(resolution,length(velocities_IC));
total_H2O2_A=zeros(resolution,1);
total_H2O2_B=zeros(resolution,1);
    
parfor j=1:resolution
    decay=decay_listA(j)
    y1=y1_listA(j)
    y2=y2_listA(j)
    
   [total_concentration_A(j,:),total_velocity_A(j,:),total_H2O2_A(j)]=H2O2_UQ(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);

   decay=decay_listB(j)
    y1=y1_listB(j)
    y2=y2_listB(j)
    
   [total_concentration_B(j,:),total_velocity_B(j,:),total_H2O2_B(j)]=H2O2_UQ(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
end 
for i=1:dimension
    matrix_C=matrix_B;
    matrix_C(:,i)=matrix_A(:,i);
    decay_listC=decay_min+(decay_max-decay_min).*matrix_C(:,1);
    y1_listC=y1_min+(y1_max-y1_min).*matrix_C(:,2);
    y2_listC=y2_min+(y2_max-y2_min).*matrix_C(:,3);
    total_concentration_C=zeros(resolution,length(concentration_IC));
    total_velocity_C=zeros(resolution,length(velocities_IC));
    total_H2O2_C=zeros(resolution,1);
 parfor j=1:resolution
    decay=decay_listC(j)
    y1=y1_listC(j)
    y2=y2_listC(j)
    [total_concentration_C(j,:),total_velocity_C(j,:),total_H2O2_C(j)]=H2O2_UQ(tau,velocities_IC,T,dt,sigma,decay,time_to_detox,response_rate,critical_concentrations,l,type,concentration_IC,max_concentration,y1,y2,y3,beta,alcohol_met,Adj_mat,mM_drink,perturb);
 end
 %{
f0=1/resolution*(sum(total_concentration_A));
sensitivity_indicies_concentrations(i,:)=((1/resolution)*dot(total_concentration_A,total_concentration_C)-f0.^2)./((1/resolution)*dot(total_concentration_A,total_concentration_A)-f0.^2);
total_effect_indicies_concentrations(i,:)=1-(((1/resolution)*dot(total_concentration_B,total_concentration_C)-f0.^2)./((1/resolution)*dot(total_concentration_A,total_concentration_A)-f0.^2));

f0=1/resolution*(sum(total_velocity_A));
sensitivity_indicies_velocities(i,:)=((1/resolution)*dot(total_velocity_A,total_velocity_C)-f0.^2)./((1/resolution)*dot(total_velocity_A,total_velocity_A)-f0.^2);
total_effect_indicies_velocities(i,:)=1-(((1/resolution)*dot(total_velocity_B,total_velocity_C)-f0.^2)./((1/resolution)*dot(total_velocity_A,total_velocity_A)-f0.^2));

f0=1/resolution*(sum(total_H2O2_A));
sensitivity_indicies_H2O2(i,:)=((1/resolution)*dot(total_H2O2_A,total_H2O2_C)-f0.^2)./((1/resolution)*dot(total_H2O2_A,total_H2O2_A)-f0.^2);
total_effect_indicies_H2O2(i,:)=1-(((1/resolution)*dot(total_H2O2_B,total_H2O2_C)-f0.^2)./((1/resolution)*dot(total_H2O2_A,total_H2O2_A)-f0.^2));
 %}

 [sensitivity_indicies_concentrations(i,:), total_effect_indicies_concentrations(i,:)] = ...
     sobol_indices(total_concentration_A, total_concentration_B, total_concentration_C);
 [sensitivity_indicies_velocities(i,:), total_effect_indicies_velocities(i,:)] = ...
     sobol_indices(total_velocity_A, total_velocity_B, total_velocity_C);
 [sensitivity_indicies_H2O2(i,:), total_effect_indicies_H2O2(i,:)] = ...
     sobol_indices(total_H2O2_A, total_H2O2_B, total_H2O2_C);

end 


%% Display results
"row 1: k_d; row 2: y_1; row 3: y_2"
"Sensitivity Indicies H_2O_2"
disp(sensitivity_indicies_H2O2)
"Total effect indicies H_2O_2"
disp(total_effect_indicies_H2O2)
"Sensitivity Indicies velocities"
disp(sensitivity_indicies_velocities)
"Total effect indicies velocities"
disp(total_effect_indicies_velocities)
"Sensitivity Indicies concentrations"
disp(sensitivity_indicies_concentrations)
"Total effect indicies concentrations"
disp(total_effect_indicies_concentrations)

%% Figure: Net velocity changes

figure()
t=tiledlayout(2,2);
title(t,"Net velocity changes")
ax1 = nexttile;
scatter3(decay_listA,y2_listA,y1_listA,40,total_velocity_A(:,1),'filled')
xlabel(ax1,"k_d")
ylabel(ax1,"y_2")
zlabel(ax1,"y_1")
title(ax1,"GSH cycle")
colorbar()
clim([-7000,7000])
ax2=nexttile;
scatter3(decay_listA,y2_listA,y1_listA,40,total_velocity_A(:,2),'filled')
xlabel(ax2,"k_d")
ylabel(ax2,"y_2")
zlabel(ax2,"y_1")
title(ax2, "MET cycle")
colorbar()
clim([-7000,7000])
ax3=nexttile;
scatter3(decay_listA,y2_listA,y1_listA,40,total_velocity_A(:,3),'filled')
xlabel(ax3,"k_d")
ylabel(ax3,"y_2")
zlabel(ax3,"y_1")
title(ax3,"Urea cycle")
colorbar()
clim([-7000,7000])
ax4=nexttile;
scatter3(decay_listA,y2_listA,y1_listA,40,total_velocity_A(:,4),'filled')
xlabel(ax4,"k_d")
ylabel(ax4,"y_2")
zlabel(ax4,"y_1")
title(ax4,"FABs cycle")
colorbar()
clim([-7000,7000])
linkaxes([ax1,ax3],'x');
linkaxes([ax1,ax2],'y');
linkaxes([ax2,ax4],'x');
linkaxes([ax3,ax4],'y');

%% Figure: Check on 

figure()
t=tiledlayout(2,2);
title(t,"Net velocity changes: increase (1) or decrease (0)")
ax1 = nexttile;
scatter3(decay_listA,y2_listA,y1_listA,40,total_velocity_A(:,1)>0,'filled')
xlabel(ax1,"k_d")
ylabel(ax1,"y_2")
zlabel(ax1,"y_1")
title(ax1,"GSH cycle")
colorbar()
clim([0,1])
ax2=nexttile;
scatter3(decay_listA,y2_listA,y1_listA,40,total_velocity_A(:,2)>0,'filled')
xlabel(ax2,"k_d")
ylabel(ax2,"y_2")
zlabel(ax2,"y_1")
title(ax2, "MET cycle")
colorbar()
clim([0,1])
ax3=nexttile;
scatter3(decay_listA,y2_listA,y1_listA,40,total_velocity_A(:,3)>0,'filled')
xlabel(ax3,"k_d")
ylabel(ax3,"y_2")
zlabel(ax3,"y_1")
title(ax3,"Urea cycle")
colorbar()
clim([0,1])
ax4=nexttile;
scatter3(decay_listA,y2_listA,y1_listA,40,total_velocity_A(:,4)>0,'filled')
xlabel(ax4,"k_d")
ylabel(ax4,"y_2")
zlabel(ax4,"y_1")
title(ax4,"FABs cycle")
colorbar()
clim([0,1])
linkaxes([ax1,ax3],'x');
linkaxes([ax1,ax2],'y');
linkaxes([ax2,ax4],'x');
linkaxes([ax3,ax4],'y');

%% ROS
figure()
title("Total ROS")
scatter3(decay_listA,y2_listA,y1_listA,40,total_H2O2_A(:,1)>1000,'filled')
xlabel("k_d")
ylabel("y_2")
zlabel("y_1")
colorbar()
%% functions
function [S, ST] = sobol_indices(fA, fB, fC)
% fA, fB, fC: N x k output matrices (one column per output)
% C = B with column i from A, so fA,fC share x_i; fB,fC differ only in x_i
V  = var([fA; fB], 1);                    % pooled variance, 1 x k
S  = mean(fA .* (fC - fB)) ./ V;          % Saltelli 2010 first order
ST = 0.5 * mean((fB - fC).^2) ./ V;       % Jansen total effect
end