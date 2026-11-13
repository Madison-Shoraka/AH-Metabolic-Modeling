A0=load("results_normal_gridSearch.mat");
A1=load("results_example1.mat");
A2=load("results_example2.mat");
A3=load("results_example3.mat");

figure
t=tiledlayout(2,2);
title(t,"Net Velocity Changes")
ax1=nexttile;
hold on
title(ax1,"normal structure")
group=repmat(["GSH","MET","Urea","FABs"],length(A0.total_velocities_case1),1);
group=categorical(reshape(group,length(A0.total_velocities_case1)*4,1));
data=reshape(A0.total_velocities_case1,length(A0.total_velocities_case1)*4,1);
c = zeros(length(data), 3); 
c(data > 0, :) = repmat([.25 0 .5], sum(data > 0), 1);
s0 = scatter(ax1,group, data, 50,c,'filled', 'jitter', 'on', 'JitterAmount', 0.15);
s0.MarkerFaceAlpha = 0.5; % Set transparency
yline(ax1,0, 'r'); % Horizontal line at y=0
ylabel(ax1,"net velocity change (mM/min)")
xlabel(ax1,"module")
hold off;

ax2=nexttile;
hold on
title(ax2,"partial decoupling")
group=repmat(["GSH","MET","Urea","FABs"],length(A1.total_velocities_case1),1);
group=categorical(reshape(group,length(A1.total_velocities_case1)*4,1));
data=reshape(A1.total_velocities_case1,length(A1.total_velocities_case1)*4,1);
c = zeros(length(data), 3); 
c(data > 0, :) = repmat([.25 0 .5], sum(data > 0), 1);
s1 = scatter(ax2,group, data, 50,c,'filled', 'jitter', 'on', 'JitterAmount', 0.15);
s1.MarkerFaceAlpha = 0.5; % Set transparency
yline(ax2,0, 'r'); % Horizontal line at y=0
ylabel(ax2,"net velocity change (mM/min)")
xlabel(ax2,"module")
hold off;

ax3=nexttile;
hold on
title(ax3,"full decoupling")
group=repmat(["GSH","MET","Urea","FABs"],length(A2.total_velocities_case1),1);
group=categorical(reshape(group,length(A2.total_velocities_case1)*4,1));
data=reshape(A2.total_velocities_case1,length(A2.total_velocities_case1)*4,1);
c = zeros(length(data), 3); 
c(data > 0, :) = repmat([.25 0 .5], sum(data > 0), 1);
s2 = scatter(ax3,group, data, 50,c,'filled', 'jitter', 'on', 'JitterAmount', 0.15);
s2.MarkerFaceAlpha = 0.5; % Set transparency
yline(ax3,0, 'r'); % Horizontal line at y=0
ylabel(ax3,"net velocity change (mM/min)")
xlabel(ax3,"module")
hold off;

ax4=nexttile;
hold on
title(ax4,"one-way chain")
group=repmat(["GSH","MET","Urea","FABs"],length(A3.total_velocities_case1),1);
group=categorical(reshape(group,length(A3.total_velocities_case1)*4,1));
data=reshape(A3.total_velocities_case1,length(A3.total_velocities_case1)*4,1);
c = zeros(length(data), 3); 
c(data > 0, :) = repmat([.25 0 .5], sum(data > 0), 1);
s3 = scatter(ax4,group, data, 50,c,'filled', 'jitter', 'on', 'JitterAmount', 0.15);
s3.MarkerFaceAlpha = 0.5; % Set transparency
yline(ax4,0, 'r'); % Horizontal line at y=0
ylabel(ax4,"net velocity change (mM/min)")
xlabel(ax4,"module")
hold off;


disp("**Percent failure**")
disp("-- normal structure--")
disp(sum(A0.total_H2O2_case1>0)/length(A0.total_H2O2_case1))
disp("-- partially decoupled structure--")
disp(sum(A1.total_H2O2_case1>0)/length(A1.total_H2O2_case1))
disp("-- fully decoupled structure--")
disp(sum(A2.total_H2O2_case1>0)/length(A2.total_H2O2_case1))
disp("-- one-way chain structure--")
disp(sum(A3.total_H2O2_case1>0)/length(A3.total_H2O2_case1))


disp("**Percent MET downregulation**")
disp("-- normal structure--")
disp(sum(A0.total_velocities_case1(:,2)<0)/length(A0.total_velocities_case1(:,2)))
disp("-- partially decoupled structure--")
disp(sum(A1.total_velocities_case1(:,2)<0)/length(A1.total_velocities_case1(:,2)))
disp("-- fully decoupled structure--")
disp(sum(A2.total_velocities_case1(:,2)<0)/length(A2.total_velocities_case1(:,2)))
disp("-- one-way chain structure--")
disp(sum(A3.total_velocities_case1(:,2)<0)/length(A3.total_velocities_case1(:,2)))

disp("**Percent Urea downregulation**")
disp("-- normal structure--")
disp(sum(A0.total_velocities_case1(:,3)<0)/length(A0.total_velocities_case1(:,3)))
disp("-- partially decoupled structure--")
disp(sum(A1.total_velocities_case1(:,3)<0)/length(A1.total_velocities_case1(:,3)))
disp("-- fully decoupled structure--")
disp(sum(A2.total_velocities_case1(:,3)<0)/length(A2.total_velocities_case1(:,3)))
disp("-- one-way chain structure--")
disp(sum(A3.total_velocities_case1(:,3)<0)/length(A3.total_velocities_case1(:,3)))

disp("**Percent FABs upregulation**")
disp("-- normal structure--")
disp(sum(A0.total_velocities_case1(:,4)>0)/length(A0.total_velocities_case1(:,4)))
disp("-- partially decoupled structure--")
disp(sum(A1.total_velocities_case1(:,4)>0)/length(A1.total_velocities_case1(:,4)))
disp("-- fully decoupled structure--")
disp(sum(A2.total_velocities_case1(:,4)>0)/length(A2.total_velocities_case1(:,4)))
disp("-- one-way chain structure--")
disp(sum(A3.total_velocities_case1(:,4)>0)/length(A3.total_velocities_case1(:,4)))